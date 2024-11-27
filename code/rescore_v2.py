"""

Not used: constrain_chainbreak, freeze_atom
"""

import gzip
import pyrosetta
import sys
import traceback
import pickle
from Bio import SeqIO
import re, os
from pathlib import Path
from pebble import ProcessPool, ProcessFuture
from typing import List
import pandas as pd

import random
import pyrosetta_help as ph

# export PATH="$HOME2/crysalin-redux/repo/code:$PATH";
sys.path.append('/opt/xchem-fragalysis-2/mferla/crysalin-redux/repo/code')
prc = pyrosetta.rosetta.core
prp = pyrosetta.rosetta.protocols
pru = pyrosetta.rosetta.utility
from thread_tune_v2 import init_pyrosetta, fix_starts, appraise_itxns, SETTINGS, get_max_cores
from common_pyrosetta import superpose

def dump_pdbgz(pose: pyrosetta.Pose, path: Path):
    with gzip.open(path, 'wt') as f:
        f.write(ph.get_pdbstr(pose))

def read_pdbgz(path) -> pyrosetta.Pose:
    with gzip.open(path, 'r') as gfh:
        pose = pyrosetta.Pose()
        prc.import_pose.pose_from_pdbstring(pose, gfh.read())
    return pose

def make_tagless(pose):
    cloned = pose.clone()
    tag_start = cloned.sequence().find('WSHPQFEKRP') + 1
    assert tag_start != 0, 'No streptag'
    prp.grafting.delete_region(cloned, tag_start, tag_start+10)
    return cloned



def analyse_design(pdb_path: Path, folder_path: Path, template_path: Path):
    design_name = pdb_path.stem.split('.')[0]
    summary_path = (folder_path / 'summary' / f'{design_name}.pkl')
    if summary_path.exists():
        return pickle.load(summary_path.open('rb'))
    info = {'name': design_name,
            'status': 'ok',
            'is_control': design_name[-1] == 'Ø'}
    try:
        trb_path = folder_path / f'{design_name[:-1]}.trb'
        if not trb_path.exists():
            raise ValueError(f'{trb_path} does not exist')
        trb = pickle.load(trb_path.open('rb'))
        mod_idx1s = []
        for idx0, b in enumerate(trb['inpaint_seq']):
            if b:
                continue
            mod_idx1s.append(idx0 + 1)
        info['mod_idxs'] = mod_idx1s

        init_pyrosetta()
        pose = read_pdbgz(pdb_path)
        fix_starts(pose, chain_letters=SETTINGS['chain_letters'], start_seqs=SETTINGS['start_seqs'])
        info['sequence'] = pose.chain_sequence(1)
        aln_map = dict(zip(trb['complex_con_ref_idx0'], trb['complex_con_hal_idx0']))
        scorefxn = pyrosetta.get_fa_scorefxn()
        ref_pose = pyrosetta.pose_from_file(template_path.as_posix())
        info['rmsd'] = superpose(ref=ref_pose, mobile=pose, aln_map=aln_map, zero_based=True)
        tagless = pose.clone()
        make_tagless(tagless)
        tagless_monomer = tagless.split_by_chain(1)
        info['tagless_complex_dG'] = scorefxn(tagless)
        info['tagless_monomer_dG_unmin'] = scorefxn(tagless_monomer)
        n_clashing, n_warning_stretch = appraise_itxns(tagless,
                                                       max_clashes=float('inf'),
                                                       clash_dist_cutoff=float('inf'),
                                                       bond_dist_cutoff=float('inf')
                                                       )
        info['N_clashes_end'] = n_clashing  # 0. pro forma
        info['N_warning_stretches_end'] = n_warning_stretch
        chainA_sele = prc.select.residue_selector.ChainSelector('A')
        pose2pdb = tagless.pdb_info().pose2pdb
        for distance in (1, 2, 3, 4, 5, 6, 8, 10, 12):
            v = prc.select.residue_selector.NeighborhoodResidueSelector(chainA_sele, distance, False).apply(tagless)
            close_residues = prc.select.get_residues_from_subset(v)
            info[f'N_close_residues_{distance}'] = len(close_residues)
            info[f'close_residues_{distance}'] = [pose2pdb(r) for r in close_residues]
        res_scores = []
        monomer = pose.split_by_chain(1)
        for i in range(1, monomer.total_residue() + 1):
            v = pru.vector1_bool(pose.total_residue())
            v[i] = 1
            # score only monomer residues, but on oligomer
            res_scores.append(scorefxn.get_sub_score(pose, v))
        info['per_res_score'] = res_scores
        info['max_per_res_scores'] = max(res_scores)


        # make selector
        if mod_idx1s:
            selector = prc.select.residue_selector.ResidueIndexSelector()
            for idx1 in mod_idx1s:
                selector.append_index(idx1)

            hphob_sasa = prc.simple_metrics.metrics.SasaMetric(selector,
                                                               prc.scoring.sasa.SasaMethodHPMode.HYDROPHOBIC_SASA
                                                               ).calculate(tagless_monomer)
            all_sasa = prc.simple_metrics.metrics.SasaMetric(selector,
                                                             prc.scoring.sasa.SasaMethodHPMode.HYDROPHOBIC_SASA
                                                             ).calculate(tagless_monomer)
            info['span_hydrophobic_sasa'] = hphob_sasa
            info['span_all_sasa'] = all_sasa
        else:
            pass
        selector = prc.select.residue_selector.ChainSelector('A')
        hphob_sasa = prc.simple_metrics.metrics.SasaMetric(selector,
                                                           prc.scoring.sasa.SasaMethodHPMode.HYDROPHOBIC_SASA
                                                           ).calculate(tagless_monomer)
        all_sasa = prc.simple_metrics.metrics.SasaMetric(selector,
                                                         prc.scoring.sasa.SasaMethodHPMode.HYDROPHOBIC_SASA
                                                         ).calculate(tagless_monomer)
        info['monomer_hydrophobic_sasa'] = hphob_sasa
        info['monomer_all_sasa'] = all_sasa



        seq_path = folder_path / f'seqs/{design_name[:-1]}.fa'
        for letter, seq_record in zip('ØBCDEFGHIJKLMNOPQRSTUVWXYZ', SeqIO.parse(seq_path, 'fasta')):
            if letter == design_name[-1]:
                # ## Read metadata
                metadata = {k: float(v) for k, v in re.findall(r'([\w_]+)=([\d.]+)', seq_record.description)}
                info['original_seq'] = str(seq_record.seq)
                info.update(metadata)
                break
        info['trb'] = trb
        if info['tagless_complex_dG'] < -2774:
            relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, 3)
            relax.apply(tagless_monomer)
            info['tagless_monomer_dG'] = scorefxn(tagless_monomer)
            if not (folder_path / 'tuned_monomer_pdbs').exists():
                os.mkdir(folder_path / 'tuned_monomer_pdbs')
            dump_pdbgz(tagless_monomer, folder_path / f'tuned_monomer_pdbs/{design_name}.pdb.gz')
            info['status'] = 'ok'
        else:
            info['status'] = 'rubbish'
        if not summary_path.parent.exists():
            os.mkdir(summary_path.parent)
        with summary_path.open('wb') as fh:
            pickle.dump(info, fh)
    except Exception as e:
        info['error_type'] = e.__class__.__name__
        info['error'] = str(e)
        info['traceback'] = traceback.format_exception(e)
        info['status'] = 'error'
    return info

def main(target_folder: Path, template_path: Path, timeout: int):
    pdb_paths = list(target_folder.glob('tuned_pdbs/*.pdb.gz'))
    random.shuffle(pdb_paths)

    futures: List[ProcessFuture] = []
    future_names: List[str] = []
    results: List[dict] = []
    with ProcessPool(max_workers=get_max_cores() - 1, max_tasks=0) as pool:
        # queue jobs
        for pdb_path in pdb_paths:
            # send each seq to the pool
            target_name = pdb_path.name.replace('.pdb.gz', '')
            future: ProcessFuture = pool.schedule(analyse_design,
                                                      kwargs=dict(pdb_path=pdb_path,
                                                                  folder_path=target_folder,
                                                                  template_path=template_path,
                                                                  ),
                                                      timeout=timeout)

            future_names.append(target_name)
            futures.append(future)
        print(f'Submitted {len(futures)} processes')
        # ## Get results
        for name, future in zip(future_names, futures):
            try:
                result = future.result()  # blocks until results are ready
                print(result['name'], result['status'])
                results.append(result)
            except Exception as error:
                error_msg = str(error)
                result = dict(target_name=name, error=str(error), error_type=error.__class__.__name__, )
                results.append(result)
                if isinstance(error, TimeoutError):
                    print(f'Function took longer than {timeout} seconds {error}')
                else:
                    print(f"Function raised {error}")
                    traceback.print_tb(error.__traceback__)  # traceback of the function
    df = pd.DataFrame(results)
    df.to_pickle(target_folder / 'tuned_v2.pkl.gz')


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Thread and tune a sequence onto a template')
    parser.add_argument('target_folder', type=str, help='A folder with seqs to thread')
    parser.add_argument('parent_pdb', type=str, help='A PDB file with the template')
    parser.add_argument('--timeout', type=int, default=20*60, help='Timeout for each process')

    args = parser.parse_args()
    main(target_folder=Path(args.target_folder),
         template_path=Path(args.parent_pdb),
         timeout=args.timeout)
    print('Done')