# -*- coding: utf-8 -*-
"""

first record is input:

    >beta_0-checkpoint, score=2.4184, global_score=1.5778, fixed_chains=['B', 'C', 'D'], designed_chains=['A'], model_name=v_48_020, git_hash=unknown, seed=37

herein I will call the input sample zero and use it as a negative control
others are:

    >T=0.1, sample=1, score=1.0129, global_score=1.5088, seq_recovery=0.0312
"""

import gzip
import json
import os
import pickle
import random
import re
import time
import traceback
from concurrent.futures import TimeoutError
from pathlib import Path
from types import ModuleType

import numpy as np
import pandas as pd
import pyrosetta
import pyrosetta_help as ph
from Bio import SeqIO
from pebble import ProcessPool, ProcessFuture

prc: ModuleType = pyrosetta.rosetta.core
prp: ModuleType = pyrosetta.rosetta.protocols
pru: ModuleType = pyrosetta.rosetta.utility  # noqa
prn: ModuleType = pyrosetta.rosetta.numeric
prs: ModuleType = pyrosetta.rosetta.std  # noqa
pr_conf: ModuleType = pyrosetta.rosetta.core.conformation
pr_scoring: ModuleType = pyrosetta.rosetta.core.scoring
pr_options: ModuleType = pyrosetta.rosetta.basic.options
from common_pyrosetta import (thread, extract_coords, create_design_tf,
                              freeze_atom, init_pyrosetta, constrain_chainbreak,
                              fix_starts, steal_frozen,
                              appraise_itxns)
from typing import List, Any, Dict, Tuple
CTypeIdx = int  # 0-indexed
FTypeIdx = int  # 1-indexed

strep_seq = 'MEAGI'

SETTINGS = {
    'timeout': 60 * 60 * 24,  # 24 hours
    'exception_to_catch': Exception,
    'relax_cycles': 5,
    'design_cycles': 15,
    'clash_dist_cutoff': 1.5,
    'bond_dist_cutoff': 1.7, # N-C ideally is 1.32 Å
    'atom_pair_constraint_weight': 3,
    'coordinate_constraint_weight': 1,
    'initial_max_clashes': 3,  # a clash or two is fine for now
    'chain_letters': 'ACDEFGB',
    # full AHIR, strep x2, CTD AHIR, strep x2, AHIR snippet
    'start_seqs': ['MKIYY', strep_seq, strep_seq, 'GEFAR', strep_seq, strep_seq, 'FKDET']
}

# --------------------------------
# Define functions

def write_log(info, log_path=Path('log.jsonl')):
    with log_path.open('a') as f:
        f.write(json.dumps(info) + '\n')
    return info

def dump_pdbgz(pose: pyrosetta.Pose, path: Path):
    with gzip.open(path, 'wt') as f:
        f.write(ph.get_pdbstr(pose))

def get_log(target_name: str, log_path=Path('log.jsonl')):
    if not log_path.exists():
        return None
    with log_path.open('r') as f:
        logs = [json.loads(line) for line in f]
    for log in logs:
        if log['name'] == target_name:
            return log
    else:
        return None

def safe_late_stop(pose, info):
    if (time.time() - info['start']) * 0.90 > SETTINGS['timeout']:
        info['status'] = 'close to timeout'
        appraise_itxns(pose, max_clashes=0, clash_dist_cutoff=SETTINGS['clash_dist_cutoff'], bond_dist_cutoff=SETTINGS['bond_dist_cutoff'])
        raise ValueError('Close to timeout')

def run_process(target_folder: Path,
                target_name: str,
                target_sequence: str,
                hallucination_name: str,
                parent_filename: str,
                metadata: dict):
    """

    :param target_folder: folder with the RFdiffusion results
    :param target_name: name of the construct with a given sequence (eg. 'beta_0A')
    :param target_sequence: ProteinMPNN sequence
    :param hallucination_name: RFdiffusion output PDB 'hallucination'  (eg. 'beta_0')
    :param parent_filename: the template pdb used for the hallucination
    :param metadata: a dictionary of useless info
    :return:
    """
    # start info
    info = metadata.copy()
    info['folder'] = target_folder.as_posix()
    info['target_name'] = target_name
    info['target_sequence'] = target_sequence  # the seq wanted by proteinMPNN not the seq of the hallucination
    info['hallucination_name'] = hallucination_name
    info['name'] = target_name
    info['start'] = time.time()
    info['status'] = 'ongoing'
    info['is_already_done'] = False
    try:
        init_pyrosetta()
        prior = get_log(target_name, log_path=target_folder / 'log.jsonl')
        if prior:
            prior['is_already_done'] = True
            return prior
        write_log(info, log_path=target_folder / 'log.jsonl')
        # ## housekeeping
        # gz --> `dump_pdbgz`
        raw_out_path = target_folder / 'unrelaxed_pdbs' / f'{target_name}.pdb.gz'
        relaxed_out_path = target_folder / 'relaxed_pdbs' / f'{target_name}.pdb.gz'
        tuned_out_path = target_folder / 'tuned_pdbs' / f'{target_name}.pdb.gz'
        os.makedirs(raw_out_path.parent, exist_ok=True)
        os.makedirs(relaxed_out_path.parent, exist_ok=True)
        os.makedirs(tuned_out_path.parent, exist_ok=True)
        # ## read hallucination
        hallucination_pdb_path = target_folder / f'{hallucination_name}.pdb'
        hallucination: pyrosetta.Pose = pyrosetta.pose_from_file(hallucination_pdb_path.as_posix())
        # ## read metadata
        trb_path =  target_folder / f'{hallucination_name}.trb'
        assert trb_path.exists(), f'{trb_path} does not exist'
        trb: Dict[str, Any] = pickle.load(trb_path.open('rb'))
        # ## Fix chains
        parent = pyrosetta.pose_from_file(str(parent_filename))
        fix_starts(hallucination, chain_letters=SETTINGS['chain_letters'], start_seqs=SETTINGS['start_seqs'])
        rmsd, tem2hal_idx1s = steal_frozen(hallucination, parent, trb, move_acceptor=False)
        info['parent_hallucination_RMSD'] = rmsd
        info['N_conserved_parent_hallucination_atoms'] = len(tem2hal_idx1s)
        info['status'] = 'sidechain_fixed'
        n_clashing, n_warning_stretch = appraise_itxns(hallucination,
                                                       max_clashes=SETTINGS['initial_max_clashes'],
                                                       clash_dist_cutoff=SETTINGS['clash_dist_cutoff'],
                                                       bond_dist_cutoff=SETTINGS['bond_dist_cutoff'])
        info['N_clashes_start'] = n_clashing
        info['N_warning_stretches_start'] = n_warning_stretch
        info['status'] = 'initial_checks_passed'
        # ## thread
        fix_starts(hallucination, chain_letters=SETTINGS['chain_letters'], start_seqs=SETTINGS['start_seqs'])
        hal_block = ph.get_pdbstr(hallucination)
        # the seq from proteinMPNN is only chain A
        # using the full sequence slows down the threading from 16s to 1m 41s
        full_target_seq = target_sequence + hallucination.sequence()[len(hallucination.chain_sequence(1)):]
        threaded = thread(hal_block, full_target_seq, target_name, hallucination_name)
        fix_starts(threaded, chain_letters=SETTINGS['chain_letters'], start_seqs=SETTINGS['start_seqs'])
        dump_pdbgz(threaded, raw_out_path)
        info['status'] = 'threaded'
        write_log(info, log_path=target_folder / 'log.jsonl')
        # freeze the conserved atoms
        for idx0 in trb['complex_con_hal_idx0']:
            if idx0 == 0:
                continue  # pointless/troublesome constrain to self, kind of
            freeze_atom(pose=threaded, frozen_index=idx0+1, ref_index=1)  # coordinate constraint
        # enforce abolition of chainbreaks in Chain A
        for i in range(1, threaded.sequence().find('KDETET') + 1):
            constrain_chainbreak(threaded, i)
        # this is redundant with freezing all conserved atoms
        # ref_index = threaded.chain_begin(2)
        # freeze_atom(pose=threaded, frozen_index=threaded.chain_end(1), ref_index=ref_index)
        # freeze_atom(pose=threaded, frozen_index=threaded.chain_end(2), ref_index=ref_index)
        # for i in range(3, threaded.num_chains() + 1):
        #     freeze_atom(pose=threaded, frozen_index=threaded.chain_begin(i), ref_index=ref_index)
        #     freeze_atom(pose=threaded, frozen_index=threaded.chain_end(i), ref_index=ref_index)
        con_scorefxn: pr_scoring.ScoreFunction = pyrosetta.get_fa_scorefxn()
        vanilla_scorefxn: pr_scoring.ScoreFunction = pyrosetta.get_fa_scorefxn()
        con_scorefxn.set_weight(pr_scoring.ScoreType.atom_pair_constraint, SETTINGS['atom_pair_constraint_weight'])
        con_scorefxn.set_weight(pr_scoring.ScoreType.coordinate_constraint, SETTINGS['coordinate_constraint_weight'])
        vanilla_relax = pyrosetta.rosetta.protocols.relax.FastRelax(con_scorefxn, SETTINGS['relax_cycles'])
        movemap = pyrosetta.MoveMap()
        resi_sele = prc.select.residue_selector.ResidueIndexSelector()
        for idx0, b in enumerate(trb['inpaint_seq']):
            if b:
                continue
            resi_sele.append_index(idx0 + 1)
        neigh_sele = prc.select.residue_selector.NeighborhoodResidueSelector(resi_sele, 8, True)
        v: pru.vector1_bool = neigh_sele.apply(threaded)
        movemap.set_bb(v)
        movemap.set_chi(v)
        movemap.set_jump(False)
        vanilla_relax.set_movemap(movemap)
        vanilla_relax.apply(threaded)
        info['dG'] = vanilla_scorefxn(threaded)
        info['dG_monomer'] = vanilla_scorefxn(threaded.split_by_chain(1))
        info['status'] = 'relaxed'
        write_log(info, log_path=target_folder / 'log.jsonl')
        dump_pdbgz(threaded, relaxed_out_path)
        relaxed = threaded
        safe_late_stop(relaxed, info)

        # design
        ref_seq = ''
        for b, resn in zip(trb['inpaint_seq'], relaxed.sequence()):
            ref_seq += resn if b else '-'
        previous_design = relaxed
        current_design = relaxed
        previous_complex_dG = vanilla_scorefxn(previous_design)
        previous_mono_dG = vanilla_scorefxn(previous_design.split_by_chain(1))
        for cycle in range(SETTINGS['design_cycles']):
            task_factory: prc.pack.task.TaskFactory = create_design_tf(previous_design, design_sele=resi_sele, distance=0)
            relax = pyrosetta.rosetta.protocols.relax.FastRelax(vanilla_scorefxn, 1)  # one cycle at a time
            relax.set_enable_design(True)
            relax.set_task_factory(task_factory)
            current_design = previous_design.clone()
            relax.apply(current_design)
            current_complex_dG = vanilla_scorefxn(current_design)
            chain = current_design.split_by_chain(1)
            current_mono_dG = vanilla_scorefxn(chain)
            info[f'design_cycle{cycle}_seq'] = chain.sequence()
            info[f'design_cycle{cycle}_dG_complex'] = current_complex_dG
            info[f'design_cycle{cycle}_dG_monomer'] = current_mono_dG
            if any([have != expected and expected != '-' for have, expected in zip(current_design.sequence(), ref_seq)]):
                print('Mismatch happened: reverting!')  # this is a weird glitch in Relax
                current_design = previous_design
                info[f'design_cycle{cycle}_outcome'] = 'mismatched'
            elif current_complex_dG > previous_complex_dG:
                print('Design is worse: reverting!')
                current_design = previous_design
                info[f'design_cycle{cycle}_outcome'] = 'worse complex'
            elif current_mono_dG > previous_mono_dG:
                print('Design is worse: reverting!')
                current_design = previous_design
                info[f'design_cycle{cycle}_outcome'] = 'worse monomer'
            else:
                info[f'design_cycle{cycle}_outcome'] = 'success'
                previous_design = current_design
                previous_complex_dG = current_complex_dG
                dump_pdbgz(current_design, tuned_out_path)
                write_log(info, log_path=target_folder / 'log.jsonl')
                safe_late_stop(current_design, info)
        designed = current_design
        n_clashing, n_warning_stretch = appraise_itxns(designed,
                                                       max_clashes=0,
                                                       clash_dist_cutoff=SETTINGS['clash_dist_cutoff'],
                                                       bond_dist_cutoff=SETTINGS['bond_dist_cutoff'])
        info['N_clashes_end'] = n_clashing  # 0. pro forma
        info['N_warning_stretches_end'] = n_warning_stretch
        chainA_sele = prc.select.residue_selector.ChainSelector('A')
        pose2pdb = designed.pdb_info().pose2pdb
        for distance in (1, 2, 3, 4, 5, 6, 8, 10, 12):
            v = prc.select.residue_selector.NeighborhoodResidueSelector(chainA_sele, distance, False).apply(designed)
            close_residues = prc.select.get_residues_from_subset(v)
            info[f'N_close_residues_{distance}'] = len(close_residues)
            info[f'close_residues_{distance}'] = [pose2pdb(r) for r in close_residues]
        res_scores = []
        monomer = designed.split_by_chain(1)
        for i in range(1, monomer.total_residue() + 1):
            v = pru.vector1_bool(designed.total_residue())
            v[i] = 1
            # score only monomer residues, but on oligomer
            res_scores.append(vanilla_scorefxn.get_sub_score(designed, v))
        info['per_res_score'] = res_scores
        info['max_per_res_scores'] = max(res_scores)

        # hurray:
        info['status'] = 'complete'
    except SETTINGS['exception_to_catch'] as e:
        info['error_type'] = e.__class__.__name__
        info['error'] = str(e)
        info['traceback'] = traceback.format_exception(e)
        info['status'] = 'error'
    info['end'] = time.time()
    write_log(info, log_path=target_folder / 'log.jsonl')
    return info

def get_novels(target_folder, log_path):
    """
    Get novel sequences from the seqs folder that are not marked in log_threaded.txt

    :param target_folder:
    :return:
    """
    seqs_folder = target_folder / 'seqs'
    log_block = log_path.read_text() if log_path.exists() else ''
    seq_paths = []
    for path in seqs_folder.glob('*.fa'):
        if path.stem in log_block:  # it is mentioned?
            continue
        # ## PDB check
        # does it have a PDB?
        pdb_path = target_folder / f'{path.stem}.pdb'
        if not pdb_path.exists():
            print(f'{path.stem} is missing ({pdb_path})', flush=True)
            continue
        # is it done already? yet not logged?!
        if (target_folder / 'unrelaxed_pdbs' / f'{path.stem}Ø.pdb.gz').exists():
            print(f'{path.stem} done already')
            continue
        seq_paths.append(path)
    random.shuffle(seq_paths)
    return seq_paths

def get_max_cores():
    """
    the number of cores to use.
    Called by main
    """
    if os.environ['SLURM_JOB_CPUS_PER_NODE']:
        return int(os.environ['SLURM_JOB_CPUS_PER_NODE'])
    else:
        return os.cpu_count()

# --------------------------------


def main(target_folder: Path,
         parent_filename: Path,
         timeout = SETTINGS['timeout']):
    seq_paths = get_novels(target_folder, log_path=Path(target_folder) / 'log.jsonl')
    futures: List[ProcessFuture] = []
    future_names: List[str] = []
    with ProcessPool(max_workers=get_max_cores() - 1, max_tasks=0) as pool:
        # queue jobs
        for path in seq_paths:
            # send each seq to the pool
            for seq_record in SeqIO.parse(path, 'fasta'):
                # ## Read metadata
                metadata = {k: float(v) for k, v in re.findall(r'([\w_]+)=([\d.]+)', seq_record.description)}
                hallucination_name: str = path.stem  # noqa
                target_sequence: str = str(seq_record.seq) # noqa
                target_name = f"{hallucination_name}{'ABCDEFGHIJKLMNOPQRSTUVWXYZØ'[int(metadata.get('sample', -1))]}"
                future: ProcessFuture = pool.schedule(run_process,
                                                      # target = seq-variant of hallucination
                                                      # parent = WT template
                                                      # hallucination = RFdiffused skeleton
                                                      kwargs=dict(target_folder=target_folder,
                                                                  target_name=target_name,
                                                                  target_sequence=target_sequence,
                                                                  hallucination_name=hallucination_name,
                                                                  parent_filename=parent_filename,
                                                                  metadata=metadata,
                                                                  ),
                                                      timeout=timeout)

                future_names.append(target_name)
                futures.append(future)
        print(f'Submitted {len(futures)} processes')
        # ## Get results
        results: List[dict] = []
        for name, future in zip(future_names, futures):
            try:
                result = future.result()  # blocks until results are ready
                print(result['name'], result['status'])
                results.append(result)
            except Exception as error:
                error_msg = str(error)
                result = dict(target_name=name, error=str(error), error_type=error.__class__.__name__,)
                results.append(result)
                if isinstance(error, TimeoutError):
                    print(f'Function took longer than {timeout} seconds {error}')
                else:
                    print(f"Function raised {error}")
                    traceback.print_tb(error.__traceback__)  # traceback of the function
    df = pd.DataFrame(results)
    df.to_pickle(target_folder / 'tuned.pkl.gz')
    return df


# monomer_block = get_ATOM_only(pdb_path.read_text())
#
# # assert len(template_seq) == len(target_seq), f'The sequences for {target_name} does not match length'
# # neg control will be Ø which looks like 0̷. Ø er ikke det sidste bogstav. But shmeh
# if (raw_out_path / f'{target_name}.pdb').exists() \
#         or (raw_out_path / f'{target_name}.dummy').exists() \
#         or (relaxed_out_path / f'{target_name}.pdb').exists():
#     continue

# --------------------------------

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Thread and tune a sequence onto a template')
    parser.add_argument('target_folder', type=str, help='A folder with seqs to thread')
    parser.add_argument('parent_pdb', type=str, help='A PDB file with the template')
    parser.add_argument('--timeout', type=int, default=SETTINGS['timeout'], help='Timeout for each process')
    parser.add_argument('--relax_cycles', type=int, default=SETTINGS['relax_cycles'], help='Number of relax cycles')
    parser.add_argument('--design_cycles', type=int, default=SETTINGS['design_cycles'], help='Number of design cycles')

    args = parser.parse_args()
    SETTINGS['timeout'] = int(args.timeout)
    SETTINGS['relax_cycles'] = int(args.relax_cycles)
    SETTINGS['design_cycles'] = int(args.design_cycles)
    main(target_folder=Path(args.target_folder),
         parent_filename=Path(args.parent_pdb))
    print('Done')


