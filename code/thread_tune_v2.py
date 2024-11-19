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
                              freeze_atom, init_pyrosetta)
from typing import List, Any, Dict, Tuple
CTypeIdx = int  # 0-indexed
FTypeIdx = int  # 1-indexed

strep_seq = 'MEAGI'

SETTINGS = {
    'exception_to_catch': Exception,
    'relax_cycles': 5,
    'design_cycles': 15,
    'clash_dist_cutoff': 1.5,
    'bond_dist_cutoff': 1.5, # N-C
    'atom_pair_constraint_weight': 3,
    'coordinate_constraint_weight': 1,
    'chain_letters': 'ACDEFGB',
    # full AHIR, strep x2, CTD AHIR, strep x2, AHIR snippet
    'start_seqs': ['MKIYY', strep_seq, strep_seq, 'GEFAR', strep_seq, strep_seq, 'FKDET']
}

# --------------------------------
# Define functions

def fix_starts(pose, chain_letters: str, start_seqs: List[str]):
    """
    Fix the chains
    In anything based on pentakaimer it is

    .. code-block:: python

        strep_seq = 'MEAGIT'
        start_seqs = ['MKIYY', strep_seq, strep_seq, 'GEFAR', strep_seq, strep_seq, 'FKDET']
        fix_starts(pose, chain_letters='ACDEFGB', start_seq=start_seq)

    :param pose:
    :param chain_letters:
    :param start_seq: Confusingly, the first is ignored: the start of the pose is the start of the first chain.
    :return:
    """
    pi = pose.pdb_info()
    seq = pose.sequence()
    seq_iter = iter(start_seqs[1:]+[None])
    chain_iter = iter(chain_letters)
    start_idx = 1
    while True:
        this_chain = next(chain_iter)
        next_seq = next(seq_iter)
        if next_seq is None:
            for i in range(start_idx, len(seq)+1):
                pi.chain(i, this_chain)
            break
        else:
            next_start = seq.find(next_seq, start_idx) + 1
            for i in range(start_idx, next_start):
                pi.chain(i, this_chain)
            start_idx = next_start
    pose.update_pose_chains_from_pdb_chains()


# whole alignment


def steal_frozen(acceptor: pyrosetta.Pose,
                 donor: pyrosetta.Pose, trb: Dict[str, Any],
                 move_acceptor: bool = False
                 ):
    """
    Copy all the conserved coordinates from the donor to the acceptor.
    These are determined by `trb` dict from RFdiffusion.

    The hallucination is the acceptor, the template is the donor.
    The RFdiffused pose is skeleton, but when imported the sidechains are added.
    The theft is done in 3 steps.

    1. A mapping of residue idx, atom idx to residue idx, atom idx is made.
    2. The hallucination is superimposed on the template if move_acceptor is True, else vice versa.
    3. The coordinates are copied from the template to the hallucination.

    The term 'ref' gets confusing.
     hallucination is fixed / mutanda, template is mobile
    fixed is called ref, but ref means template for RFdiffusion so it is flipped)


    :param acceptor:
    :param donor:
    :param trb:
    :return:
    """

    # ## Make mapping of all atoms of conserved residues
    donor2acceptor_idx1s: Dict[Tuple[FTypeIdx, FTypeIdx], Tuple[FTypeIdx, FTypeIdx, str]] = {}
    # these run off 0-based indices
    for donor_res_idx0, acceptor_res_idx0 in zip(trb['complex_con_ref_idx0'], trb['complex_con_hal_idx0']):
        donor_res_idx1 = donor_res_idx0 + 1
        acceptor_res_idx1 = acceptor_res_idx0 + 1
        acceptor_res = acceptor.residue(acceptor_res_idx1)
        donor_res = donor.residue(donor_res_idx1)
        assert donor_res.name3() == acceptor_res.name3(), f'donor {donor_res.name3()} != acceptor {acceptor_res.name3()}'
        mob_atomnames = [donor_res.atom_name(ai1) for ai1 in range(1, donor_res.natoms() + 1)]
        for fixed_atm_idx1 in range(1, acceptor_res.natoms() + 1):
            aname = acceptor_res.atom_name(fixed_atm_idx1)  # key to map one to other: overkill bar for HIE/HID
            if aname not in mob_atomnames:
                print(f'Template residue {donor_res.annotated_name()}{donor_res_idx1} lacks atom {aname}')
                continue
            donor_atm_idx1 = donor_res.atom_index(aname)
            donor2acceptor_idx1s[(donor_res_idx1, donor_atm_idx1)] = (acceptor_res_idx1, fixed_atm_idx1, aname)

    # ## Align
    atom_map = prs.map_core_id_AtomID_core_id_AtomID()
    if move_acceptor:
        mobile: pyrosetta.Pose = acceptor
        fixed: pyrosetta.Pose = donor
    else:
        mobile: pyrosetta.Pose = donor
        fixed: pyrosetta.Pose = acceptor
    for (donor_res_idx1, donor_atm_idx1), (acceptor_res_idx1, acceptor_atm_idx1, aname) in donor2acceptor_idx1s.items():
        if move_acceptor:
            mob_res_idx1, mob_atm_idx1 = acceptor_res_idx1, acceptor_atm_idx1
            fixed_res_idx1, fixed_atm_idx1 = donor_res_idx1, donor_atm_idx1
        else:
            mob_res_idx1, mob_atm_idx1 = donor_res_idx1, donor_atm_idx1
            fixed_res_idx1, fixed_atm_idx1 = acceptor_res_idx1, acceptor_atm_idx1
        if aname.strip() not in ('N', 'CA', 'C', 'O'):  # BB alignment
            continue
        fixed_atom = pyrosetta.AtomID(fixed_atm_idx1, fixed_res_idx1)
        mobile_atom = pyrosetta.AtomID(mob_atm_idx1, mob_res_idx1)
        atom_map[mobile_atom] = fixed_atom
    rmsd = prc.scoring.superimpose_pose(mod_pose=mobile, ref_pose=fixed, atom_map=atom_map)
    # I am unsure why this is not near zero but around 0.1–0.3
    assert rmsd < 1, f'RMSD {rmsd} is too high'

    # ## Copy coordinates
    to_move_atomIDs = pru.vector1_core_id_AtomID()
    to_move_to_xyz = pru.vector1_numeric_xyzVector_double_t()
    for (donor_res_idx1, donor_atm_idx1), (acceptor_res_idx1, acceptor_atm_idx1, aname) in donor2acceptor_idx1s.items():
        # if aname in ('N', 'CA', 'C', 'O'):  # BB if common
        #     continue
        # this does not stick: fixed_res.set_xyz( fixed_ai1, mob_res.xyz(mob_ai1) )
        to_move_atomIDs.append(pyrosetta.AtomID(acceptor_atm_idx1, acceptor_res_idx1))
        to_move_to_xyz.append(donor.residue(donor_res_idx1).xyz(donor_atm_idx1))

    acceptor.batch_set_xyz(to_move_atomIDs, to_move_to_xyz)
    return rmsd, donor2acceptor_idx1s

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
        logs = json.load(f)
    for log in logs:
        if log['name'] == target_name:
            return log
    else:
        return None

def appraise_itxns(pose):
    """
    Assumes the chains have been fixed already.

    :param pose:
    :return:
    """
    n_clashing = 0
    chains = pose.split_by_chain()
    xyz_model = extract_coords(chains[1])
    for i in range(1, pose.num_chains()):  # chain 0 is the full AHIR, designed
        xyz_other = extract_coords(chains[i + 1])
        distances = np.sqrt(np.sum((xyz_other[:, np.newaxis, :] - xyz_model[np.newaxis, :, :]) ** 2, axis=-1))
        # 1.5 Å is too close
        n_clashing += np.count_nonzero(distances < SETTINGS['clash_dist_cutoff'])
    if n_clashing > 0:
        raise ValueError(f'Clashes')
    # check no stretch
    for chain in chains:
        for i in range(chain.total_residue() - 1):
            d: float = chain.residue(i + 1).xyz('C').distance(chain.residue(i + 1 + 1).xyz('N'))
            if d > SETTINGS['bond_dist_cutoff']:
                raise ValueError(f'Stretch {d}')

def run_process(target_folder: Path,
                target_name: str,
                target_sequence: str,
                template_name: str,
                template_filename: str,
                metadata: dict):
    # start info
    info = metadata.copy()
    info['folder'] = target_folder.as_posix()
    info['target_name'] = target_name
    info['target_sequence'] = target_sequence
    info['template_name'] = template_name
    info['name'] = target_name
    info['start'] = time.time()
    info['status'] = 'ongoing'
    info['is_already_done'] = False
    try:
        prior = get_log(target_name, log_path=target_folder / 'log.jsonl')
        if prior:
            prior['is_already_done'] = True
            return prior
        write_log(info, log_path=target_folder / 'log.jsonl')
        # housekeeping
        # gz --> `dump_pdbgz`
        raw_out_path = target_folder / 'unrelaxed_pdbs' / f'{target_name}.pdb.gz'
        relaxed_out_path = target_folder / 'relaxed_pdbs' / f'{target_name}.pdb.gz'
        tuned_out_path = target_folder / 'tuned_pdbs' / f'{target_name}.pdb.gz'
        os.makedirs(raw_out_path.parent, exist_ok=True)
        os.makedirs(relaxed_out_path.parent, exist_ok=True)
        os.makedirs(tuned_out_path.parent, exist_ok=True)
        # read hallucination
        hallucination_pdb_path = target_folder / f'{template_name}.pdb'
        hallucination: pyrosetta.Pose = pyrosetta.pose_from_file(hallucination_pdb_path.as_posix())
        # read metadata
        trb_path =  target_folder / f'{template_name}.trb'
        assert trb_path.exists(), f'{trb_path} does not exist'
        trb: Dict[str, Any] = pickle.load(trb_path.open('rb'))
        # Fix chains
        fix_starts(hallucination, chain_letters=SETTINGS['chain_letters'], start_seqs=SETTINGS['start_seqs'])
        appraise_itxns(hallucination)
        info['status'] = 'initial_checks_passed'
        # thread
        template = pyrosetta.pose_from_file(template_filename)
        fix_starts(hallucination, chain_letters=SETTINGS['chain_letters'], start_seqs=SETTINGS['start_seqs'])
        rmsd, tem2hal_idx1s = steal_frozen(hallucination, template, trb, move_acceptor=False)
        info['template_hallucination_RMSD'] = rmsd
        info['N_conserved_template_hallucination_atoms'] = len(tem2hal_idx1s)
        info['status'] = 'sidechain_fixed'
        hal_block = ph.get_pdbstr(hallucination)
        # the seq from proteinMPNN is only chain A
        # using the full sequence slows down the threading from 16s to 1m 41s
        full_target_seq = target_sequence + hallucination.sequence()[len(hallucination.chain_sequence(1)):]
        threaded = thread(hal_block, full_target_seq, target_name, template_name)
        fix_starts(hallucination, chain_letters=SETTINGS['chain_letters'], start_seqs=SETTINGS['start_seqs'])
        dump_pdbgz(threaded, raw_out_path)
        info['status'] = 'threaded'
        write_log(info, log_path=target_folder / 'log.jsonl')
        for idx0 in trb['complex_con_hal_idx0']:
            if idx0 == 0:
                continue  # pointless/troublesome constrain to self, kind of
            freeze_atom(pose=threaded, frozen_index=idx0+1, ref_index=1)  # coordinate constraint
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
        info['status'] = 'relaxed'
        write_log(info, log_path=target_folder / 'log.jsonl')
        dump_pdbgz(threaded, relaxed_out_path)
        relaxed = threaded

        # design
        ref_seq = ''
        for b, resn in zip(trb['inpaint_seq'], relaxed.sequence()):
            ref_seq += resn if b else '-'
        previous_design = relaxed
        previous_complex_dG = vanilla_scorefxn(previous_design)
        previous_mono_dG = vanilla_scorefxn(previous_design.split_by_chain(1))
        for cycle in range(SETTINGS['design_cycles']):
            task_factory: prc.pack.task.TaskFactory = create_design_tf(relaxed, design_sele=resi_sele, distance=0)
            relax = pyrosetta.rosetta.protocols.relax.FastRelax(vanilla_relax, 1)  # one cycle at a time
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
        designed = current_design
        appraise_itxns(designed)
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
        info['status'] = 'error'
    info['end'] = time.time()
    write_log(info, log_path=target_folder / 'log.jsonl')
    return info

def get_novels(target_folder):
    """
    Get novel sequences from the seqs folder that are not marked in log_threaded.txt

    :param target_folder:
    :return:
    """
    seqs_folder = target_folder / 'seqs'
    log = Path('log_threaded.txt').read_text() if  Path('log_threaded.txt').exists() else ''
    seq_paths = []
    for path in seqs_folder.glob('*.fa'):
        if path.stem in log:
            continue
        # ## PDB check
        # does it have a PDB?
        pdb_path = target_folder / f'{path.stem}.pdb'
        if not pdb_path.exists():
            print(f'{path.stem} is missing ({pdb_path})', flush=True)
            continue
        # is it done already? yet not logged?!
        if (target_folder / 'unrelaxed_pdbs' / f'{path.stem}Ø.pdb').exists():
            print(f'P{path.stem} done already')
            with Path('log_threaded.txt').open('a') as fh:
                fh.write(f'{path.stem}\toutput exists but no data\n')
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
         template_filename: Path,
         timeout = 60 * 60 * 24):
    print('\n## Init PyRosetta\n')
    init_pyrosetta()
    seq_paths = get_novels(target_folder)
    futures: List[ProcessFuture] = []
    with ProcessPool(max_workers=get_max_cores() - 1, max_tasks=0) as pool:
        # queue jobs
        for path in seq_paths:
            # send each seq to the pool
            for seq_record in SeqIO.parse(path, 'fasta'):
                # ## Read metadata
                metadata = {k: float(v) for k, v in re.findall(r'([\w_]+)=([\d.]+)', seq_record.description)}
                template_name: str = path.stem  # noqa
                target_sequence: str = str(seq_record.seq) # noqa
                target_name = f"{template_name}{'ABCDEFGHIJKLMNOPQRSTUVWXYZØ'[int(metadata.get('sample', -1))]}"
                future: ProcessFuture = pool.schedule(run_process,
                                                      kwargs=dict(target_folder=target_folder,
                                                                  target_name=target_name,
                                                                  target_sequence=target_sequence,
                                                                  template_name=template_name,
                                                                  template_filename=template_filename,
                                                                  metadata=metadata,
                                                                  ),
                                                      timeout=timeout)
                futures.append(future)
        print(f'Submitted {len(futures)} processes')
        # ## Get results
        results: List[dict] = []
        for future in futures:
            try:
                result = future.result()  # blocks until results are ready
                print(result['name'], result['status'])
                results.append(result)
            except Exception as error:
                error_msg = str(error)
                if isinstance(error, TimeoutError):
                    print(f'Function took longer than {timeout} seconds {error}')
                else:
                    print(f"Function raised {error}")
                    traceback.print_tb(error.__traceback__)  # traceback of the function
                results.append(dict(error=error.__class__.__name__, error_msg=error_msg))
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
    args = parser.parse_args()
    main(target_folder=Path(args.target_folder), template_filename=Path(args.parent_pdb))
    print('Done')


