# -*- coding: utf-8 -*-
import os
import operator
import string
import random
import numpy as np
from Bio import SeqIO
from pathlib import Path
import re, os, traceback
from pebble import ProcessPool, ProcessFuture
from concurrent.futures import TimeoutError
import time
import json
from typing import Optional, Dict, List
from types import ModuleType
import pandas as pd
import pyrosetta
import pyrosetta_help as ph
prc: ModuleType = pyrosetta.rosetta.core
prp: ModuleType = pyrosetta.rosetta.protocols
pru: ModuleType = pyrosetta.rosetta.utility  # noqa
prn: ModuleType = pyrosetta.rosetta.numeric
prs: ModuleType = pyrosetta.rosetta.std  # noqa
pr_conf: ModuleType = pyrosetta.rosetta.core.conformation
pr_scoring: ModuleType = pyrosetta.rosetta.core.scoring
pr_options: ModuleType = pyrosetta.rosetta.basic.options
from common_pyrosetta import (add_chain, superpose, thread,get_ATOM_only,design_different, extract_streptavidins,
                                extract_not_chainA, superpose_pose_by_chain, score_interface, three_to_one,
                                get_sequence, extract_coords,
                              constrain_chainbreak, freeze_atom, relax_chainA, init_pyrosetta)

work_path = Path(os.environ.get('WORKPATH', 'output_redux2'))
template_folder = work_path / 'superposed_skeletons'
# annoying mistake:
seqs_folder = work_path / 'seqs/seqs' if (work_path / 'seqs/seqs').exists() else work_path / 'seqs'
seq_paths = list(seqs_folder.glob('*.fa'))
random.shuffle(seq_paths)
raw_out_path = work_path / 'unrelaxed_pdbs'
relaxed_out_path = work_path / 'relaxed_pdbs'
tuned_out_path = work_path / 'tuned_pdbs'
os.makedirs(raw_out_path, exist_ok=True)
os.makedirs(relaxed_out_path, exist_ok=True)
os.makedirs(tuned_out_path, exist_ok=True)


design_cycles = 15
relax_cycles = 5

# --------------------------------
print('\n## Init PyRosetta\n')
init_pyrosetta()

clasher = pyrosetta.pose_from_file('clasher.pdb')
xyz_clasher = extract_coords(clasher)
del clasher
strep = pyrosetta.pose_from_file('pentakaihemimer.relax.pdb').split_by_chain(2)
xyz_strep = extract_coords(strep)
del strep
# --------------------------------
# Define functions


def run_process(**inputs):
    target_name = inputs['target_name']
    info = inputs.copy()
    info['name'] = target_name
    info['start'] = time.time()
    (raw_out_path / f'{target_name}.dummy').write_text('hello world')
    # precheck against clashes
    template = pyrosetta.Pose()
    prc.import_pose.pose_from_pdbstring(template, info['template_block'])
    xyz_model = extract_coords(template)
    del template
    all_distances = np.sqrt(np.sum((xyz_clasher[:, np.newaxis, :] - xyz_model[np.newaxis, :, :]) ** 2, axis=-1))
    n_clashing = np.count_nonzero(all_distances < 1.)
    strep_distances = np.sqrt(np.sum((xyz_strep[:, np.newaxis, :] - xyz_model[np.newaxis, :, :]) ** 2, axis=-1))
    n_strep = np.count_nonzero(strep_distances < 3.)
    info['n_clashing'] = n_clashing
    info['n_strep'] = n_strep
    if n_clashing > 0:
        print(f'Clashing {target_name}', flush=True)
        info['status'] = 'clashing'
        info['end'] = time.time()
        with Path('log.jsonl').open('a') as f:
            f.write(json.dumps(info) + '\n')
        return info
    elif n_strep < 99:  # 122 for real, but 99 for glycine control
        print(f'Weak {target_name}', flush=True)
        info['status'] = 'weak_strep'
        info['end'] = time.time()
        with Path('log.jsonl').open('a') as f:
            f.write(json.dumps(info) + '\n')
        return info
    monomer: pyrosetta.Pose = thread(**inputs)
    monomer.dump_pdb(f'{raw_out_path}/{target_name}.pdb')
    oligomer: pyrosetta.Pose = monomer.clone()
    add_chain(oligomer, others)  # as this is in place
    oligomer: pyrosetta.Pose = monomer.clone()
    add_chain(oligomer, others)  # as this is in place
    info['original_sequence'] = monomer.sequence()   # repetition for sanity checks
    # constrain
    original_start = monomer.sequence().find('KDETET') + 1
    for i in range(original_start, oligomer.total_residue() + 1):
        freeze_atom(pose=oligomer, frozen_index=i, ref_index=1)
    scorefxn: pr_scoring.ScoreFunction = pyrosetta.get_fa_scorefxn()
    scorefxn.set_weight(pr_scoring.ScoreType.atom_pair_constraint, 3)
    scorefxn.set_weight(pr_scoring.ScoreType.coordinate_constraint, 1)
    vanilla_relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, relax_cycles)
    movemap = pyrosetta.MoveMap()
    # move front part starting from residue 2 and up to KDETET
    v = pru.vector1_bool(oligomer.total_residue())
    for i in range(2, original_start + 1):
        v[i] = 1
    movemap.set_bb(v)
    movemap.set_chi(v)
    movemap.set_jump(False)
    vanilla_relax.set_movemap(movemap)
    vanilla_relax.apply(oligomer)
    info['dG'] = scorefxn(oligomer)
    monomer: pyrosetta.Pose = oligomer.split_by_chain(1)
    monomer.dump_pdb(str(relaxed_out_path / f'{target_name}.pdb'))
    print(f'Relaxed {target_name}', flush=True)
    if info['dG'] > -1e3:
        info['status'] = 'score too high'
        info['end'] = time.time()
        with Path('log.jsonl').open('a') as f:
            f.write(json.dumps(info) + '\n')
        return info
    # design
    print(f'Designing {target_name}', flush=True)
    design_different(oligomer, ref, cycles=design_cycles, scorefxn=scorefxn)
    monomer = oligomer.split_by_chain(1)
    print(f'designed {target_name}', flush=True)
    monomer.dump_pdb(str(tuned_out_path / f'{target_name}.pdb'))
    info['complex_dG_post'] = scorefxn(oligomer)
    info['monomer_dG_post'] = scorefxn(monomer)
    info['tweaked_sequence'] = monomer.sequence()
    streptavidins = extract_streptavidins(ref, cardinality=2)
    superpose_pose_by_chain(oligomer, streptavidins, 'B', strict=False)
    # combine with streptavidins only
    minicomplex = monomer.clone()
    minicomplex.remove_constraints()
    add_chain(minicomplex, streptavidins)
    info.update(**score_interface(minicomplex, 'A_BC'))
    info['status'] = 'tuned'
    vanilla_relax.apply(oligomer)
    info['designed_dG'] = scorefxn(oligomer)
    pose2pdb = oligomer.pdb_info().pose2pdb
    chainA_sele = prc.select.residue_selector.ChainSelector('A')
    for distance in (1, 2, 3, 4, 5, 6, 8, 10, 12):
        v = prc.select.residue_selector.NeighborhoodResidueSelector(chainA_sele, distance, False).apply(oligomer)
        close_residues = prc.select.get_residues_from_subset(v)
        info[f'N_close_residues_{distance}'] = len(close_residues)
        info[f'close_residues_{distance}'] = [pose2pdb(r) for r in close_residues]
    res_scores = []
    for i in range(1, monomer.total_residue() + 1):
        v = pru.vector1_bool(monomer.total_residue())
        v[i] = 1
        # score only monomer residues, but on oligomer
        res_scores.append(scorefxn.get_sub_score(oligomer, v))
    info['per_res_score'] = res_scores
    info['max_per_res_scores'] = max(res_scores)
    info['end'] = time.time()
    with Path('log.jsonl').open('a') as f:
        f.write(json.dumps(info) + '\n')
    return info

# --------------------------------
print('\n## Running\n')

futures: List[ProcessFuture] = []
preresults: List[dict] = []
results: List[dict] = []
timeout = 60 * 60 * 24
# this is the 5 and a half mer without A. It is therefore a 4 & 1/2 mer
ref: pyrosetta.Pose = pyrosetta.pose_from_file('pentakaihemimer.relax.pdb')
others: pyrosetta.Pose = extract_not_chainA( ref )

with ProcessPool(max_workers=os.cpu_count() - 1, max_tasks=0) as pool:
    # queue jobs
    for path in seq_paths:
        if Path(f'{raw_out_path}/{path.stem}Ø.pdb').exists():
            print(f'P{path.stem} done already')
            continue
        pdb_path = Path(template_folder) / f'{path.stem}.pdb'
        if not pdb_path.exists():
            print(f'{path.stem} is missing ({pdb_path})', flush=True)
            continue
        for seq_record in SeqIO.parse(path, 'fasta'):
            # first record is input:
            # >beta_0-checkpoint, score=2.4184, global_score=1.5778, fixed_chains=['B', 'C', 'D'], designed_chains=['A'], model_name=v_48_020, git_hash=unknown, seed=37
            # herein I will call the input sample zero and use it as a negative control
            # others are:
            # >T=0.1, sample=1, score=1.0129, global_score=1.5088, seq_recovery=0.0312
            metadata = {k: float(v) for k, v in re.findall(r'([\w_]+)=([\d.]+)', seq_record.description)}
            metadata['name']: str = path.stem  # noqa
            metadata['sequence']: str = str(seq_record.seq) # noqa
            preresults.append(metadata)
            monomer_block = get_ATOM_only(pdb_path.read_text())
            template_name = metadata['name']
            template_seq = get_sequence(monomer_block)
            target_seq = metadata['sequence'][:len(template_seq)]  # not going to thread the other chains, only A.
            #assert len(template_seq) == len(target_seq), f'The sequences for {target_name} does not match length'
            # neg control will be Ø which looks like 0̷
            # okay. Ø er ikke det sidste bogstav. But shmeh
            target_name = f"{metadata['name']}{'ABCDEFGHIJKLMNOPQRSTUVWXYZØ'[int(metadata.get('sample', -1))]}"
            if (raw_out_path / f'{target_name}.pdb').exists() \
                or (raw_out_path / f'{target_name}.dummy').exists() \
                or (relaxed_out_path / f'{target_name}.pdb').exists():
                continue
            future: ProcessFuture = pool.schedule(run_process,
                                                  kwargs=dict(template_block=monomer_block, target_seq=target_seq,
                                               target_name=target_name, template_name=template_name),
                                                  timeout=timeout)
            futures.append(future)
    print(f'Submitted {len(futures)} processes')
    # get results
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
df.to_pickle('tuned.pkl.gz')


