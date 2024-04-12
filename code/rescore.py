"""
Thread tune combo pipeline had an issue: dG were with constraints!
So the interface was not scored correctly.
"""

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
                                get_sequence, extract_coords, read_jsonl,
                              constrain_chainbreak, freeze_atom, relax_chainA, init_pyrosetta)

work_path = Path(os.environ.get('WORKPATH', 'output_redux2'))
pdb_path = work_path / 'tuned_pdbs'
out_path = work_path / 'complex_pdbs'
os.makedirs(out_path, exist_ok=True)
filenames = []
rescore_logfile = work_path / 'rescored.jsonl'

# --------------------------------
print('\n## Init PyRosetta\n')
init_pyrosetta()
vanilla_scorefxn: pr_scoring.ScoreFunction = pyrosetta.get_fa_scorefxn()
ref: pyrosetta.Pose = pyrosetta.pose_from_file('pentakaihemimer.relax.pdb')
streptavidins = extract_streptavidins(ref, cardinality=2)
others: pyrosetta.Pose = extract_not_chainA( ref )
# --------------------------------
# Define functions

def interfaceB_relax(pose: pyrosetta.Pose, cycles=1):
    scorefxn: pr_scoring.ScoreFunction = pyrosetta.get_fa_scorefxn()
    rs: ModuleType = prc.select.residue_selector
    chainA_sele: rs.ResidueSelectorr = rs.ChainSelector('A')
    neigh_sele: rs.ResidueSelector = prc.select.residue_selector.NeighborhoodResidueSelector(chainA_sele, False, 5)
    neighs: pru.vector1_bool = neigh_sele.apply(pose)
    relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, cycles)
    movemap = pyrosetta.MoveMap()
    movemap.set_bb(False)
    movemap.set_chi(neighs)
    movemap.set_jump(False)
    relax.set_movemap(movemap)
    relax.apply(pose)

def run_process(target_name):
    info = dict(name=target_name, start=time.time(), status='error')
    # ## skip done
    if (out_path / f'{target_name}.pdb').exists():
        print(f'{target_name} done already')
    # ## Load pose
    monomer: pyrosetta.Pose = pyrosetta.pose_from_file(f'{pdb_path}/{target_name}.pdb').split_by_chain(1)
    oligomer: pyrosetta.Pose = monomer.clone()
    add_chain(oligomer, others)  # as this is in place
    info['sequence'] = monomer.sequence()   # repetition for sanity checks
    info['complex_dG'] = vanilla_scorefxn(oligomer)
    info['monomer_dG'] = vanilla_scorefxn(monomer)
    # combine with streptavidins only
    minicomplex = monomer.clone()
    minicomplex.remove_constraints()
    add_chain(minicomplex, streptavidins)
    minicomplex = monomer.clone()
    add_chain(minicomplex, streptavidins)
    interfaceB_relax(minicomplex)
    minicomplex.dump_pdb(f'{out_path}/{target_name}.pdb')
    info.update(**score_interface(minicomplex, 'A_BC'))
    res_scores = []
    scorefxn: pr_scoring.ScoreFunction = pyrosetta.get_fa_scorefxn()
    scorefxn(minicomplex)
    for i in range(1, monomer.total_residue() + 1):
        v = pru.vector1_bool(minicomplex.total_residue())
        v[i] = 1
        # score only monomer residues, but on oligomer
        res_scores.append(scorefxn.get_sub_score(minicomplex, v))
    info['per_res_scores'] = res_scores
    info['max_per_res_score'] = max(res_scores)
    pose2pdb = oligomer.pdb_info().pose2pdb
    chainA_sele = prc.select.residue_selector.ChainSelector('A')
    for distance in (1, 2, 3, 4, 5, 6, 8, 10, 12):
        v = prc.select.residue_selector.NeighborhoodResidueSelector(chainA_sele, distance, False).apply(oligomer)
        close_residues = prc.select.get_residues_from_subset(v)
        info[f'N_close_residues_{distance}'] = len(close_residues)
        info[f'close_residues_{distance}'] = [pose2pdb(r) for r in close_residues]
    info['end'] = time.time()
    with Path(rescore_logfile).open('a') as f:
        f.write(json.dumps(info) + '\n')
    return info

# --------------------------------
print('\n## Running\n')

futures: List[ProcessFuture] = []
preresults: List[dict] = []
results: List[dict] = []
timeout = 60 * 60 * 1
# this is the 5 and a half mer without A. It is therefore a 4 & 1/2 mer
ref: pyrosetta.Pose = pyrosetta.pose_from_file('pentakaihemimer.relax.pdb')

with ProcessPool(max_workers=os.cpu_count() - 1, max_tasks=0) as pool:
    # queue jobs
    target_names = json.loads(Path('to_rescore.json').read_text())
    random.shuffle(target_names)
    for target_name in target_names:
        future: ProcessFuture = pool.schedule(run_process,
                                              kwargs=dict(target_name=target_name),
                                              timeout=timeout)
        futures.append(future)
    print(f'Submitted {len(futures)} processes')
    # get results
    for future in futures:
        try:
            result = future.result()  # blocks until results are ready
            print(result['name'], result['status'])
        except Exception as error:
            error_msg = str(error)
            if isinstance(error, TimeoutError):
                print(f'Function took longer than {timeout} seconds {error}')
            else:
                print(f"Function raised {error}")
                traceback.print_tb(error.__traceback__)  # traceback of the function


