from common_pyrosetta import (init_pyrosetta, design_interface_onA, superpose_pose_by_chain,
                              relax_chainA,transpant_CTD,constrain_chainbreak,freeze_atom,design_different,
                              extract_not_chainA, add_chain, extract_streptavidins, score_interface)
import multiprocessing
import time

from collections import defaultdict
from pathlib import Path
import pebble
import os
import random
import json
from typing import Optional, Dict, List, Iterator
from types import ModuleType
import warnings
import pyrosetta
prc: ModuleType = pyrosetta.rosetta.core
prp: ModuleType = pyrosetta.rosetta.protocols
pru: ModuleType = pyrosetta.rosetta.utility  # noqa
prn: ModuleType = pyrosetta.rosetta.numeric
prs: ModuleType = pyrosetta.rosetta.std  # noqa
pr_conf: ModuleType = pyrosetta.rosetta.core.conformation
pr_scoring: ModuleType = pyrosetta.rosetta.core.scoring
pr_options: ModuleType = pyrosetta.rosetta.basic.options
pr_res: ModuleType = pyrosetta.rosetta.core.select.residue_selector

timeout = 24 * 60 * 60
work_path = Path(os.environ.get('WORKPATH', 'output_redux'))
print(f'path: {work_path}')
print('starting Pyrosetta...')
init_pyrosetta()
num_cores = multiprocessing.cpu_count()

# --------------------------------
print('\n## Running\n')

futures: List[pebble.ProcessFuture] = []
preresults: List[dict] = []
results: List[dict] = []
timeout = 60 * 60 * 5
# this is the 5 and a half mer without A. It is therefore a 4 & 1/2 mer
ref: pyrosetta.Pose = pyrosetta.pose_from_file('pentakaihemimer.relax.pdb')
others: pyrosetta.Pose = extract_not_chainA( ref )
threaded_in_path = work_path / 'threaded_pdbs'
relaxed_out_path = work_path / 'relaxed_pdbs'
tuned_out_path = work_path / 'tuned_pdbs'
os.makedirs(relaxed_out_path, exist_ok=True)
os.makedirs(tuned_out_path, exist_ok=True)

def run(path: Path, cycles=5):
    # check if done already
    info = dict(name=path.stem, status='crashed', start=time.time())
    if (relaxed_out_path / path.name).exists():
        info['status'] = 'done already'
        print(info, flush=True)
        return info
    else:
        print(f'starting {path.stem}', flush=True)
        (relaxed_out_path / path.name).write_text('')
    # make complex
    monomer = pyrosetta.pose_from_file(str(path))
    info['original_sequence'] = monomer.sequence()
    oligomer: pyrosetta.Pose = monomer.clone()
    add_chain(oligomer, others)  # as this is in place
    # constrain
    original_start = monomer.sequence().find('KDETET') + 1
    for i in range(original_start, oligomer.total_residue() + 1):
        freeze_atom(pose=oligomer, frozen_index=i, ref_index=1)
    scorefxn: pr_scoring.ScoreFunction = pyrosetta.get_fa_scorefxn()
    scorefxn.set_weight(pr_scoring.ScoreType.atom_pair_constraint, 3)
    scorefxn.set_weight(pr_scoring.ScoreType.coordinate_constraint, 1)
    relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, cycles)
    movemap = pyrosetta.MoveMap()
    # move front part starting from residue 2 and up to KDETET
    v = pru.vector1_bool(oligomer.total_residue())
    for i in range(2, original_start + 1):
        v[i] = 1
    movemap.set_bb(v)
    movemap.set_chi(v)
    movemap.set_jump(False)
    relax.set_movemap(movemap)
    relax.apply(oligomer)
    info['dG'] = scorefxn(oligomer)
    monomer: pyrosetta.Pose = oligomer.split_by_chain(1)
    monomer.dump_pdb(str(relaxed_out_path / path.name))
    print(f'Relaxed {path.stem}', flush=True)
    if info['dG'] > -1e3:
        info['status']='score too high'
    else:
        # design
        print(f'Designing {path.stem}', flush=True)
        design_different(oligomer, ref, cycles=15, scorefxn=scorefxn)
        monomer = oligomer.split_by_chain(1)
        print('designed {name}', flush=True)
        monomer.dump_pdb(str(tuned_out_path / path.name))
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
    info['end'] = time.time()
    with Path('log.jsonl').open('a') as f:
        f.write(json.dumps(info) + '\n')
    return info



with pebble.ProcessPool(max_workers=os.cpu_count() - 1, max_tasks=0) as pool:
    pdb_paths = list(threaded_in_path.glob('*.pdb'))
    random.shuffle(pdb_paths)
    futuremap: pebble.ProcessMapFuture = pool.map(run, pdb_paths, timeout=timeout)
    iter_future: Iterator = futuremap.result()
    while True:
        try:
            result = next(iter_future)
            print(f"{result['name']}: {result['status']}")
        except StopIteration:
            break
        except Exception as error:
            warnings.warn(f"{error.__class__.__name__}: {str(error)}")
    print('Completed successfully!')
