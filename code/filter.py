from common_pyrosetta import *

import multiprocessing
import operator
import os
import time
from concurrent.futures import TimeoutError
from pathlib import Path
from types import ModuleType
from typing import List, Dict, Union
import pandas as pd
import pyrosetta
import pyrosetta_help as ph
from Bio.Align import PairwiseAligner, Alignment
from Bio.SeqUtils import ProtParam
from pebble import ProcessPool

prc: ModuleType = pyrosetta.rosetta.core
prp: ModuleType = pyrosetta.rosetta.protocols
pru: ModuleType = pyrosetta.rosetta.utility
prn: ModuleType = pyrosetta.rosetta.numeric
pr_conf: ModuleType = pyrosetta.rosetta.core.conformation
pr_scoring: ModuleType = pyrosetta.rosetta.core.scoring
pr_options: ModuleType = pyrosetta.rosetta.basic.options

out_path_pkl = Path('scored_complexes3.pkl.gz')
out_path_csv = Path('scored_complexes3.csv')
previous_path = Path('scored_complexes2.pkl.gz')

logger = ph.configure_logger()
pyrosetta.distributed.maybe_init(extra_options=ph.make_option_string(no_optH=False,
                                                                     ex1=None,
                                                                     ex2=None,
                                                                     # mute='all',
                                                                     ignore_unrecognized_res=True,
                                                                     load_PDB_components=False,
                                                                     ignore_waters=True,
                                                                     # in=dict(in:detect_disulf=True)
                                                                     )
                                 )
pr_options.set_boolean_option('in:detect_disulf', True)

# ---- Load reference pose -----

others: pyrosetta.Pose = pyrosetta.pose_from_file('woA_crysalin_lattice.pdb')
Analyse.disulfides2alanine(others)
ref_chainB = others.split_by_chain(1)
print('chain B', ref_chainB.total_residue(), ref_chainB.total_atoms())

ref: pyrosetta.Pose = pyrosetta.pose_from_file('crysalin_lattice.pdb')
Analyse.disulfides2alanine(ref)

trihemi = pyrosetta.pose_from_pdb('trikaihemimer.pdb')
mini_others = Analyse.make_woA(trihemi)

analyse = Analyse(ref=ref, others=others)

score = {'name': 'reference', 'path': '', 'error': '', 'error_msg': '', **analyse(analyse.ref_chainA, 'reference')}
scores = [score]
seen = []
test_one_of_each = False
num_cores = multiprocessing.cpu_count()


# -------- funs using these globals


def run_process(path: Path):
    name = os.path.splitext(path.name)[0]
    score = {'name': name, 'path': path}
    tick = time.time()
    # read and align
    pose: pyrosetta.Pose = pyrosetta.pose_from_file(path.as_posix())
    rmsd = align_crop(pose, ref_chainB)
    score['rmsd'] = rmsd
    assert rmsd < 1., 'Something went wrong with alignment'
    mini: pyrosetta.Pose = pose.clone()
    Analyse.add_chain(mini, mini_others)
    if test_one_of_each:
        mini.dump_pdb(f'{name}_minicomplex.pdb')
    score.update(**analyse(pose, name))
    tock = time.time()
    score['monomer'] = ph.get_pdbstr(pose)
    score['minicomplex'] = ph.get_pdbstr(mini)
    score['time'] = tock - tick
    return score


def run_done(future):
    nan = float('nan')
    try:
        result = future.result()  # blocks until results are ready
        name = result['name']
        timed = result['time']
        print(f'{name} done in {timed}s')
    except Exception as error:
        error_msg = str(error)
        if isinstance(error, TimeoutError):
            print(f'Function took longer than 240 seconds {error}')
        else:
            print(f"Function raised {error}")
            print(error.traceback)  # traceback of the function
        result = dict(error=error.__class__.__name__, error_msg=error_msg)
    scores.append(result)


# previous run
seen_already: List[str] = []
if previous_path.exists():
    df = pd.read_pickle(previous_path)
    seen_already = df.path.to_list()

# run
with ProcessPool(max_workers=50, max_tasks=0) as pool:
    for path in Path('output').glob('*/*.pdb'):
        if path.as_posix() in seen_already:
            continue
        name = os.path.splitext(path.name)[0]
        group = name.split('_')[0]
        if 'trash' in str(path) or 'traj' in str(path):
            continue
        elif test_one_of_each and group in seen:
            continue
        elif group not in ('iota', 'eta', 'theta'):
            continue
        seen.append(group)
        future = pool.schedule(run_process, [path], timeout=360)
        future.add_done_callback(run_done)

df = pd.DataFrame(scores)
df = df.set_index('name')
df['group'] = df.index.to_series().str.split('_', expand=True)[0]
df['subgroup'] = df.index.to_series().apply(lambda name: '_'.join(name.split('_')[:-1]))

df[[c for c in df.columns if c not in ('monomer', 'minicomplex')]].to_csv(out_path_csv)
df.to_pickle(out_path_pkl)
