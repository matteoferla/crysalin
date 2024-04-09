import os
import operator
import string
from Bio import SeqIO
from pathlib import Path
import re, os, traceback, multiprocessing
import pebble
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
from common_pyrosetta import (add_chain, superpose, thread, superpose_pose_by_chain,
                              constrain_chainbreak, freeze_atom, relax_chainA, init_pyrosetta)
#input_folder = 'output_redux/filtered'
input_folders = ['output/omicron', 'output/rho', 'output/sigma']
work_path = Path(os.environ.get('WORKPATH', 'output_redux'))
raw_out_path = work_path / 'superposed_skeletons'
os.makedirs(raw_out_path, exist_ok=True)

# --------------------------------
print('\n## Init PyRosetta\n')
init_pyrosetta()
num_cores = multiprocessing.cpu_count()
futures: List[pebble.ProcessFuture] = []
timeout = 60 * 60 * 5
ref: pyrosetta.Pose = pyrosetta.pose_from_file('pentakaihemimer.relax.pdb')

# --------------------------------

def run_process(path):
    pose: pyrosetta.Pose = pyrosetta.pose_from_file(str(path))
    rmsd = superpose_pose_by_chain(pose, ref, chain='B', strict=True)
    pose.split_by_chain(1).dump_pdb(f'{raw_out_path}/{path.stem}.pdb')
    print(f'{path.stem} done ({rmsd:.2f} A)')

# --------------------------------
print('\n## Running\n')

n = 10  # <-- positive number: for testing
paths = [path for input_folder in input_folders for path in Path(input_folder).glob('*.pdb')]
with pebble.ProcessPool(max_workers=os.cpu_count() - 1, max_tasks=0) as pool:
    # queue jobs
    for path in paths:
        if Path(f'{raw_out_path}/{path.stem}.pdb').exists():
            print(f'{path.stem} done already')
            continue
        n -= 1
        if n == 0:
            print('break')
            break
        future: pebble.ProcessFuture = pool.schedule(run_process,
                                              kwargs=dict(path=path),
                                              timeout=timeout)
        futures.append(future)
    print(f'Submitted {len(futures)} processes')
    # get results
    for future in futures:
        try:
            result = future.result()  # blocks until results are ready
        except Exception as error:
            error_type = type(error).__name__
            error_msg = str(error)
            if isinstance(error, TimeoutError):
                print(f'Function took longer than {timeout} seconds {error}')
            else:
                print(f"Function raised {error_type} {error_msg}")
                traceback.print_tb(error.__traceback__)
