from common_pyrosetta import (init_pyrosetta, design_interface_onA, superpose_pose_by_chain,
                              relax_chainA,
                              extract_not_chainA, add_chain, extract_streptavidins, score_interface)
import multiprocessing
import time
from pathlib import Path
import pebble
import os
import random
import json
from typing import Iterator
import warnings
import pyrosetta

timeout = 4 * 60 * 60
outfolder = 'output_scoring'
os.makedirs(f'{outfolder}/tuned', exist_ok=True)
print('starting Pyrosetta...')
init_pyrosetta()
num_cores = multiprocessing.cpu_count()


def parse(pdb_path: Path) -> dict:
    name = pdb_path.stem  # same as `str(pdb_path).split('/')[-1].split('.')[0]`
    if Path(f'{outfolder}/tuned/{pdb_path.stem}.pdb').exists():
        return dict(name=name, status='done already')
    Path(f'{outfolder}/tuned/{pdb_path.stem}.pdb').touch()
    init_pyrosetta()
    scorefxn = pyrosetta.get_fa_scorefxn()
    info = json.loads(Path(f'{outfolder}/info/{name}.json').read_text())
    info['status'] = 'crashed_mysteriously'
    pdbblock: str = Path(pdb_path).read_text()
    if not pdbblock.strip():
        info['error'] = 'empty'
        info['end'] = time.time()
        print(name, 'empty')
        return info
    monomer: pyrosetta.Pose = pyrosetta.pose_from_file(str(pdb_path)).split_by_chain(1)
    ref = pyrosetta.pose_from_file('pentakaihemimer.relax.pdb')
    woA = extract_not_chainA(ref)
    complex = monomer  # for semantics
    add_chain(complex, woA)  # as this is in place
    relax_chainA(complex, cycles=5)
    monomer = complex.split_by_chain(1)
    monomer.dump_pdb(f'{outfolder}/relaxed/{pdb_path.stem}.pdb')
    streptavidins = extract_streptavidins(ref, cardinality=2)
    superpose_pose_by_chain(complex, streptavidins, 'B')
    # combine with streptavidins only
    minicomplex = monomer.clone()
    add_chain(minicomplex, streptavidins)
    info['minicomplex_dG_pre'] = scorefxn(minicomplex)
    design_interface_onA(minicomplex, distance=7, cycles=15)
    superpose_pose_by_chain(minicomplex, streptavidins, 'B')
    monomer = minicomplex.split_by_chain(1)
    monomer.dump_pdb(f'{outfolder}/tuned/{pdb_path.stem}.pdb')
    info['minicomplex_dG_post'] = scorefxn(minicomplex)
    info['tweaked_sequence'] = monomer.sequence()
    info.update(**score_interface(minicomplex, 'A_BC'))
    info['end'] = time.time()
    info['status'] = 'done'
    Path(f'{outfolder}/info/{name}.json').write_text(json.dumps(info))
    return info


with pebble.ProcessPool(max_workers=num_cores - 1, max_tasks=0) as pool:
    names = json.loads(Path('shortlist.json').read_text())
    pdbs_paths = [Path(f'{outfolder}/relaxed/{name}.pdb') for name in names]
    random.shuffle(pdbs_paths)
    futuremap: pebble.ProcessMapFuture = pool.map(parse, pdbs_paths, timeout=timeout)
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
