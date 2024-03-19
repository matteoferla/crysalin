"""
Given a shortlist, run FastDesign and scoring on each of them.
"""


from common_pyrosetta import (init_pyrosetta, design_interface_onA, superpose_pose_by_chain,
                              relax_chainA,transpant_CTD,constrain_chainbreak,freeze_atom,design_different,
                              extract_not_chainA, add_chain, extract_streptavidins, score_interface)
import multiprocessing
import time
from pathlib import Path
import pebble
import os
import random
import json
from typing import Iterator
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
    oligomer: pyrosetta.Pose = transpant_CTD(monomer, ref, 'KDETET')  # for semantics
    add_chain(oligomer, woA)  # as this is in place
    for i in range(1, monomer.sequence().find('KDETET') + 1):
        constrain_chainbreak(oligomer, i)
    ref_index = oligomer.chain_begin(2)
    freeze_atom(pose=oligomer, frozen_index=oligomer.chain_end(1), ref_index=ref_index)
    freeze_atom(pose=oligomer, frozen_index=oligomer.chain_end(2), ref_index=ref_index)
    for i in range(3, oligomer.num_chains() + 1):
        freeze_atom(pose=oligomer, frozen_index=oligomer.chain_begin(i), ref_index=ref_index)
        freeze_atom(pose=oligomer, frozen_index=oligomer.chain_end(i), ref_index=ref_index)
    scorefxn: pr_scoring.ScoreFunction = pyrosetta.get_fa_scorefxn()
    scorefxn.set_weight(pr_scoring.ScoreType.atom_pair_constraint, 3)
    scorefxn.set_weight(pr_scoring.ScoreType.coordinate_constraint, 5)
    relax_chainA(oligomer, cycles=5, distance=0, scorefxn=scorefxn)  # only chain A
    info['complex_dG_pre'] = scorefxn(oligomer)
    monomer = oligomer.split_by_chain(1)
    info['monomer_dG_pre'] = scorefxn(monomer)
    monomer.dump_pdb(f'{outfolder}/relaxed/{pdb_path.stem}.pdb')
    design_different(oligomer, ref, cycles=15, scorefxn=scorefxn)
    monomer = oligomer.split_by_chain(1)
    monomer.dump_pdb(f'{outfolder}/tuned/{pdb_path.stem}.pdb')
    info['complex_dG_post'] = scorefxn(oligomer)
    info['monomer_dG_post'] = scorefxn(monomer)
    info['tweaked_sequence'] = monomer.sequence()
    streptavidins = extract_streptavidins(ref, cardinality=2)
    superpose_pose_by_chain(oligomer, streptavidins, 'B')
    # combine with streptavidins only
    minicomplex = monomer.clone()
    add_chain(minicomplex, streptavidins)
    info.update(**score_interface(minicomplex, 'A_BC'))
    info['end'] = time.time()
    info['status'] = 'done'
    Path(f'{outfolder}/info/{name}.json').write_text(json.dumps(info))
    return info

if __name__ == '__main__':
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
