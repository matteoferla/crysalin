# -*- coding: utf-8 -*-
# global:
do_final_relax = False  # do a final relax and backrub — takes hours
folder_names = ['output_MPNN/pdbs','output_MPNN2/pdbs']
subset = 0
timeout = 4 * 60 * 60

import itertools
import json
import multiprocessing
import time
from pathlib import Path
from types import ModuleType
from typing import Iterator
import pebble
import pyrosetta
import pyrosetta_help as ph

from common_pyrosetta import init_pyrosetta, extract_not_chainA, extract_wo_strep, add_chain, movement, get_pose_break, \
    extract_streps

prc: ModuleType = pyrosetta.rosetta.core
prp: ModuleType = pyrosetta.rosetta.protocols
pru: ModuleType = pyrosetta.rosetta.utility  # noqa
prn: ModuleType = pyrosetta.rosetta.numeric
pr_conf: ModuleType = pyrosetta.rosetta.core.conformation
pr_scoring: ModuleType = pyrosetta.rosetta.core.scoring
pr_options: ModuleType = pyrosetta.rosetta.basic.options
pr_res: ModuleType = pyrosetta.rosetta.core.select.residue_selector


def parse(pdb_path: Path):
    # it is already aligned.
    # align_pose_to_ref_by_chain(pose, ref, 'B')
    # it is already chain A only
    # chainA = pose.split_by_chain(1)
    name = str(pdb_path).split('/')[-1].split('.')[0]
    def checkpoint():
        with open(f'output_MPNN/info/{name}.json', 'w') as f:
            json.dump(info, f)

    if Path(f'output_MPNN/info/{pdb_path.stem}.json').exists():
        return dict(name=name, status='done already')
    init_pyrosetta()
    scorefxn = pyrosetta.get_fa_scorefxn()
    info = dict(name=name, path=str(pdb_path), start=time.time())
    info['status'] = 'crashed_mysteriously'
    checkpoint()
    pdbblock: str = Path(pdb_path).read_text()
    if not pdbblock.strip():
        info['error'] = 'empty'
        info['end'] = time.time()
        print(name, 'empty')
        return info
    monomer: pyrosetta.Pose = pyrosetta.pose_from_file(str(pdb_path)).split_by_chain(1)
    # sanitize pdb info
    pi: prc.pose.PDBInfo = monomer.pdb_info()
    if pi is None or pi.chain(1) != 'A':
        for i in range(1, monomer.total_residue() + 1):
            pi.set_resinfo(res=i, chain_id='A', pdb_res=i)
    pi.obsolete(False)
    # store sequence
    info['sequence'] = monomer.sequence()
    info['mono_dG'] = scorefxn(monomer)
    info['status'] = 'loaded'
    info['chainbreak'] = get_pose_break(monomer)
    if info['chainbreak']:
        info['status'] = 'chainbreak'
        checkpoint()
        print(name, 'chainbreak')
        return info
    else:
        checkpoint()
    ref = pyrosetta.pose_from_file('pentakaihemimer.relax.pdb')
    woA = extract_not_chainA(ref)
    # combine with other parts
    complex = monomer.clone()
    add_chain(complex, woA)
    chainA_sele = prc.select.residue_selector.ChainSelector('A')
    relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, 5)
    movemap = pyrosetta.MoveMap()
    movemap.set_bb(False)
    movemap.set_chi(allow_chi=chainA_sele.apply(complex))
    relax.set_movemap(movemap)
    relax.apply(complex)
    complex.split_by_chain(1).dump_pdb(f'output_relaxed/{name}.pdb')
    info['complex_dG'] = scorefxn(complex)
    woS = extract_wo_strep(complex)
    info['wo_strep_dG'] = scorefxn(woS)
    info['strep_dG'] = scorefxn(extract_streps(complex))
    # get number of close contacts
    pose2pdb = complex.pdb_info().pose2pdb
    for distance in (3, 4, 5, 6, 8, 10, 12):
        v = prc.select.residue_selector.NeighborhoodResidueSelector(chainA_sele, distance, False).apply(complex)
        close_residues = prc.select.get_residues_from_subset(v)
        info[f'N_close_residues_{distance}'] = len(close_residues)
        info[f'close_residues_{distance}'] = [pose2pdb(r) for r in close_residues]
    # per residue
    residue_scores = []
    for i in ph.pose_range(monomer):
        v = pru.vector1_bool(monomer.total_residue())
        v[i] = 1
        residue_scores.append(scorefxn.get_sub_score(monomer, v))
    info['per_residue'] = residue_scores
    # pdbblock = ammend_pdbblock(pdbblock, info)
    info['end'] = time.time()
    info['interface'] = info['complex_dG'] - info['wo_strep_dG'] - 318.5261163127267
    if max(residue_scores) > 10:
        info['status'] = 'high_residue_score'
        checkpoint()
        print(name, 'high_residue_score')
        return info
    elif info['interface'] > -400:
        info['status'] = 'weaker_interface'
        checkpoint()
        print(name, 'weaker_interface')
        return info
    else:
        checkpoint()
    # relax & backrub
    if not do_final_relax:
        return info
    relaxed_monomer = monomer.clone()
    pyrosetta.rosetta.protocols.relax.FastRelax.register_options()  # noqa
    relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, 5)
    relax.apply(relaxed_monomer)
    info['relaxed_monomer_dG'] = scorefxn(relaxed_monomer)
    relax.apply(complex)
    info['relaxed_complex_dG'] = scorefxn(complex)
    info['relaxed_complex_dG_chainA'] = scorefxn(complex.split_by_chain(1))
    info['relaxed_complex_dG_noA'] = scorefxn(extract_not_chainA(complex))
    info['relaxed_complex_dG_woS'] = scorefxn(extract_wo_strep(complex))
    relax.apply(woS)
    info['relaxed_wo_strep_dG'] = scorefxn(woS)
    info['monomer_backrub_movement'] = movement(relaxed_monomer)
    info['end'] = time.time()
    info['status'] = 'complete'
    checkpoint()
    print(name, 'complete')
    return info


# ---------------------------------------------------------------------------------
import random

folders = [name for name in folder_names]
pdbs_paths = itertools.chain.from_iterable([Path(folder).glob('*.pdb') for folder in folders])
pdbs_paths = list(pdbs_paths)
random.shuffle(pdbs_paths)
if subset:
    pdbs_paths = pdbs_paths[:1000]

print('starting Pyrosetta...')
init_pyrosetta()
num_cores = multiprocessing.cpu_count()

with pebble.ProcessPool(max_workers=num_cores - 1, max_tasks=0) as pool:
    print(f'starting pool... {pdbs_paths}')
    futuremap: pebble.ProcessMapFuture = pool.map(parse, pdbs_paths, timeout=timeout)
    iter_future: Iterator = futuremap.result()
    while True:
        try:
            result = next(iter_future)
            print(result['name'], result['status'])
        except StopIteration:
            break
        except Exception as error:
            print(error.__class__.__name__, str(error))
print('Completed successfully!')
