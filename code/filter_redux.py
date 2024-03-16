"""
Filter polyGlycine templates against terrible clashes.

Align chain B (streptavidin) to chain B of the template
Extract chain A and replace in template
assess scored the complex
for each other streptavidin chain
    align chain B (streptavidin) to  of the template, copy over chain A
align CTD of chain A to chain F of the template
assess scored the complex
"""

from common_pyrosetta import (init_pyrosetta, get_pose_break, superpose_pose_by_chain, add_chain, rename_chain,
                              superpose_pose_by_alt_chains, superpose_by_seq_alignment,
                              extract_not_chainA
                              )

import multiprocessing
import operator
import os
import random
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
res_sele: ModuleType = pyrosetta.rosetta.core.select.residue_selector

output_folder = Path('output_redux/filtered')
out_path_pkl = Path('scored_complexes_redux.pkl.gz')
os.makedirs(output_folder, exist_ok=True)

logger = ph.configure_logger()
pyrosetta.distributed.maybe_init(extra_options=ph.make_option_string(no_optH=False,
                                                                     ex1=None,
                                                                     ex2=None,
                                                                     # mute='all',
                                                                     ignore_unrecognized_res=True,
                                                                     load_PDB_components=False,
                                                                     ignore_waters=True,
                                                                     # in=dict(in:detect_disulf=False)
                                                                     )
                                 )
pr_options.set_boolean_option('in:detect_disulf', False)


# ---- Load reference pose -----

# clasher.pdb is pentakaihemimer.relax.pdb but with two extra streptavidin chains, no paired chain of chain A
# basically aligning chain B to each of the 6 streptavidin chains will result in a corner.


# -------- funs using these globals

class Scorecard:
    notes_path = Path('notes.tsv')
    clasher = pyrosetta.pose_from_pdb('clasher.pdb')
    ref_chainB = clasher.split_by_chain(1)
    ref_woA = extract_not_chainA(pyrosetta.pose_from_file('pentakaihemimer.relax.pdb'))
    complex_delta_cutoff = 10_000
    rep_delta_cutoff = 4
    max_residue_score = 10

    def __init__(self, path: Path):
        self.name = path.stem
        self.path = path
        self.start_time = time.time()
        self.verdict = 'started'
        self.rmsd = float('nan')
        self.monomer_score = float('nan')
        self.complex_score = float('nan')
        self.residue_scores = []
        self.typesplit_complex_scores = {}
        self.neighbors = {}

    def note(self, verdict):
        # print(self.to_dict(), flush=True)
        with self.notes_path.open('a') as fh:
            self.verdict = verdict
            fh.write('\t'.join(map(str, self.to_dict().values())))

    @classmethod
    def get_seen(cls):
        # name\t
        return [l.split('\t')[0] for l in cls.notes_path.read_text().split('\n')]

    @property
    def time_taken(self):
        return time.time() - self.start_time

    def to_dict(self):
        return {'name': self.name,
                'path': str(self.path),
                'verdict': self.verdict,
                'monomer_score': self.monomer_score,
                'rmsd': self.rmsd,
                'complex_score': self.complex_score,
                'per_residue': self.residue_scores,
                'typesplit_complex_scores': self.typesplit_complex_scores,
                'neighbors': self.neighbors,
                'time': self.time_taken}


def run_process(path: Path, save_complex=False) -> dict:
    scorecard = Scorecard(path)
    if scorecard.name in scorecard.get_seen():
        scorecard.verdict = 'done_already'
        return scorecard.to_dict()
    try:
        # ---------------------------------
        # ## Load
        # ---------------------------------
        scorefxn: pr_scoring.ScoreFunction = pyrosetta.get_fa_scorefxn()
        # read and align
        pose: pyrosetta.Pose = pyrosetta.pose_from_file(path.as_posix())
        # ---------------------------------
        # ## Check for streptavidin
        # ---------------------------------
        seq: str = pose.sequence()
        # reject not streptavidin
        if 'HPQFEKR' not in seq:
            scorecard.note('no_streptavidin')
            return scorecard.to_dict()
        # ---------------------------------
        # ## Check for discontinuity
        # ---------------------------------
        disconnected_idx = get_pose_break(pose)
        if disconnected_idx:
            scorecard.note(f'disconnected {disconnected_idx}')
            return scorecard.to_dict()
        # ---------------------------------
        # ## Prep monomer (will be saved)
        # ---------------------------------
        scorecard.rmsd = superpose_by_seq_alignment(pose, scorecard.ref_chainB)
        if scorecard.rmsd > 1.:  # Something went wrong with alignment
            scorecard.note(f'bad_alignment {scorecard.rmsd}')
            return scorecard.to_dict()
        monomer: pyrosetta.Pose = pose.split_by_chain(1)
        scorecard.monomer_score = scorefxn(monomer)
        # ---------------------------------
        # ## Make complex
        # ---------------------------------
        complexo: pyrosetta.Pose = scorecard.clasher.clone()  # complex is a builtin!
        pi = complexo.pdb_info()
        # letter to number
        chain_ids: Dict[str, int] = {pi.chain(scorecard.clasher.chain_begin(i)): i for i in
                                     range(1, 1 + scorecard.clasher.num_chains())}
        for strep_chain, new_chain in zip('BCDEFG', 'AHIJKLMNO'):
            new = pose.clone()
            ref_chain = scorecard.clasher.split_by_chain(chain_ids[strep_chain])
            superpose_by_seq_alignment(new, ref_chain)
            neomono = new.split_by_chain(1)
            rename_chain(neomono, new_chain)
            add_chain(complexo, neomono)
        scorecard.complex_score = scorefxn(complexo)
        if save_complex:
            complexo.dump_pdb(str(output_folder / f'{scorecard.name}_complex.pdb'))
        # ---------------------------------
        # ## Check for clashes (global ddG)
        # ---------------------------------
        delta = scorecard.complex_score - (scorecard.monomer_score * 6)
        if delta > scorecard.complex_delta_cutoff:  # the complex is not perfectly made (glycine control scores 818)
            scorecard.note(f'clash {delta}')
            return scorecard.to_dict()
        # ---------------------------------
        # ## Check for high Lenard-Jones repulsion
        # ---------------------------------
        for st in scorefxn.get_nonzero_weighted_scoretypes():
            scorecard.typesplit_complex_scores[st.name] = scorefxn.score_by_scoretype(complexo, st)
        delta = scorecard.typesplit_complex_scores['fa_rep'] / (6 * monomer.total_residue()) * 0.55
        # 1.65 for glycine control
        if delta > scorecard.rep_delta_cutoff:
            scorecard.note(f'high_rep {delta}')
            return scorecard.to_dict()
        # ---------------------------------
        # ## Check for funky residues
        # ---------------------------------
        scorecard.residue_scores = []
        pi = complexo.pdb_info()
        chain_ids: Dict[str, int] = {pi.chain(complexo.chain_begin(i)): i for i in range(1, 1 + complexo.num_chains())}
        offset = complexo.chain_begin(chain_ids['A'])
        for i in range(monomer.total_residue()):
            v = pru.vector1_bool(complexo.total_residue())
            v[i + offset] = 1
            scorecard.residue_scores.append(scorefxn.get_sub_score(complexo, v))
        if max(scorecard.residue_scores) > scorecard.max_residue_score:
            scorecard.note(f'Disruptive residue')
            return scorecard.to_dict()

        # ---------------------------------
        # ## Save
        # ---------------------------------
        monomer.dump_pdb(str(output_folder / f'{scorecard.name}.pdb'))
        scorecard.note('saved')
        pi = complexo.pdb_info()
        chainA_sele = res_sele.ChainSelector('A')
        for distance in (1, 2, 4, 5, 6, 8, 10, 12):
            sele = res_sele.NeighborhoodResidueSelector(chainA_sele, distance, False)
            close_residues = prc.select.get_residues_from_subset( sele.apply(complexo))
            scorecard.neighbors[f'close_residues_{distance}'] = [pi.pose2pdb(r) for r in close_residues]
        scorecard.note('complete')
        return scorecard.to_dict()
    except Exception as error:
        scorecard.note(f'crashed {error.__class__.__name__}: {error}')
        return scorecard.to_dict()


scores = []


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
            print(f"Function raised {error.__class__.__name__}: {error_msg}")
        result = dict(error=error.__class__.__name__, error_msg=error_msg)
    scores.append(result)


if __name__ == '__main__':
    group_seen_already = []
    test_one_of_each = False
    # run
    num_cores = multiprocessing.cpu_count()
    with ProcessPool(max_workers=50, max_tasks=0) as pool:
        paths = list(Path('output').glob('*/*.pdb'))
        random.shuffle(paths)
        for path in paths:
            name = path.stem
            group = name.split('_')[0]
            if 'trash' in str(path) or 'traj' in str(path):
                continue
            elif test_one_of_each and group in group_seen_already:
                continue
            elif group not in ('theta', 'eta', 'iota', 'kappa', 'lambda', 'mu', 'nu'):
                continue
            group_seen_already.append(group)
            future = pool.schedule(run_process, [path], timeout=360)
            future.add_done_callback(run_done)

    df = pd.DataFrame(scores)
    df = df.set_index('name')
    df.to_pickle(out_path_pkl)
    print('Completed successfully')
