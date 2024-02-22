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


def align_crop(pose: pyrosetta.Pose, ref: pyrosetta.Pose, verbose=False) -> float:
    """
    Superpose ``pose`` on ``ref`` based on Pairwise alignment and superposition of CA
    Then crops all residues starting from it. In place.
    """
    # pad with '-' to make it faux-local alignment and deal with Fortran counting does not work '-' is a match not gap
    # hence the silly +1s and the PairwiseAligner settings
    aligner = PairwiseAligner()
    aligner.internal_gap_score = -10
    aligner.extend_gap_score = -0.01
    aligner.end_gap_score = -0.01
    # pose is longer and right does not matter. left aligned!
    aligner.target_right_gap_score = 0.
    aligner.target_right_extend_gap_score = 0.
    ref_seq: str = ref.sequence()
    pose_seq: str = pose.sequence()
    aln: Alignment = aligner.align(ref_seq, pose_seq)[0]
    if verbose:
        print(aln)
    aln_map = {t: q for t, q in zip(aln.indices[0], aln.indices[1]) if
               q != -1 and t != -1 and ref_seq[t] == pose_seq[q]}
    if verbose:
        print(aln_map)
    # ## make pyrosetta map
    atom_map = pyrosetta.rosetta.std.map_core_id_AtomID_core_id_AtomID()
    for r, m in aln_map.items():
        ref_atom = pyrosetta.AtomID(ref.residue(r + 1).atom_index("CA"), r + 1)
        mobile_atom = pyrosetta.AtomID(pose.residue(m + 1).atom_index("CA"), m + 1)
        atom_map[mobile_atom] = ref_atom
    # return RMSD
    rmsd = prc.scoring.superimpose_pose(mod_pose=pose, ref_pose=ref, atom_map=atom_map)
    pyrosetta.rosetta.protocols.grafting.delete_region(pose, next(iter(aln_map.values())), pose.total_residue())
    return rmsd


class Analyse:
    def __init__(self, ref: pyrosetta.Pose, others: pyrosetta.Pose):
        self.others = others
        self.disulfides2alanine(others)
        self.disulfides2alanine(ref)
        pi = ref.pdb_info()
        self._ref_info: List[str] = [''] + [pi.pose2pdb(i) for i in range(1, 1 + ref.total_residue())]
        self.ref_chainA: pyrosetta.Pose = ref.split_by_chain(1)
        self._ref_chainA_len: int = self.ref_chainA.total_residue()

    @staticmethod
    def disulfides2alanine(pose):
        """
        Replaces disulfides with alanines.
        """
        disulfide_res = prc.select.residue_selector.ResidueNameSelector("CYS:disulfide")
        for resi in prc.select.get_residues_from_subset(disulfide_res.apply(pose)):
            prp.simple_moves.MutateResidue(target=resi, new_res='ALA').apply(pose)

    @staticmethod
    def add_chain(built: pyrosetta.Pose, new: pyrosetta.Pose) -> None:
        """
        Add a chain ``new`` to a pose ``built`` preserving the residue numbering.
        """
        offset: int = built.total_residue()
        for chain in new.split_by_chain():
            pyrosetta.rosetta.core.pose.append_pose_to_pose(built, chain, new_chain=True)
            built_pi = built.pdb_info()
            chain_pi = chain.pdb_info()
            for r in range(1, chain.total_residue() + 1):
                built_pi.set_resinfo(res=r + offset, chain_id=chain_pi.chain(r), pdb_res=chain_pi.number(r))

    @classmethod
    def make_woA(cls, original_pose: pyrosetta.Pose) -> pyrosetta.Pose:
        chains: pru.vector1_std_shared_ptr_core_pose_Pose_t = original_pose.split_by_chain()
        others: pyrosetta.Pose = chains[2].clone()
        for c in range(3, len(chains) + 1):
            cls.add_chain(others, chains[c])
        return others

    def get_ref_info(self, resi: int, new_len: int) -> str:
        """
        Converts the residue number to the PDB numbering of the reference.
        """
        return self._ref_info[resi + new_len - self._ref_chainA_len]

    def __call__(self, pose: pyrosetta.Pose, name: str) -> Dict[str, Union[int, float, List[Dict[str, str]]]]:
        """
        Assumes pose is aligned already.
        """
        scorefxn = pyrosetta.get_fa_scorefxn()
        new_len: int = pose.total_residue()
        sequence = pose.sequence()
        info = {'name': name, 'length': new_len,
                'monomer_score': scorefxn(pose), 'sequence': sequence,
                'pI': ProtParam.ProteinAnalysis(sequence).isoelectric_point(),
                'error': '', 'error_msg': ''}
        complex: pyrosetta.Pose = pose.clone()
        self.add_chain(complex, self.others)
        self.disulfides2alanine(complex)
        info['complex_score'] = scorefxn(complex)
        res_sele = prc.select.residue_selector
        chainA = res_sele.ChainSelector('A')
        for distance in (4, 5, 6, 8, 10, 12):
            close_residues = prc.select.get_residues_from_subset(
                res_sele.NeighborhoodResidueSelector(chainA, distance, False).apply(complex))
            info[f'N_close_residues_{distance}'] = len(close_residues)
            info[f'close_residues_{distance}'] = [self.get_ref_info(r, new_len) for r in close_residues]
        # info['complex'] = test
        # complex.dump_pdb(f'/data/outerhome/tmp/complexes/{name}_complex.pdb')
        return info


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
df['group'] = df.index.to_series().str.split('_').apply(operator.itemgetter(0))
df['subgroup'] = df.index.to_series().apply(lambda name: '_'.join(name.split('_')[:-1]))

df[[c for c in df.columns if c not in ('monomer', 'minicomplex')]].to_csv(out_path_csv)
df.to_pickle(out_path_pkl)

