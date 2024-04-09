import re
from pathlib import Path
from types import ModuleType
from typing import List, Dict, Union, Optional, Sequence
import numpy as np
import pandas as pd
import pyrosetta
import pyrosetta_help as ph
from Bio.Align import PairwiseAligner, Alignment
from Bio.SeqUtils import ProtParam

prc: ModuleType = pyrosetta.rosetta.core
prp: ModuleType = pyrosetta.rosetta.protocols
pru: ModuleType = pyrosetta.rosetta.utility  # noqa
prn: ModuleType = pyrosetta.rosetta.numeric
prs: ModuleType = pyrosetta.rosetta.std  # noqa
pr_conf: ModuleType = pyrosetta.rosetta.core.conformation
pr_scoring: ModuleType = pyrosetta.rosetta.core.scoring
pr_options: ModuleType = pyrosetta.rosetta.basic.options
pr_res: ModuleType = pyrosetta.rosetta.core.select.residue_selector


# --------------------------------

def init_pyrosetta(detect_disulf=False):
    logger = ph.configure_logger()
    pyrosetta.distributed.maybe_init(extra_options=ph.make_option_string(no_optH=False,
                                                                         ex1=None,
                                                                         ex2=None,
                                                                         # mute='all',
                                                                         ignore_unrecognized_res=True,
                                                                         load_PDB_components=False,
                                                                         ignore_waters=True,
                                                                         # in=dict(in:detect_disulf=detect_disulf)
                                                                         )
                                     )
    pr_options.set_boolean_option('in:detect_disulf', detect_disulf)
    return logger


def add_chain(built: pyrosetta.Pose, new: pyrosetta.Pose, reset: bool = False) -> None:
    """
    Add a chain ``new`` to a pose ``built`` preserving the residue numbering.

    :param built: this is the pyrosetta.Pose that will be built into...
    :param new: the addendum
    :param reset: resets the PDBInfo for the chain present to A
    :return:
    """
    built_pi = built.pdb_info()
    if built_pi is None or reset:
        built_pi = prc.pose.PDBInfo(built)
        built.pdb_info(built_pi)
        for r in range(1, built.total_residue() + 1):
            built_pi.set_resinfo(res=r, chain_id='A', pdb_res=r)
    for chain in new.split_by_chain():
        offset: int = built.total_residue()
        pyrosetta.rosetta.core.pose.append_pose_to_pose(built, chain, new_chain=True)
        chain_pi = chain.pdb_info()
        for r in range(1, chain.total_residue() + 1):
            built_pi.set_resinfo(res=r + offset, chain_id=chain_pi.chain(r), pdb_res=chain_pi.number(r))
    built_pi.obsolete(False)



def fix_pdb(pose: pyrosetta):
    """
    Very crude fix...

    :param pose:
    :return:
    """
    pi: prc.pose.PDBInfo = pose.pdb_info()
    previous_max = 0
    previous_chain = 'A'
    for r in range(1, pose.total_residue() + 1):
        chain_letter = chr(64 + pose.chain(r))
        if chain_letter != previous_chain:
            previous_max = r - 1
            previous_chain = chain_letter
        pi.chain(r, chain_letter)
        pi.number(r, r - previous_max)
    pi.obsolete(False)

def rename_chain(pose: pyrosetta.Pose, chain: str) -> None:
    pi: prc.pose.PDBInfo = pose.pdb_info()
    for r in range(1, pose.total_residue() + 1):
        pi.chain(r, chain)
    pi.obsolete(False)

# ------------------------------------------------------------------------------------
# Superposition

def align_for_atom_map(mobile: pyrosetta.Pose, ref: pyrosetta.Pose) -> Dict[int, int]:
    """
    Pairwise alignment of the sequences of the poses.
    return  (ref_index, mobile_index)
    :param mobile:
    :param ref:
    :return:
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
    pose_seq: str = mobile.sequence()
    aln: Alignment = aligner.align(ref_seq, pose_seq)[0]
    return {t: q for t, q in zip(aln.indices[0], aln.indices[1]) if
               q != -1 and t != -1 and ref_seq[t] == pose_seq[q]}


def superpose_pose_by_chain(pose, ref, chain: str, strict: bool=True) -> float:
    """
    superpose by PDB chain letter

    :param pose:
    :param ref:
    :param chain:
    :return:
    """
    atom_map = prs.map_core_id_AtomID_core_id_AtomID()
    chain_sele: pr_res.ResidueSelector = pr_res.ChainSelector(chain)
    for r, m in zip(pr_res.selection_positions(chain_sele.apply(ref)),
                    pr_res.selection_positions(chain_sele.apply(pose))
                    ):
        if strict:
            assert pose.residue(m).name3() == ref.residue(r).name3(), 'Mismatching residue positions!'
        ref_atom = pyrosetta.AtomID(ref.residue(r).atom_index("CA"), r)
        mobile_atom = pyrosetta.AtomID(pose.residue(m).atom_index("CA"), m)
        atom_map[mobile_atom] = ref_atom
    return prc.scoring.superimpose_pose(mod_pose=pose, ref_pose=ref, atom_map=atom_map)

def superpose_pose_by_alt_chains(pose, ref, pose_chain: str, ref_chain: str) -> float:
    """
    superpose by PDB chain letter

    :param pose:
    :param ref:
    :param chain:
    :return:
    """
    atom_map = prs.map_core_id_AtomID_core_id_AtomID()
    pose_chain_sele: pr_res.ResidueSelector = pr_res.ChainSelector(pose_chain)
    ref_chain_sele: pr_res.ResidueSelector = pr_res.ChainSelector(ref_chain)
    for r, m in zip(pr_res.selection_positions(ref_chain_sele.apply(ref)),
                    pr_res.selection_positions(pose_chain_sele.apply(pose))
                    ):
        assert pose.residue(m) == ref.residue(r), 'Mismatching residue positions!'
        ref_atom = pyrosetta.AtomID(ref.residue(r).atom_index("CA"), r)
        mobile_atom = pyrosetta.AtomID(pose.residue(m).atom_index("CA"), m)
        atom_map[mobile_atom] = ref_atom
    return prc.scoring.superimpose_pose(mod_pose=pose, ref_pose=ref, atom_map=atom_map)

def superpose_by_seq_alignment(mobile: pyrosetta.Pose, ref: pyrosetta.Pose) -> float:
    """
    Superpose ``pose`` on ``ref`` based on Pairwise alignment and superposition of CA

    :param mobile:
    :param ref:
    :param verbose:
    :return:
    """

    aln_map = align_for_atom_map(mobile, ref)
    rmsd: float = superpose(ref=ref, mobile=mobile, aln_map=aln_map)
    return rmsd

def superpose(ref: pyrosetta.Pose, mobile: pyrosetta.Pose, aln_map: Optional[Dict[int, int]] = None, zero_based=True) -> float:
    """
    Superpose ``mobile`` on ``ref`` based on CA of indices in ``aln_map`` (ref_indices, mobile_indices).
    Indices are 0-based.
    :param ref:
    :param mobile:
    :param aln_map:
    :return:
    """
    offset = 1 if zero_based else 0
    if aln_map is None:
        aln_map = dict(zip(range(ref.total_residue()), range(mobile.total_residue())))
    # ## make pyrosetta map
    atom_map = prs.map_core_id_AtomID_core_id_AtomID()
    for r, m in aln_map.items():
        ref_atom = pyrosetta.AtomID(ref.residue(r + offset).atom_index("CA"), r + offset)
        mobile_atom = pyrosetta.AtomID(mobile.residue(m + offset).atom_index("CA"), m + offset)
        atom_map[mobile_atom] = ref_atom
    # return RMSD
    return prc.scoring.superimpose_pose(mod_pose=mobile, ref_pose=ref, atom_map=atom_map)

def get_ATOM_only(pdbblock: str) -> str:
    return '\n'.join([line for line in pdbblock.splitlines() if line.startswith('ATOM')])

three_to_one = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',
    'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
    'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
    'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}
def get_sequence(pdbblock: str) -> str:
    sequence = ''
    residues_seen = set()
    for line in pdbblock.splitlines():
        if line.startswith("ATOM") and " CA " in line:
            res_info = line[17:26]  # Residue name and number for uniqueness
            if res_info not in residues_seen:
                residues_seen.add(res_info)
                res_name = line[17:20].strip()
                sequence += three_to_one.get(res_name, '?')
    return sequence


def unused_relax(pose: pyrosetta.Pose, others, cycles=1):
    # unused?
    scorefxn: pr_scoring.ScoreFunction = pyrosetta.get_fa_scorefxn()
    dG_others = scorefxn(others)
    rs: ModuleType = prc.select.residue_selector
    chainA_sele: rs.ResidueSelectorr = rs.ChainSelector('A')
    chainA: pru.vector1_bool = chainA_sele.apply(pose)
    neigh_sele: rs.ResidueSelector = prc.select.residue_selector.NeighborhoodResidueSelector(chainA_sele, True, 5)
    neighs: pru.vector1_bool = neigh_sele.apply(pose)
    relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, cycles)
    movemap = pyrosetta.MoveMap()
    movemap.set_bb(chainA)
    movemap.set_chi(neighs)
    movemap.set_jump(chainA)
    relax.set_movemap(movemap)
    relax.apply(pose)
    return scorefxn(pose) - dG_others

# ------------------------------------------------------------------------------------
# streptavidin extraction

def extract_wo_streptavidins(pose):
    chains = pose.split_by_chain()
    wo = chains[1]
    for i, c in enumerate(chains):
        if i == 0:
            continue
        if 'MEAGIT' in c.sequence():
            continue
        add_chain(wo, c)
    return wo


def extract_streptavidins(pose, cardinality=4):
    """
    This will extract the tetramer (``cardinality=4``) or dimer (``cardinality=4``)
    """
    chains = pose.split_by_chain()
    streptavidins: Union[None, pyrosetta.Pose] = None
    for i, c in enumerate(chains):
        if cardinality == 0:
            break
        elif 'MEAGIT' not in c.sequence():
            continue
        elif streptavidins is None:
            streptavidins = c
            cardinality -= 1
        else:
            cardinality -= 1
            add_chain(streptavidins, c)
    return streptavidins

# ------------------------------------------------------------------------------------
# other ops


def combine_and_relax(pose: pyrosetta.Pose, others, cycles=1):
    add_chain(pose, others)
    scorefxn: pr_scoring.ScoreFunction = pyrosetta.get_fa_scorefxn()
    dG_others = scorefxn(others)
    relax_chainA(pose, cycles)
    return scorefxn(pose) - dG_others

def relax_chainA(pose: pyrosetta.Pose, cycles=1, distance=5, scorefxn=None):
    if scorefxn is None:
        scorefxn: pr_scoring.ScoreFunction = pyrosetta.get_fa_scorefxn()
        scorefxn.set_weight(pr_scoring.ScoreType.atom_pair_constraint, 3)
        scorefxn.set_weight(pr_scoring.ScoreType.coordinate_constraint, 5)
    relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, cycles)
    movemap = pyrosetta.MoveMap()
    rs: ModuleType = prc.select.residue_selector
    chainA_sele: rs.ResidueSelector = rs.ChainSelector('A')
    chainA: pru.vector1_bool = chainA_sele.apply(pose)
    if distance > 0:
        neigh_sele: rs.ResidueSelector = rs.NeighborhoodResidueSelector(chainA_sele, True, distance)
        neighs: pru.vector1_bool = neigh_sele.apply(pose)
        movemap.set_chi(neighs)
    else:
        movemap.set_chi(chainA)
    movemap.set_bb(chainA)
    movemap.set_jump(False)
    relax.set_movemap(movemap)
    relax.apply(pose)
    return scorefxn(pose)


def thread(template_block, target_seq, target_name, template_name,
           temp_folder='/data/outerhome/tmp'):
    # load template
    template = pyrosetta.Pose()
    prc.import_pose.pose_from_pdbstring(template, template_block)
    # thread
    aln_filename = f'{temp_folder}/{template_name}-{target_name}.grishin'
    ph.write_grishin(target_name=target_name,
                     target_sequence=target_seq,
                     template_name=template_name,
                     template_sequence=template.sequence(),
                     outfile=aln_filename
                     )
    aln: prc.sequence.SequenceAlignment = prc.sequence.read_aln(format='grishin', filename=aln_filename)[1]
    threaded: pyrosetta.Pose
    threader: prp.comparative_modeling.ThreadingMover
    threadites: pru.vector1_bool
    threaded, threader, threadites = ph.thread(target_sequence=target_seq,
                                               template_pose=template,
                                               target_name=target_name,
                                               template_name=template_name,
                                               align=aln
                                               )
    # no need to superpose. It is already aligned
    # ...
    # fix pdb info
    n = threaded.total_residue()
    pi = prc.pose.PDBInfo(n)
    for i in range(1, n + 1):
        pi.number(i, i)
        pi.chain(i, 'A')
    threaded.pdb_info(pi)
    superpose(template, threaded)
    return threaded


# -------

def create_design_tf(pose:pyrosetta.Pose, design_sele: pr_res.ResidueSelector, distance:int) -> prc.pack.task.TaskFactory:
    """
    Create design task factory for relax.
    Designs the ``design_sele`` and repacks around ``distance`` of it.

    Remember to do

    relax.set_enable_design(True)
    relax.set_task_factory(task_factory)
    """
    #residues_to_design = design_sele.apply(pose)
    # this is default:
    # design_ops = prc.pack.task.operation.OperateOnResidueSubset(????, residues_to_design)
    no_cys = pru.vector1_std_string(1)
    no_cys[1] = 'CYS'
    no_cys_ops =  prc.pack.task.operation.ProhibitSpecifiedBaseResidueTypes(no_cys)
    # No design, but repack
    repack_sele = pr_res.NeighborhoodResidueSelector(design_sele, distance, False)
    residues_to_repack = repack_sele.apply(pose)
    repack_rtl = prc.pack.task.operation.RestrictToRepackingRLT()
    repack_ops = prc.pack.task.operation.OperateOnResidueSubset(repack_rtl, residues_to_repack)
    # No repack, no design
    frozen_sele = pr_res.NotResidueSelector(pr_res.OrResidueSelector(design_sele, repack_sele))
    residues_to_freeze = frozen_sele.apply(pose)
    prevent_rtl = prc.pack.task.operation.PreventRepackingRLT()
    frozen_ops = prc.pack.task.operation.OperateOnResidueSubset(prevent_rtl, residues_to_freeze)
    # pyrosetta.rosetta.core.pack.task.operation.RestrictAbsentCanonicalAASRLT
    # pyrosetta.rosetta.core.pack.task.operation.PreserveCBetaRLT
    task_factory = prc.pack.task.TaskFactory()
    task_factory.push_back(no_cys_ops)
    task_factory.push_back(repack_ops)
    task_factory.push_back(frozen_ops)
    return task_factory


def design_interface_onA(pose: pyrosetta.Pose, distance: int, cycles = 5, scorefxn=None):
    interface_onA = pr_res.AndResidueSelector(
                                                pr_res.NeighborhoodResidueSelector(pr_res.OrResidueSelector(pr_res.ChainSelector('B'),
                                                                                                            pr_res.ChainSelector('C')),
                                                                                   distance, False),
                                                pr_res.ChainSelector('A')
                                            )
    task_factory: prc.pack.task.TaskFactory = create_design_tf(pose, design_sele=interface_onA, distance=distance)
    if scorefxn is None:
        scorefxn = pyrosetta.get_fa_scorefxn()
    relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, cycles)
    relax.set_enable_design(True)
    relax.set_task_factory(task_factory)
    relax.apply(pose)

def design_different(pose: pyrosetta.Pose, ref: pyrosetta.Pose, cycles = 5, scorefxn=None):
    ref = ref.split_by_chain(1)
    ref2pose: dict = align_for_atom_map(pose.split_by_chain(1), ref)
    conserved = list(ref2pose.values())
    idx_sele = pr_res.ResidueIndexSelector()
    for i in range(1, len(pose.chain_sequence(1))):
        if i not in conserved:
            idx_sele.append_index(i)
    print(idx_sele.apply(pose))
    task_factory: prc.pack.task.TaskFactory = create_design_tf(pose, design_sele=idx_sele, distance=0)
    if scorefxn is None:
        scorefxn = pyrosetta.get_fa_scorefxn()
    relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, cycles)
    relax.set_enable_design(True)
    relax.set_task_factory(task_factory)
    relax.apply(pose)

def score_interface(complex: Union[pyrosetta.Pose, Sequence[pyrosetta.Pose]], interface: str):
    if isinstance(complex, Sequence):
        _complex = complex[0].clone()
        for c in complex[1:]:
            add_chain(_complex, c)
        complex = _complex
    ia = pyrosetta.rosetta.protocols.analysis.InterfaceAnalyzerMover(interface)
    ia.apply(complex)
    return {'complex_energy': ia.get_complex_energy(),
            'separated_interface_energy': ia.get_separated_interface_energy(),
            'complexed_sasa': ia.get_complexed_sasa(),
            'crossterm_interface_energy': ia.get_crossterm_interface_energy(),
            'interface_dG': ia.get_interface_dG(),
            'interface_delta_sasa': ia.get_interface_delta_sasa()}


def rosetta_pdb_to_df(pdbblock: str) -> pd.DataFrame:
    parsable = False
    _data = []
    for line in pdbblock.split('\n'):
        if '#BEGIN_POSE_ENERGIES_TABLE' in line:
            parsable = True
            continue
        elif '#END_POSE_ENERGIES_TABLE' in line:
            break
        elif not parsable:
            continue
        parts = line.strip().split()
        if parts[0] == 'label':
            _data.append(parts)
        elif parts[0] == 'weights':
            _data.append([parts[0]] + list(map(float, parts[1:-1])) + [float('nan')])
        else:
            _data.append([parts[0]] + list(map(float, parts[1:])))
    data = pd.DataFrame(_data)
    data.columns = data.iloc[0]
    data = data.iloc[1:].copy()
    return data

# -------- constraints
def constrain_chainbreak(pose, chain_break, x0_in=1.334, sd_in=0.2, tol_in=0.02):
    AtomPairConstraint = pr_scoring.constraints.AtomPairConstraint  # noqa
    fore_c = pyrosetta.AtomID(atomno_in=pose.residue(chain_break).atom_index('C'),
                                rsd_in=chain_break)
    aft_n = pyrosetta.AtomID(atomno_in=pose.residue(chain_break + 1).atom_index('N'),
                              rsd_in=chain_break + 1)
    fun = pr_scoring.func.FlatHarmonicFunc(x0_in=x0_in, sd_in=sd_in, tol_in=tol_in)
    con = AtomPairConstraint(fore_c, aft_n, fun)
    pose.add_constraint(con)

def freeze_atom(pose: pyrosetta.Pose, frozen_index: int, ref_index: int, x0_in=0., sd_in=0.01):
    ref_ca = pyrosetta.AtomID(atomno_in=pose.residue(ref_index).atom_index('CA'), rsd_in=ref_index)
    frozen_ca = pyrosetta.AtomID(atomno_in=pose.residue(frozen_index).atom_index('CA'), rsd_in=frozen_index)
    frozen_xyz = pose.residue(frozen_index).xyz(frozen_ca.atomno())
    fun = pr_scoring.func.HarmonicFunc(x0_in=x0_in, sd_in=sd_in)
    con = pr_scoring.constraints.CoordinateConstraint(a1=frozen_ca, fixed_atom_in=ref_ca, xyz_target_in=frozen_xyz, func=fun, scotype=pr_scoring.ScoreType.coordinate_constraint)
    pose.add_constraint(con)

# -------- other issues
def extract_coords(pose: pyrosetta.Pose) -> np.ndarray:
    # this seems to be present in the docs but not in my version?
    # pyrosetta.toolbox.extract_coords_pose.pose_coords_as_row
    return np.array([list(pose.xyz(pyrosetta.AtomID(a, r))) for r in range(1, 1+pose.total_residue()) for a in range(1, 1+pose.residue(r).natoms())])


def get_pose_break(pose: pyrosetta.Pose) -> int:
    """
    If the job failed, the pose will have a large distance between C and N atoms.
    Only first chain is considered.
    Returns the Fortran-style residue numbering of the first break found.
    Else returns 0

    :param pose:
    :return: index of the first break or zero
    """
    for i in range(1, pose.split_by_chain(1).total_residue() - 1):
        n_xyz: prn.xyzVector_double_t = pose.residue(i).xyz('C')
        c_xyz: prn.xyzVector_double_t = pose.residue(i+1).xyz('N')
        d: float = n_xyz.distance(c_xyz)
        if d > 1.5:
            return i
    else:
        return 0


def extract_not_chainA(pose: pyrosetta.Pose) -> pyrosetta.Pose:
    others: pyrosetta.Pose
    for ci, c in enumerate(pose.split_by_chain()):
        if ci == 0:
            continue
        elif ci == 1:
            others = c
        else:
            add_chain(others, c)
    return others

def transpant_CTD(target: pyrosetta.Pose, addendum: pyrosetta.Pose, start_seq='FKDETET') -> pyrosetta.Pose:
    addendum = addendum.split_by_chain(1)
    pyrosetta.rosetta.protocols.grafting.delete_region(addendum, 1, addendum.sequence().find(start_seq) + 1)
    junction_idx = target.sequence().find(start_seq)
    collage = pyrosetta.rosetta.protocols.grafting.return_region(target, 1, junction_idx + 1)
    pyrosetta.rosetta.core.pose.append_pose_to_pose(collage, addendum, False)
    # PDB needs resetting
    pi = collage.pdb_info()
    for r in range(1, collage.total_residue()+1):
        pi.chain(r, 'A')
        pi.number(r, r)
    return collage

def get_chainA_per_residue(pose: pyrosetta.Pose) -> List[float]:
    """
    Get the score per residue...
    This is not fixed for Hbond halving, but that is actually good in this case as it is a negative filter
    """
    scorefxn: pyrosetta.ScoreFunction = pyrosetta.get_fa_scorefxn()
    chainA: pyrosetta.Pose = pose.split_by_chain(1)
    residue_scores: List[float] = []
    for i in ph.pose_range(chainA):
        v = pru.vector1_bool(pose.total_residue())
        v[i] = 1
        residue_scores.append(scorefxn.get_sub_score(pose, v))
    return residue_scores



def movement(original: pyrosetta.Pose,
             trials: int = 100, temperature: int = 1.0, replicate_number: int = 20) -> List[float]:
    """
    This method is adapted from my one in pyrosetta-help
    """
    scorefxn = pyrosetta.get_fa_scorefxn()
    backrub = pyrosetta.rosetta.protocols.backrub.BackrubMover()
    monégasque = pyrosetta.rosetta.protocols.monte_carlo.GenericMonteCarloMover(maxtrials=trials,
                                                                                max_accepted_trials=trials,
                                                                                # gen.max_accepted_trials() = 0
                                                                                task_scaling=5,
                                                                                # gen.task_scaling()
                                                                                mover=backrub,
                                                                                temperature=temperature,
                                                                                sample_type='low',
                                                                                drift=True)
    monégasque.set_scorefxn(scorefxn)
    # find most deviant
    rs = []
    for i in range(replicate_number):
        variant = original.clone()
        monégasque.apply(variant)
        if monégasque.accept_counter() > 0:
            rs.append(pr_scoring.bb_rmsd(variant, original))
        else:
            rs.append(float('nan'))
    return rs



# ------------------------------------------------------------------------------------
# Metadata test

def ammend_pdbblock(pdbblock, info) -> str:
    """
    An attempt at storing metadata in the PDB

    :param pdbblock:
    :param info:
    :return:
    """
    lines = pdbblock.split('\n')
    lines.insert(1, f'TITLE     {info["name"]}')
    lines.insert(3, f'REMARK   2    group: {info["group"]}')
    lines.insert(4, f'REMARK   3 subgroup: {info["subgroup"]}')
    i = 5
    while True:
        if 'ATOM' in lines[i]:
            break
        i += 1
    for key in ('complex_dG', 'chainA_dG', 'other_dG',):
        lines.insert(i, f'REMARK 250 {key: >12}: {info[key]:.1f} kcal/mol')
        i += 1
    lines.insert(i, f'REMARK 250  max_residue: {max(info["per_residue"]):.1f} kcal/mol')
    lines.insert(i + 1, f'REMARK 250  wo_strep: {info["wo_strep"]:.1f} kcal/mol')
    return '\n'.join(lines)


def read_metadata(pdbblock: str) -> dict:
    if 'TITLE ' not in pdbblock:
        return {}
    info = {}
    info['name'] = re.search(r'TITLE\s+(.*)', pdbblock).group(1).strip()
    for k, v in re.findall(r'REMARK\s+\d+\s+(\w+):\s+([\d.\-+]+)', pdbblock):
        info[k] = float(v)
    return info

# ------------------------------------------------------------------------------------


def get_pdbblock(name):
    """
    A problem arising from iterations...

    :param name:
    :return:
    """
    path = Path(f'output_MPNN/pdbs/{name}.pdb')
    path2 = Path(f'output_MPNN2/pdbs/{name}.pdb')
    if path.exists():
        return path.read_text()
    elif path2.exists():
        return path2.read_text()
    else:
        raise Exception(name)
