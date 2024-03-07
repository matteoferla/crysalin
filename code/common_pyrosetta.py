import re
from pathlib import Path
from types import ModuleType
from typing import List, Dict, Union, Optional, Sequence

import pandas as pd
import pyrosetta
import pyrosetta_help as ph
from Bio.Align import PairwiseAligner, Alignment
from Bio.SeqUtils import ProtParam

prc: ModuleType = pyrosetta.rosetta.core
prp: ModuleType = pyrosetta.rosetta.protocols
pru: ModuleType = pyrosetta.rosetta.utility
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

def get_pose_break(pose) -> int:
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


def add_chain(built: pyrosetta.Pose, new: pyrosetta.Pose, reset: bool = False) -> None:
    """
    Add a chain ``new`` to a pose ``built`` preserving the residue numbering.
    """
    built_pi = built.pdb_info()
    if built_pi is None or reset:
        built_pi = prc.pose.PDBInfo(built)
        built.pdb_info(built_pi)
        for r in range(1, built.total_residue() + 1):
            built_pi.set_resinfo(res=r, chain_id='A', pdb_res=r)
    built_pi.obsolete(False)
    for chain in new.split_by_chain():
        offset: int = built.total_residue()
        pyrosetta.rosetta.core.pose.append_pose_to_pose(built, chain, new_chain=True)
        chain_pi = chain.pdb_info()
        for r in range(1, chain.total_residue() + 1):
            built_pi.set_resinfo(res=r + offset, chain_id=chain_pi.chain(r), pdb_res=chain_pi.number(r))

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
        pi.set_chain(r, chain_letter)
        pi.set_number(r, r - previous_max)
    pi.obsolete(False)

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


def get_chainA_per_residue(pose: pyrosetta.Pose) -> List[float]:
    """
    Get the score per residue...
    This is not fixed for Hbond halving, but that is actually good
    """
    scorefxn: pyrosetta.ScoreFunction = pyrosetta.get_fa_scorefxn()
    chainA: pyrosetta.Pose = pose.split_by_chain(1)
    residue_scores: List[float] = []
    for i in ph.pose_range(chainA):
        v = pru.vector1_bool(pose.total_residue())
        v[i] = 1
        residue_scores.append(scorefxn.get_sub_score(pose, v))
    return residue_scores


def ammend_pdbblock(pdbblock, info) -> str:
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


def get_pdbblock(name):
    path = Path(f'output_MPNN/pdbs/{name}.pdb')
    path2 = Path(f'output_MPNN2/pdbs/{name}.pdb')
    if path.exists():
        return path.read_text()
    elif path2.exists():
        return path2.read_text()
    else:
        raise Exception(name)


def align_pose_to_ref_by_chain(pose, ref, chain: str) -> float:
    atom_map = pyrosetta.rosetta.std.map_core_id_AtomID_core_id_AtomID()
    chain_sele: pr_res.ResidueSelector = pr_res.ChainSelector(chain)
    for r, m in zip(pr_res.selection_positions(chain_sele.apply(ref)),
                    pr_res.selection_positions(chain_sele.apply(pose))
                    ):
        ref_atom = pyrosetta.AtomID(ref.residue(r).atom_index("CA"), r)
        mobile_atom = pyrosetta.AtomID(pose.residue(m).atom_index("CA"), m)
        atom_map[mobile_atom] = ref_atom
    return prc.scoring.superimpose_pose(mod_pose=pose, ref_pose=ref, atom_map=atom_map)


def extract_wo_strep(pose):
    chains = pose.split_by_chain()
    wo = chains[1]
    for i, c in enumerate(chains):
        if i == 0:
            continue
        if 'MEAGIT' in c.sequence():
            continue
        add_chain(wo, c)
    return wo


def extract_streps(pose):
    """
    This is just for the ref value: -318.5261163127267 (not minimized)
    """
    chains = pose.split_by_chain()
    w = chains[2]
    for i, c in enumerate(chains):
        if i == 1:
            continue
        if 'MEAGIT' not in c.sequence():
            continue
        add_chain(w, c)
    return w


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
        for chain in new.split_by_chain():
            pyrosetta.rosetta.core.pose.append_pose_to_pose(built, chain, new_chain=True)
            built_pi = built.pdb_info()
            chain_pi = chain.pdb_info()
            for r in range(1, chain.total_residue() + 1):
                offset: int = built.total_residue()
                built_pi.set_resinfo(res=r + offset, chain_id=chain_pi.chain(r), pdb_res=chain_pi.number(r))
        built_pi.obsolete(False)

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



def superpose(ref: pyrosetta.Pose, mobile: pyrosetta.Pose, aln_map: Optional[Dict[int, int]] = None) -> float:
    if aln_map is None:
        aln_map = dict(zip(range(1, ref.total_residue() + 1), range(1, mobile.total_residue() + 1)))
    # ## make pyrosetta map
    atom_map = prs.map_core_id_AtomID_core_id_AtomID()
    for r, m in aln_map.items():
        ref_atom = pyrosetta.AtomID(ref.residue(r + 1).atom_index("CA"), r + 1)
        mobile_atom = pyrosetta.AtomID(mobile.residue(m + 1).atom_index("CA"), m + 1)
        atom_map[mobile_atom] = ref_atom
    # return RMSD
    return prc.scoring.superimpose_pose(mod_pose=mobile, ref_pose=ref, atom_map=atom_map)


def relax(pose: pyrosetta.Pose, others, cycles=1):
    add_chain(pose, others)
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


def thread(template_block, target_seq, target_name, template_name, temp_folder='/data/outerhome/tmp'):
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
    aln: prc.sequence.SequenceAlignment = \
    pyrosetta.rosetta.core.sequence.read_aln(format='grishin', filename=aln_filename)[1]
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
    residues_to_design = design_sele.apply(pose)
    # this is default:
    # design_ops = prc.pack.task.operation.OperateOnResidueSubset(????, residues_to_design)
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
    task_factory.push_back(frozen_ops)
    task_factory.push_back(repack_ops)
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

def score_interface(complex: pyrosetta.Pose, interface: str):
    if isinstance(complex, Sequence):
        complex = complex[0].clone()
        for c in complex[1:]:
            add_chain(complex, c)
    ia = pyrosetta.rosetta.protocols.analysis.InterfaceAnalyzerMover(interface)
    ia.apply(complex)
    return {'complex_energy': ia.get_complex_energy(),
            'separated_interface_energy': ia.get_separated_interface_energy(),
            'complexed_sasa': ia.get_complexed_sasa(),
            'crossterm_interface_energy': ia.get_crossterm_interface_energy(),
            'interface_dG': ia.get_interface_dG(),
            'interface_delta_sasa': ia.get_interface_delta_sasa()}