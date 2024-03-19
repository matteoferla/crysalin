# -*- coding: utf-8 -*-

import operator
import string
from Bio import SeqIO
from pathlib import Path
import re, os
from Bio.SeqUtils import ProtParam
from pebble import ProcessPool, ProcessFuture
from concurrent.futures import TimeoutError
import time
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


template_path = Path('output_redux/filtered')
seq_paths = list(Path('output_redux/seqs').glob('*.fa'))
out_path = Path('output_redux/pdbs')
csv_path = Path('threading_redux.csv')
apo_score = -337.2997472992847
split_size = 250

# --------------------------------
print('\n## Init PyRosetta\n')

logger = ph.configure_logger()
pyrosetta.distributed.maybe_init(extra_options=ph.make_option_string(no_optH=False,
                                                                     ex1=None,
                                                                     ex2=None,
                                                                     # mute='all',
                                                                     ignore_unrecognized_res=True,
                                                                     load_PDB_components=False,
                                                                     ignore_waters=True,
                                                                    )
                                 )
pr_options.set_boolean_option('in:detect_disulf', False)

# --------------------------------
# Define functions

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


def run_process(target_name, **kvargs):
    tick = time.time()
    (out_path / f'{target_name}.dummy').write_text('hello world')
    threaded: pyrosetta.Pose = thread(target_name=target_name, **kvargs)
    add_chain(pose, others)
    dG = relax(threaded, others)
    threaded.dump_pdb(f'{out_path}/{target_name}.pdb')
    print(f'{target_name}: {dG} kcal/mol ({tick - time.time()}s)')
    return dict(target_name=target_name, dG=dG)

# --------------------------------
print('\n## Load the data\n')

# df1 = pd.read_pickle('scored_complexes.pkl.gz').set_index('name')
# df2 = pd.read_pickle('scored_complexes2.pkl.gz') #.set_index('name')
# df = pd.concat([df2, df1])
df = pd.read_pickle(data_path)
df = df.reset_index().drop_duplicates('name', keep='first').set_index('name')
if split_size:
    df = df.sample(split_size).copy()
df['dG_bind'] = df.complex_score - df.monomer_score - (apo_score)
df['group'] = df.index.to_series().str.split('_').apply(operator.itemgetter(0))
df['subgroup'] = df.index.to_series().apply(lambda name: '_'.join(name.split('_')[:-1]))

# --------------------------------
print('\n## Running\n')

futures: List[ProcessFuture] = []
preresults: List[dict] = []
results: List[dict] = []
t = 60 * 60 * 1
# this is the 5 and a half mer without A. It is therefore a 4 & 1/2 mer
others: pyrosetta.Pose = pyrosetta.pose_from_file('woA_pentakaihemimer.pdb')

with ProcessPool(max_workers=os.cpu_count() - 1, max_tasks=0) as pool:
    # queue jobs
    n = 0
    for path in seq_paths:
        for seq_record in SeqIO.parse(path, 'fasta'):
            n += 1
            # first record is input:
            # >beta_0-checkpoint, score=2.4184, global_score=1.5778, fixed_chains=['B', 'C', 'D'], designed_chains=['A'], model_name=v_48_020, git_hash=unknown, seed=37
            # herein I will call the input sample zero and use it as a negative control
            # others are:
            # >T=0.1, sample=1, score=1.0129, global_score=1.5088, seq_recovery=0.0312
            info = {k: float(v) for k, v in re.findall(r'([\w_]+)=([\d.]+)', seq_record.description)}
            info['name']: str = os.path.splitext(path.name)[0]  # noqa
            info['sequence']: str = str(seq_record.seq) # noqa
            info['pI'] = ProtParam.ProteinAnalysis(seq_record.seq).isoelectric_point()
            preresults.append(info)
            if info['name'] not in df.index:
                with open('missing.txt', 'a') as fh:
                    fh.write(f'{info["name"]}\n')
                continue
            row: pd.Series = df.loc[info['name']]
            # remove LINK... which there should not be anyway
            template_block = re.sub(r'(LINK[^\n]+\n)', '', row.monomer)
            template_block = re.sub(r'(SSBOND[^\n]+\n)', '', template_block)
            template_name = info['name']
            target_seq = info['sequence']
            # neg control will be Ø which looks like 0̷
            # okay. Ø er ikke det sidste bogstav, Å er. But shmeh
            target_name = f"{info['name']}{'ABCDEFGHIJKLMNOPQRSTUVWXYZØ'[int(info.get('sample', -1))]}"
            if (out_path / f'{target_name}.pdb').exists() \
                or (out_path / f'{target_name}.dummy').exists() \
                or (old_path / f'{target_name}.pdb').exists():
                continue
            assert len(row.sequence) == len(info['sequence']), f'The sequences for {target_name} does not match length'
            future: ProcessFuture = pool.schedule(run_process,
                                   kwargs=dict(template_block=template_block, target_seq=target_seq,
                                               target_name=target_name, template_name=template_name),
                                   timeout=t)
    print(f'Submitted {n} processes')
    # get results
    for future in futures:
        try:
            result = future.result() # blocks until results are ready
            print(result)
            results.append(result)
        except Exception as error:
            error_msg = str(error)
            if isinstance(error, TimeoutError):
                print(f'Function took longer than {t} seconds {error}')
            else:
                print(f"Function raised {error}")
                print(error.__traceback__)  # traceback of the function
            results.append(dict(error=error.__class__.__name__, error_msg=error_msg))

df = pd.DataFrame(results)
if csv_path.exists():
    original = pd.read_csv(csv_path)
    df = pd.concat([original, df])
df.to_csv(csv_path)

