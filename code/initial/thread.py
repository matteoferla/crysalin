# -*- coding: utf-8 -*-
import os
import operator
import string
from Bio import SeqIO
from pathlib import Path
import re, os, traceback
from Bio.SeqUtils import ProtParam
from pebble import ProcessPool, ProcessFuture
from concurrent.futures import TimeoutError
import time
import json
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
from common_pyrosetta import add_chain, superpose, constrain_chainbreak, freeze_atom, relax_chainA

work_path = Path(os.environ.get('WORKPATH', 'output_redux'))
template_folder = work_path / 'filtered'
# annoying mistake:
seqs_folder = work_path / 'seqs/seqs' if (work_path / 'seqs/seqs').exists() else work_path / 'seqs'
seq_paths = list(seqs_folder.glob('*.fa'))
raw_out_path = work_path / 'unrelaxed_pdbs'
os.makedirs(raw_out_path, exist_ok=True)
relaxed_out_path = work_path / 'relaxed_pdbs'
os.makedirs(relaxed_out_path, exist_ok=True)
i = 0
while True:
    csv_path = Path(f'threading_{i}.csv')
    if not csv_path.exists():
        break
    else:
        i+=1

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

def get_ATOM_only(pdbblock: str) -> str:
    return '\n'.join([line for line in pdbblock.splitlines() if line.startswith('ATOM')])
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
    superpose(template, threaded)
    return threaded


def run_process(target_name, **kvargs):
    tick = time.time()
    (raw_out_path / f'{target_name}.dummy').write_text('hello world')
    monomer: pyrosetta.Pose = thread(target_name=target_name, **kvargs)
    monomer.dump_pdb(f'{raw_out_path}/{target_name}.pdb')
    oligomer: pyrosetta.Pose = monomer.clone()
    add_chain(oligomer, others)  # as this is in place
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
    dG = relax_chainA(oligomer, cycles=5, distance=0, scorefxn=scorefxn)  # only chain A
    monomer = oligomer.split_by_chain(1)
    monomer.dump_pdb(f'{relaxed_out_path}/{target_name}.pdb')
    info = dict(target_name=target_name, dG=dG)
    print(f'{target_name}: {dG} kcal/mol ({tick - time.time()}s)')
    with (work_path / 'threading.jsonl').open('a') as fh:
        fh.write(f'{json.dumps(info)}\n')
    return info

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

# --------------------------------
print('\n## Running\n')

futures: List[ProcessFuture] = []
preresults: List[dict] = []
results: List[dict] = []
timeout = 60 * 60 * 5
# this is the 5 and a half mer without A. It is therefore a 4 & 1/2 mer
others: pyrosetta.Pose = pyrosetta.pose_from_file('woA_pentakaihemimer.pdb')

with ProcessPool(max_workers=os.cpu_count() - 1, max_tasks=0) as pool:
    # queue jobs
    for path in seq_paths:
        for seq_record in SeqIO.parse(path, 'fasta'):
            # first record is input:
            # >beta_0-checkpoint, score=2.4184, global_score=1.5778, fixed_chains=['B', 'C', 'D'], designed_chains=['A'], model_name=v_48_020, git_hash=unknown, seed=37
            # herein I will call the input sample zero and use it as a negative control
            # others are:
            # >T=0.1, sample=1, score=1.0129, global_score=1.5088, seq_recovery=0.0312
            info = {k: float(v) for k, v in re.findall(r'([\w_]+)=([\d.]+)', seq_record.description)}
            info['name']: str = path.stem  # noqa
            info['sequence']: str = str(seq_record.seq) # noqa
            info['pI'] = ProtParam.ProteinAnalysis(seq_record.seq).isoelectric_point()
            preresults.append(info)
            pdb_path = Path(template_folder / (path.stem + '.pdb'))
            if not pdb_path.exists():
                with open('missing_redux.txt', 'a') as fh:
                    print(f'{path.stem} is missing')
                    fh.write(f'{info["name"]}\n')
                continue
            monomer_block = pdb_path.read_text()
            # remove LINK... which there should not be anyway
            # template_block = re.sub(r'(LINK[^\n]+\n)', '', row.monomer)
            # template_block = re.sub(r'(SSBOND[^\n]+\n)', '', template_block)
            template_name = info['name']
            template_seq = get_sequence(monomer_block)
            target_seq = info['sequence'][:len(template_seq)]  # not going to thread the other chains, only A.
            #assert len(template_seq) == len(target_seq), f'The sequences for {target_name} does not match length'
            # neg control will be Ø which looks like 0̷
            # okay. Ø er ikke det sidste bogstav, Å er. But shmeh
            target_name = f"{info['name']}{'ABCDEFGHIJKLMNOPQRSTUVWXYZØ'[int(info.get('sample', -1))]}"
            if (raw_out_path / f'{target_name}.pdb').exists() \
                or (raw_out_path / f'{target_name}.dummy').exists() \
                or (relaxed_out_path / f'{target_name}.pdb').exists():
                continue
            future: ProcessFuture = pool.schedule(run_process,
                                                  kwargs=dict(template_block=monomer_block, target_seq=target_seq,
                                               target_name=target_name, template_name=template_name),
                                                  timeout=timeout)
            futures.append(future)
    print(f'Submitted {len(futures)} processes')
    # get results
    for future in futures:
        try:
            result = future.result()  # blocks until results are ready
            print(result)
            results.append(result)
        except Exception as error:
            error_msg = str(error)
            if isinstance(error, TimeoutError):
                print(f'Function took longer than {timeout} seconds {error}')
            else:
                print(f"Function raised {error}")
                traceback.print_tb(error.__traceback__)  # traceback of the function
            results.append(dict(error=error.__class__.__name__, error_msg=error_msg))

df = pd.DataFrame(results)
if csv_path.exists():
    original = pd.read_csv(csv_path)
    df = pd.concat([original, df])
df.to_csv(csv_path)

