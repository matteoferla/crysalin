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

data_path = Path('scored_complexes3.pkl.gz')
seqs_path = Path('output_MPNN2/seqs')
out_path = Path('output_MPNN2/pdbs')
csv_path = Path('threading3.csv')
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
                                                                    #in=dict(in:detect_disulf=True)
                                                                    )
                                 )
pr_options.set_boolean_option('in:detect_disulf', True)

# --------------------------------
# Define functions



def run_process(target_name, **kvargs):
    tick = time.time()
    (out_path / f'{target_name}.dummy').write_text('hello world')
    threaded: pyrosetta.Pose = thread(target_name=target_name, **kvargs)
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
    for path in Path(seqs_path).glob('*.fa'):
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
            template_block = row.monomer
            template_name = info['name']
            target_seq = info['sequence']
            # neg control will be Ø which looks like 0̷
            # okay. Ø er ikke det sidste bogstav, Å er. But shmeh
            target_name = f"{info['name']}{'ABCDEFGHIJKLMNOPQRSTUVWXYZØ'[int(info.get('sample', -1))]}"
            if (out_path / f'{target_name}.pdb').exists() or (out_path / f'{target_name}.dummy'):
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

