"""
$INSERTION_LEN
$EXPERIMENT ins139_141 ins59_61 ins36_41
"""


# check input is correct
import os

exp_names = ['ins139_141', 'ins59_61', 'ins36_41']
exp_name = os.environ.get('EXPERIMENT', '')
if exp_name not in exp_names:
    raise ValueError(f'{exp_name} not in {exp_names}')
else:
    print(f'Running {exp_name}')

# --------------------------------------------
# ## Hardcoded settings
n_trials = 1_000
quick_and_dirty = False
generic_aa = 'G'
dr_cycles = 3  # default 3
max_linear_chainbreak = 0.2  # default 0.07
redesign_loop_neighborhood = False
outfolder = 'output_remodel'
timeout = 60 * 60
n_insertions = int(os.environ.get('INSERTION_LEN', 4))

os.makedirs(outfolder, exist_ok=True)
os.makedirs(f'{outfolder}/polyG', exist_ok=True)
os.makedirs(f'{outfolder}/broken', exist_ok=True)
os.makedirs(f'{outfolder}/valid', exist_ok=True)

# --------------------------------------------

from pathlib import Path
import pyrosetta
import pyrosetta_help as ph
from types import ModuleType
from typing import Generator
from collections import Counter
import pebble
import multiprocessing
from concurrent.futures import TimeoutError
from common_pyrosetta import get_pose_break, align_pose_to_ref_by_chain, fix_pdb, score_interface

prc: ModuleType = pyrosetta.rosetta.core
prp: ModuleType = pyrosetta.rosetta.protocols
pru: ModuleType = pyrosetta.rosetta.utility
prn: ModuleType = pyrosetta.rosetta.numeric
pr_conf: ModuleType = pyrosetta.rosetta.core.conformation
pr_scoring: ModuleType = pyrosetta.rosetta.core.scoring
pr_options: ModuleType = pyrosetta.rosetta.basic.options

logger = ph.configure_logger()
pyrosetta.distributed.maybe_init(extra_options=ph.make_option_string(no_optH=False,
                                                                     ex1=None,
                                                                     ex2=None,
                                                                     # mute='all',
                                                                     ignore_unrecognized_res=True,
                                                                     load_PDB_components=False,
                                                                     ignore_waters=True)
                                 )



# load pose
original: pyrosetta.Pose = pyrosetta.pose_from_file('pentakaihemimer.relax.pdb')  # trimer-renumbered.pdb
blue = ph.Blueprinter.from_pose(original)

# make specific bluprint
if exp_name == 'ins139_141':
    blue[139:141] = 'NATAA'
    for i in range(n_insertions):
        blue.insert_after(139, 'NATAA')
elif exp_name == 'ins59_61':
    blue[59:61] = 'NATAA'
    blue.rows[59 - 1][2] = 'H'
    blue.pick_native(61)
    # extra turn
    for i in range(3):
        blue.insert_before(60, 'NATAA')
        blue.rows[60 + i - 1][2] = 'H'
    # return
    for i in range(n_insertions):
        blue.insert_before(60, 'NATAA')
elif exp_name == 'ins36_41':
    blue[36:42] = 'NATAA'
    blue.rows[36 - 1][2] = 'H'
    blue.rows[37 - 1][2] = 'H'
    for i in range(n_insertions):
        blue.insert_after(37, 'NATAA')
else:
    raise ValueError(f'{exp_name} unknown')

blue.set(f'{exp_name}-{n_insertions}.blu')
pr_options.set_boolean_option('remodel:quick_and_dirty', quick_and_dirty)
pr_options.set_string_option('remodel:generic_aa', generic_aa)
rm = prp.forge.remodel.RemodelMover()
rm.register_options()
rm.dr_cycles(dr_cycles)
rm.max_linear_chainbreak(max_linear_chainbreak)
rm.redesign_loop_neighborhood(redesign_loop_neighborhood)


# run
def run(i):
    scorefxn: pr_scoring.ScoreFunction = pyrosetta.get_fa_scorefxn()
    run_name = f'{exp_name}-{n_insertions}_{i:0>3}'
    pose: pyrosetta.Pose = original.clone()
    rm.apply(pose)
    print(f'remodel #{i} complete')
    if 'GGG' in pose.sequence():
        pose.dump_pdb(f'{outfolder}/polyG/{run_name}.pdb')
        raise ValueError('GGG found')
    chain_break = get_pose_break(pose)
    if chain_break:
        pose.dump_pdb(f'{outfolder}/broken/{run_name}.pdb')
    else:
        pose.dump_pdb(f'{outfolder}/valid/{run_name}.pdb')
    # deal with chain-break
    if chain_break:
        pose.dump_pdb(f'{outfolder}/broken/{run_name}.pdb')
        AtomPairConstraint = pr_scoring.constraints.AtomPairConstraint  # noqa
        fore_c = pyrosetta.AtomID(atomno_in=pose.residue(chain_break).atom_index('C'),
                                    rsd_in=chain_break)
        aft_n = pyrosetta.AtomID(atomno_in=pose.residue(chain_break + 1).atom_index('N'),
                                  rsd_in=chain_break + 1)
        fun = pr_scoring.func.HarmonicFun(x0_in=1.334, sd_in=0.2)
        con = AtomPairConstraint(fore_c, aft_n, fun)
        pose.add_constraint(con)
        scorefxn.set_weight(pr_scoring.ScoreType.atom_pair_constraint, 3)
    # ## Relax
    relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, 3)
    movemap = pyrosetta.MoveMap()
    chain_A = prc.select.residue_selector.ChainSelector(1).apply(pose)
    movemap.set_bb(chain_A)
    movemap.set_chi(chain_A)
    relax.set_movemap(movemap)
    relax.apply(pose)
    pose.dump_pdb(f'{outfolder}/valid/{run_name}.pdb')
    fix_pdb(pose)
    pose.pdb_info().obsolete(False)
    align_pose_to_ref_by_chain(pose, original, 'B')
    monomer = pose.split_by_chain(1)
    monomer.dump_pdb(f'{outfolder}/relaxed/{run_name}.pdb')
    holo = scorefxn(pose)
    score_interface(pose, 'A_BC')
    apo = sum(map(scorefxn, pose.split_by_chain()))
    return dict(run_name=run_name, seq=pose.sequence(), holo_score=holo, apo_scores=apo, delta_score=holo - apo,
                error='')


def run_done(future):
    nan = float('nan')
    try:
        result = future.result()  # blocks until results are ready
    except Exception as error:
        error_msg = error.args[1]
        if isinstance(error, TimeoutError):
            print("Function took longer than %d seconds" % error_msg)
        else:
            print("Function raised %s" % error)
            print(error.__traceback__)  # traceback of the function
        result = dict(seq='',
                      holo_score=nan, apo_scores=nan, delta_score=nan,
                      error=error.__class__.__name__, error_msg=error_msg)
    with open(f'{outfolder}/{exp_name}.tsv', 'a') as fh:
        fh.write('\t'.join(result.values()) + '\n')


num_cores = multiprocessing.cpu_count()

with pebble.ProcessPool(max_workers=num_cores - 1, max_tasks=0) as pool:
    futuremap: pebble.ProcessMapFuture = pool.map(run, range(1, n_trials), timeout=timeout)
    iter_future: Generator = futuremap.result()
    results = []
    while True:
        try:
            result = next(iter_future)
            results.append(result)
            print(result['run_name'])
        except StopIteration:
            break
        except TimeoutError as error:
            print(error.__class__.__name__, str(error))
        except Exception as error:
            print(error.__class__.__name__, str(error))
