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
prefix = 'output/remodel_'
timeout = 60 * 60

# --------------------------------------------

from pathlib import Path
import pyrosetta
import pyrosetta_help as ph
from types import ModuleType
from collections import Counter
from pebble import ProcessPool
import multiprocessing
from concurrent.futures import TimeoutError

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
original: pyrosetta.Pose = pyrosetta.pose_from_file('trimer-renumbered.pdb')
blue = ph.Blueprinter.from_pose(original)

# make specific bluprint
if exp_name == 'ins139_141':
    blue[139:141] = 'NATAA'
    for i in range(4):
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
    for i in range(6):
        blue.insert_before(60, 'NATAA')
elif exp_name == 'ins36_41':
    blue[36:42] = 'NATAA'
    blue.rows[36 - 1][2] = 'H'
    blue.rows[37 - 1][2] = 'H'
    for i in range(6):
        blue.insert_after(37, 'NATAA')
else:
    raise ValueError(f'{exp_name} unknown')

blue.set(f'{exp_name}.blu')
pr_options.set_boolean_option('remodel:quick_and_dirty', quick_and_dirty)
pr_options.set_string_option('remodel:generic_aa', generic_aa)
rm = prp.forge.remodel.RemodelMover()
rm.register_options()
rm.dr_cycles(dr_cycles)
rm.max_linear_chainbreak(max_linear_chainbreak)
rm.redesign_loop_neighborhood(redesign_loop_neighborhood)


# run
def run(i):
    run_name = f'{prefix}{exp_name}_{i:0>3}'
    pose: pyrosetta.Pose = original.clone()
    rm.apply(pose)
    pose.dump_pdb(f'{run_name}.pdb')
    print(f'remodel #{i} complete')
    if 'GGG' in pose.sequence():
        raise ValueError('GGG found')
    scorefxn: pr_scoring.ScoreFunction = pyrosetta.get_fa_scorefxn()
    relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, 3)
    relax.apply(pose)
    pose.dump_pdb(f'{run_name}.relaxed.pdb')
    holo = scorefxn(pose)
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
            print(error.traceback)  # traceback of the function
        result = dict(seq='',
                      holo_score=nan, apo_scores=nan, delta_score=nan,
                      error=error.__class__.__name__, error_msg=error_msg)
    with open(f'{prefix}{exp_name}.tsv', 'a') as fh:
        fh.write('\t'.join(result.values()) + '\n')


num_cores = multiprocessing.cpu_count()

with ProcessPool(max_workers=num_cores - 1, max_tasks=0) as pool:
    for index in range(1, n_trials):
        future = pool.schedule(run, [index], timeout=timeout)
        future.add_done_callback(run_done)
