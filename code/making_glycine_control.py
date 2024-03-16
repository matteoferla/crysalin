"""
This snippet was used to make the glycine control. Ie. pentakaihemimer.relax.pdb but with polyglycine NTD.
The exception is the streptag sequence, which is kept as is.
"""

# boilerplate imports
from pathlib import Path
import pyrosetta
import pyrosetta_help as ph
from types import ModuleType
from collections import Counter
from typing import Dict
prc: ModuleType = pyrosetta.rosetta.core
prp: ModuleType = pyrosetta.rosetta.protocols
pru: ModuleType = pyrosetta.rosetta.utility  # noqa
prn: ModuleType = pyrosetta.rosetta.numeric
prs: ModuleType = pyrosetta.rosetta.std  # noqa
pr_conf: ModuleType = pyrosetta.rosetta.core.conformation
pr_scoring: ModuleType = pyrosetta.rosetta.core.scoring
pr_options: ModuleType = pyrosetta.rosetta.basic.options
pr_res: ModuleType = pyrosetta.rosetta.core.select.residue_selector


logger = ph.configure_logger()
pyrosetta.distributed.maybe_init(extra_options=ph.make_option_string(no_optH=False,
                                                                     ex1=None,
                                                                     ex2=None,
                                                                     # mute='all',
                                                                     ignore_unrecognized_res=True,
                                                                     load_PDB_components=False,
                                                                     ignore_waters=True)
                                 )

# ---------------------

path = Path('pentakaihemimer.relax.pdb')
ctd_start_seq = 'KDETET'
strep_tag_seq = 'HPQFEKR'

pose = pyrosetta.pose_from_file(path.as_posix())
junction_idx = pose.sequence().find(ctd_start_seq) + 1  # this is the start of the CTD

streptag_idx = pose.sequence().find(strep_tag_seq)
prp.simple_moves.MutateResidue(1, 'GLY:NtermProteinFull').apply(pose)
for i in range(2, junction_idx):
    if i < streptag_idx or i> streptag_idx + len(strep_tag_seq):
        prp.simple_moves.MutateResidue(i, 'GLY').apply(pose)


pose.dump_pdb('glycine_control.pdb')