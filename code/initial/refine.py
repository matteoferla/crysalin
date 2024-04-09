import pandas as pd
import shutil
import re
import json
import os
import pandas as pd
from pathlib import Path
from typing import Dict, List, Any, Sequence
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
pr_res: ModuleType = pyrosetta.rosetta.core.select.residue_selector
from common_pyrosetta import init_pyrosetta, design_interface_onA



# --------------

init_pyrosetta()
pose: pyrosetta.Pose = pyrosetta.pose_from_file('pentakaihemimer.relax.pdb')
design_interface_onA(pose, 3)