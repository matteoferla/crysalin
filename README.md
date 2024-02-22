# crysalin
Engineering crysalin lattice

NB. Many of the scripts were run in a reverse–port-forwarded Jupyter Lab notebook running in the cluster
([tutorial](https://www.blopig.com/blog/2023/10/ssh-the-boss-fight-level-jupyter-notebooks-from-compute-nodes/)).
And for archiving they are presented here as scripts with variables hardcoded in the top...
A few step generate a godzillion files in a single directory, which is bad for the FS.

## Scheme

* Experiments were given Greek letters
* Replicates of RFdiffusion were given numbers
* Replicates of ProteinMPNN were given letters, with the letter Ø for the original polyglycine (control)

However, as the experiments had minor variations, suffices to the greek letters were used,
making the end result a bit of a mess and like this `eta_comboplus_104F.pdb`.

Steps:

1. RFdiffusion
2. Filtering in PyRosetta for clashes in full complex
3. ProteinMPNN
4. Threading & relaxing in PyRosetta and filtering for clashes
5. Ranking
6. AlphaFold2 validation
7. Thermal tempering validation

## Installation

> [install.sh](code/install.sh) :construction: :warning: TODO FINISH COPYING

Installation of RFdiffusion is curious, but straightforward albeit manual ([install.sh](code/install.sh))
—it needs some cleaning up lest one is fine with stuff dumped around.

## RFdiffusion

RFdiffusion has a lot of settings. [rfdiffusion](code/job_RFdiffusion.sh) was run in the cluster
with the variable `APPTAINERENV_EXPERIMENT` controlling the experiment,
see script for each one. 

> :construction: :warning: TODO add a table of experiments and their settings.

> :construction: :warning: TODO add a picture

## Rosetta Remodel

In Parallel Rosetta Remodel was used in PyRosetta in [remodel.py](code/remodel.py).
For the sake of sanity, the constructs were renumbered PDB -> Pose numbering
as Remodel is odd when multichain PDBs and PDB numbering are used.

The blueprints are in misc, but were generated dynamically in the script 
(cf [tutorial](https://blog.matteoferla.com/2021/04/remodel-in-pyrosetta.html) for working).

Three sub-experiments were run:

* ins36_41 — 6 any SS residues with redesigned flanks (with helical preference)
* ins59_61 — 6 helical residue insertion and 3 any SS residues with redesigned flanks
* ins139_141 — 4 residue insertion with redesigned flanks

> :construction: :warning: TODO add a picture

## Filtering

See [filtering](code/filter.py) for the filtering scripts

## ProteinMPNN

The helper scripts for ProteinMPNN do not allow variable length protein and are a bit awkward to use.
The [functional_proteinMPNN_helper.py](code/functional_proteinMPNN_helper.py) is a refactored version called by 
[prep_MPNN.py](code/prep_MPNN.py).

## Threading & relaxing

In [thread.py](code/thread.py) I used the threading module from RosettaCM in PyRosetta as opposed to 
brutally substituting each residue sequentially as in original paper.
I then placed the monomer in a subset complex and relaxed the chain.

> :construction: :warning: TODO graphs galore



