> This are sections that are being rewritten. They are not included in the final version of the notes.



## ProteinMPNN

> Script: [code/prep_MPNN.py](code/prep_MPNN_v1.py)

The helper scripts for ProteinMPNN do not allow variable length protein and are a bit awkward to use.
The [functional_proteinMPNN_helper.py](code/functional_proteinMPNN_helper.py) is a refactored version called by
[prep_MPNN.py](code/prep_MPNN_v1.py).

## Threading, relaxing & tuning

> Script: [code/thread_tune.py](code/thread_tune_v1.py)

I used the threading module from RosettaCM in PyRosetta as opposed to
brutally substituting each residue sequentially as in original paper.
I then placed the monomer in a subset complex and relaxed the chain.

> :construction: :warning: TODO graphs galore

## StrepTag

Initially, I removed the streptag from the designs, but this was determined to be unnecessary and unhelpful.

The StrepTag is an 8+ residue tag that binds to Streptavidin.
Some constructs retain the tag, others don't.
I did not remove it from the templates as I did not want its rediscovey to be a factor in the ranking.
However, it causes issues.
The sequence between 149:A and 159:A (inclusives) is `NWSHPQFEKRP`.

* Complex minimised: -1974.9 kcal/mol
* Monomer extracted: -773.9 kcal/mol
* Tag extracted: +9.5 kcal/mol
* Monomer w/o tag: -768.0 kcal/mol
* Complex w/o tag: -1949.2 kcal/mol

The tag only add -15.5 kcal/mol to the monomer's score, but -19.7 kcal/mol to the complex's score.

## AF2 validation

This is problematic.
Based on visual inspection, the following was found. AF2 predicts:

* Steptagged WT AHIR + streptavidin in single-sequence mode: AHIR is split up, binds side of streptavidin
* Steptagged WT AHIR + streptavidin in MSA mode: streptag is roughly bound right, angle is off between two protein
* WT AHIR Rossmann fold in single-sequence mode: mostly predicted correctly
* Streptavidin with streptag: in correct area, but poor

Template mode could be used, but that would be useless to filter out monomers that are broken,
if the streptag is kept in the template, then things might work for dimer inference.
However, the single streptavidin has an exposed binding interface.


## Thermal tempering

> :construction: :warning: TODO ...