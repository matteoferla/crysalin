# crysalin

Aim:  engineer a better crysalin lattice.

![interface-design](images/interface-design.jpg)

## Background

Crysalins are an engineered complex that forms a lattice for encapsulating protein for X-ray crystallography.
The complex consists of AHIR and streptavidin.
AHIR consists of N-terminal Rossmann fold knob with strep-tag and a C-terminal helical bundle that knots with others,
forming a dodecamer, resulting in a cubic lattice with each steptavidin tetramers forming the vertices,
leaving a void for cargo.

The Rossmann fold binds weakly to Streptavidin beyond the tag held via a floppy linker.
Fortuitously, the complex is stabilised by a sugar, but that is a crystallisation artefact.
Crystallography >3Å, while AHIR ~1.6Å, Streptavidin ~1.0Å.
Previous attempts were point-mutants, but need major remodelling.
Therefore, can new designs be made that bind more strongly and are more rigid?

(![parts](images/parts.jpg))

## People involved
Project members (Newcastle): Martin Noble, Mathew Martin, Rhianna Rowland
Project members (Oxford): Frank von Delft, Michael Fairhead, me (Matteo Ferla)

## Aim

Run both RFdiffusion and Rosetta Remodel to generate designs.

## Methods

> A subset of this project was written up as a blog post:
> [PyRosetta for RFdiffusion](https://www.blopig.com/blog/2024/05/pyrosetta-for-rfdiffusion/)

### Summary
Rosetta Remodel is useful for generating parametrically designed proteins.
However, very few designs were successful and were outcompeted by RFdiffusion designs which tweaked the Rossmann fold.

RFdiffusion-AA was tried, but due to its early release status the results were more problematic than RFdiffusion.

Many iterations were performed.
Whereas RFdiffusion is for de novo design of whole domains, 
tweaking the Rossmann fold (primarily resi 41-46,64-66,144-146) was far more successful.

## Steps

1. RFdiffusion
2. Superposition to template in full complex: filtering ref2015 score (initial) or by Euclidean distance
3. ProteinMPNN x5
4. Threading
5. Complex with other AHIR dimer's stalk domain
6. Relaxing only chain A
8. FastDesign on differing residues
9. Interface scoring & interactions with Streptavidin counted
9. Ranking

Initially, these three extra steps were planned:

1. AlphaFold2 validation
2. Thermal tempering validation
3. Backrub drift

But various test failed their validity.

## Iterations

Two major design cycles were done.
The first was when Iris cluster still was on CentOS 7, so Singularity needed to be run.
The second was when it was fine.

## Template
The complex is huge.
The final template used in diffusion was the pentakaihemimer (5 and a half chains: (1.5x AHIR, 4x Streptavidin).

![pentakaihemimer](images/pentakaihemimer.png)

This is because collisions were happening, for example here are the first 5 models from experiment Alpha,
which was a dimer...

![alpha](images/alpha.png)

Another example, from theta-hot experiment, which uses the trikaihemimer (1.5x AHIR, 2x Streptavidin):
![img.png](images/thetahot.png)

## Results

### Iteration 1
> For iteration 1. See [iteration1](iteration_1.md)

### Iteration 2

> For iteration 2. See [iteration2](iteration_2.md)



## Footnotes

see [former_notes.md](former_notes.md) for more details.

### helical

RFdiffusion is said to be good at helices and sheet, but poor at loops.
Here is a helical interaction for example:

![img.png](images/helical.png)

### Open-closed loop choice

![open-closed](images/open-close.png)

PDB:1SWE vs PDB:1SWB

> :construction: :warning: TODO ...

### I remember seeing that

I have not run AF2 validation en masse, however I am concerned by some tests:
AF2 makes models which better match the template than Rosetta energy minimised.
That is to say, AF2 is remembering seeing the fold, which I don't like,
hence the threading+relaxing step
and leaving AF2 as a 'redocking test'.

![test-zero](images/test0.jpg)

### Caveat

In initial steps, many of the scripts were run in a reverse–port-forwarded Jupyter Lab notebook running in the cluster
([tutorial](https://www.blopig.com/blog/2023/10/ssh-the-boss-fight-level-jupyter-notebooks-from-compute-nodes/)).
And for archiving they are presented here as scripts with variables hardcoded in the top...
A few step generate a godzillion files in a single directory, which is bad for the FS.

### Streptavidin binders

Does Streptavin have natural binders?

Guilt by association, especially in thermophiles, is a good way to catch functional intertwined protein.
I just checked Steptomyces avidinii but not its relatives and it's not in a operon.
Its neighbours are:

* WP_189973519.1: Aryl-sulfate sulfotransferase (AssT) is a dimer in _E. coli_. In operon w/ ABC transporter
* WP_189973633.1: Cytochrome P450. Suppressed entry, but homologous to other actinomycetian P450

* Thematically similar, but not a smoking gun

![img.png](images/operon.png)


### RFdiffusion installation

> [install.sh](code/install.sh) :construction: :warning: TODO FINISH COPYING

Installation of RFdiffusion is curious, but straightforward albeit manual ([install.sh](code/install.sh))
—it needs some cleaning up lest one is fine with stuff dumped around.
