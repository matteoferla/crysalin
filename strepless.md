```python
import pandas as pd

scorefxn: pr_scoring.ScoreFunction = pyrosetta.get_fa_scorefxn()
pose: pyrosetta.Pose = pyrosetta.pose_from_file('pentakaihemimer.renumbered.pdb')
print(scorefxn(pose))
energies: npt.NDArray[np.float64] = pose.energies().residue_total_energies_array()
df = pd.DataFrame(energies)
df.iloc[143-5:155]
```

```python
print('WT', scorefxn(pose))
pi: prc.pose.PDBInfo = pose.pdb_info()
ri = pi.pdb2pose(chain='A', res=150)
assert pose.residue(ri).name3() == 'TRP'
for i in range(3):
    prp.simple_moves.MutateResidue(target=ri+i, new_res='GLY').apply(pose)
print('W150G_S151G_H152G mutant', scorefxn(pose))
strep_sele: pr_res.ResidueSelector = pr_res.ResidueSpanSelector(pi.pdb2pose(chain='A', res=148), pi.pdb2pose(chain='A', res=154))
neigh_sele: pr_res.ResidueSelector = pr_res.NeighborhoodResidueSelector(strep_sele, True, 5)
neighs: pru.vector1_bool = neigh_sele.apply(pose)
relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, 15)
movemap = pyrosetta.MoveMap()
movemap.set_bb(strep_sele.apply(pose))
movemap.set_chi(neigh_sele.apply(pose))
relax.set_movemap(movemap)
relax.apply(pose)

pose.energies().clear_energies()
print(scorefxn(pose))
energies: npt.NDArray[np.float64] = pose.energies().residue_total_energies_array()
set_bfactor(pose, energies['total_score'])
pose.dump_pdb('pentakaihemimer.W150G_S151G_H152G.pdb')
```