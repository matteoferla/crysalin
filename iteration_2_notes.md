## Fixes

## Commands run

```python
import yaml

with open('experiments.yaml', 'r') as file:
    experiment_definitions = yaml.safe_load(file)
print('export NUMDESIGNS=1000;')
print('cd $HOME2/crysalin-redux')
for definition in experiment_definitions.values():
    print('export EXPERIMENT={codename}; export HOTSPOTS=\'{hotspots}\';export CONTIGMAP=\'{contig}\'; sbatch rf.slurm.sh'.format(**definition))
```
The names are words, see experiment.yaml codename, a digit for the iteration (last = 4),
underscore number (That is the RFdiffusion design) followed by letter,
where Ø is the polyglycine reference.


### Tag optimisation
The reference pose has Streptag+Streptavidin, which have very high energy.
See [first_postmortem](first_postmortem/README.md) for details.

These are possibly 'tangled' (cf. James Holton's work).
I have not tested if this is the case.






```python
# show hotspots
from collections import defaultdict

hotspots = '[A203,A204,A207,B591,B525,B13,B15,B31,B33,B34,B35,B37,B40,B42,B43,B67,B68,B69,B70,B71,B72,B73,B74,B75,B76,B78,B80,B96,B98,B100,B124,B125,B126,B146,B148,B170,B171,B173,B174,B175,B176,B177,B178,B189,B207,B209,B410,B411]'

perchain = defaultdict(list)
for hot in hotspots[1:-1].split(','):
    perchain[hot[0]].append(hot[1:])

' or '.join([f'(chain {c} and resi '+'+'.join(rs)+')' for c, rs in perchain.items()])
```

```pymol
set cartoon_gap_cutoff, 0
select hotspots, (chain A and resi 203+204+207) or (chain B and resi 591+525+13+15+31+33+34+35+37+40+42+43+67+68+69+70+71+72+73+74+75+76+78+80+96+98+100+124+125+126+146+148+170+171+173+174+175+176+177+178+189+207+209+410+411)
select bed, chain A and (resi 9-24)
select pocket, chain A and (resi 27-30 or resi 50-55 or resi 83-86 or resi 134-137 or resi 268-271)
select mid, chain A and (resi 42-43)
select side, chain A and (resi 64-66)
select inner, chain A and (resi 145-146)
select posttag, chain A and (resi 161-167)

color gray80
color turquoise, hotspots
color coral, bed or mid or side or inner or posttag
color atomic, not element C

remove chain A and resi 148-159
alter chain A, resv-=4
sort
```



## Show nice overall view

```python
from pathlib import Path
import pymol2
import re
import pickle
from typing import List, Dict, Tuple, Any
from collections import defaultdict


root = '/Users/user/Coding/crysalin'
colormap = {'bed4': '0xFF7F50', 
            'mid4': '0xFFD700', 
            'pocketfill4': '0x40E0D0',
            'posttag4': '0x800080', 
            'side4': '0x32CD32'}

hotspots = {'bed4': '[B37,B39,B41,B43,B70,B71]', 
            'mid4': '[B43,B68,B70,B75,B147,B170,B173,B174]', 
            'pocketfill4': '[A205,A325,B568,B572,B576,B579,B587,B588]',
            'posttag4': '[A135,A139,A148,A149,A151,A158,A194,A258,B37,B39,B507,B564,B568]', 
            'side4': '[B124,B125,B126,B127,B171,B173,B174,B175,B176,B177,B178,B178,B182,B207,B209]'}


with pymol2.PyMOL() as pymol:
    pymol.cmd.load(f'{root}/pentakaihemimer_renumbered.pdb', 'ref')
    pymol.cmd.color('0xDCDCDC', 'chain A')
    pymol.cmd.color('0x808080', 'not chain A')
    # add hot colors
    for group, hots in hotspots.items():
        perchain = defaultdict(list)
        for hot in hots[1:-1].split(','):
            perchain[hot[0]].append(hot[1:])
        hot_sele = ' or '.join([f'(chain {c} and resi '+'+'.join(rs)+')' for c, rs in perchain.items()])
        pymol.cmd.select(f'hot_{group}', hot_sele)
        print(f"cmd.color('{colormap[group]}', 'hot_{group} and sidechain and element C'); cmd.show('sticks', hot_{group}')")
        #pymol.cmd.color(colormap[group], f'not element H and {hot_sele} and sidechains')
        #pymol.cmd.show('sticks', f'not element H and {hot_sele}')
    # add designs
    for path in Path(f'{root}/checking').glob('*/*.pdb'):
        if int(re.search(r'_(\d+).', path.name).group(1)) >= 20:
            continue
        pymol.cmd.load(path.as_posix(), path.stem)
        pymol.cmd.align(f'{path.stem} and chain B', 'ref and chain B')
        trb_path = path.parent / (path.stem + '.trb')
        assert trb_path.exists(), f'{trb_path} does not exist'
        trb: Dict[str, Any] = pickle.load(trb_path.open('rb'))
        altered = [resi+1 for resi, is_kept in enumerate(trb['inpaint_str']) if not is_kept]
        altered_sele = "resi " + '+'.join(map(str, altered))
        pymol.cmd.remove(f'{path.stem} and not {altered_sele}')
        group = path.stem.split('_')[0]
        pymol.cmd.color(colormap[group], path.stem)
        pymol.cmd.set('cartoon_transparency', 0.5, path.stem)
        pymol.cmd.group(group, path.stem)
    pymol.cmd.zoom('ref')
    pymol.cmd.bg_color('white')
    pymol.cmd.set('cartoon_gap_cutoff', 0)
    pymol.cmd.color('atomic', 'not element C')
    pymol.cmd.save('collage_test.pse')
```

Cutout
```python


    #
    # proximals = {}
    # n_strep = 0
    # strep_chain_idxs: List[CTypeIdx] = [i for i, chain in enumerate(chains) if chain.sequence()[:6] == strep_seq]
    # pi = hallucination.pdb_info()
    # i: CTypeIdx
    # for i in range(1, hallucination.num_chains()):  # chain 0 is the full AHIR, designed
    #     chain_letter = pi.chain(hallucination.chain_begin(i + 1))
    #     xyz_strep = extract_coords(chains[i + 1])
    #     distances = np.sqrt(np.sum((xyz_strep[:, np.newaxis, :] - xyz_model[np.newaxis, :, :]) ** 2, axis=-1))
    #     # 1.5 Å is too close
    #     n_clashing += np.count_nonzero(distances < 1.5)
    #     proximals[chain_letter] = np.count_nonzero(distances < 4. )
    #     if i in strep_chain_idxs:
    #         n_strep += np.count_nonzero(distances < 4.)
    # store info
    info['n_clashing'] = n_clashing
    for chain_letter, tally in proximals.items():
        info[f'close_{chain_letter}'] = tally
    info['n_strep'] = n_strep
    if n_clashing > 0:
        print(f'Clashing {target_name}', flush=True)
        info['status'] = 'clashing'
        info['end'] = time.time()
        with Path('log.jsonl').open('a') as f:
            f.write(json.dumps(info) + '\n')
        return info
    elif n_strep < 99:  # 122 for real, but 99 for glycine control
        print(f'Weak {target_name}', flush=True)
        info['status'] = 'weak_strep'
        info['end'] = time.time()
        with Path('log.jsonl').open('a') as f:
            f.write(json.dumps(info) + '\n')
        return info
    monomer: pyrosetta.Pose = thread(**inputs)
    monomer.dump_pdb(f'{raw_out_path}/{target_name}.pdb')
    oligomer: pyrosetta.Pose = monomer.clone()
    add_chain(oligomer, others)  # as this is in place
    oligomer: pyrosetta.Pose = monomer.clone()
    add_chain(oligomer, others)  # as this is in place
    info['original_sequence'] = monomer.sequence()   # repetition for sanity checks



    # constrain
    # KDETET is the start of the oligomerisation domain
    original_start = monomer.sequence().find('KDETET') + 1
    for i in range(original_start, oligomer.total_residue() + 1):
        freeze_atom(pose=oligomer, frozen_index=i, ref_index=1)
    con_scorefxn: pr_scoring.ScoreFunction = pyrosetta.get_fa_scorefxn()
    vanilla_scorefxn: pr_scoring.ScoreFunction = pyrosetta.get_fa_scorefxn()
    con_scorefxn.set_weight(pr_scoring.ScoreType.atom_pair_constraint, 3)
    con_scorefxn.set_weight(pr_scoring.ScoreType.coordinate_constraint, 1)
    vanilla_relax = pyrosetta.rosetta.protocols.relax.FastRelax(con_scorefxn, relax_cycles)
    movemap = pyrosetta.MoveMap()
    # move front part starting from residue 2 and up to KDETET
    v = pru.vector1_bool(oligomer.total_residue())
    for i in range(2, original_start + 1):
        v[i] = 1
    movemap.set_bb(v)
    movemap.set_chi(v)
    movemap.set_jump(False)
    vanilla_relax.set_movemap(movemap)
    vanilla_relax.apply(oligomer)
    info['dG'] = vanilla_scorefxn(oligomer)
    monomer: pyrosetta.Pose = oligomer.split_by_chain(1)
    monomer.dump_pdb(str(relaxed_out_path / f'{target_name}.pdb'))
    print(f'Relaxed {target_name}', flush=True)
    if info['dG'] > -1e3:
        info['status'] = 'score too high'
        info['end'] = time.time()
        with Path('log.jsonl').open('a') as f:
            f.write(json.dumps(info) + '\n')
        return info

    # design
    print(f'Designing {target_name}', flush=True)
    design_different(oligomer, ref, cycles=design_cycles, scorefxn=con_scorefxn)
    monomer = oligomer.split_by_chain(1)
    print(f'designed {target_name}', flush=True)
    monomer.dump_pdb(str(tuned_out_path / f'{target_name}.pdb'))
    info['complex_dG_post'] = vanilla_scorefxn(oligomer)
    info['monomer_dG_post'] = vanilla_scorefxn(monomer)
    info['tweaked_sequence'] = monomer.sequence()
    streptavidins = extract_streptavidins(ref, cardinality=2)
    superpose_pose_by_chain(oligomer, streptavidins, 'B', strict=False)

    # combine with streptavidins only
    minicomplex = monomer.clone()
    minicomplex.remove_constraints()
    add_chain(minicomplex, streptavidins)
    info.update(**score_interface(minicomplex, 'A_BC'))
    info['status'] = 'tuned'
    vanilla_relax.apply(oligomer)
    info['designed_dG'] = vanilla_scorefxn(oligomer)
    pose2pdb = oligomer.pdb_info().pose2pdb
    chainA_sele = prc.select.residue_selector.ChainSelector('A')
    for distance in (1, 2, 3, 4, 5, 6, 8, 10, 12):
        v = prc.select.residue_selector.NeighborhoodResidueSelector(chainA_sele, distance, False).apply(oligomer)
        close_residues = prc.select.get_residues_from_subset(v)
        info[f'N_close_residues_{distance}'] = len(close_residues)
        info[f'close_residues_{distance}'] = [pose2pdb(r) for r in close_residues]
    res_scores = []
    for i in range(1, monomer.total_residue() + 1):
        v = pru.vector1_bool(oligomer.total_residue())
        v[i] = 1
        # score only monomer residues, but on oligomer
        res_scores.append(vanilla_scorefxn.get_sub_score(oligomer, v))
    info['per_res_score'] = res_scores
    info['max_per_res_scores'] = max(res_scores)
    info['end'] = time.time()
    with Path('log.jsonl').open('a') as f:
        f.write(json.dumps(info) + '\n')
    return info
```