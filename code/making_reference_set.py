import pickle
from pathlib import Path
import yaml

with open('experiments-new.yaml', 'r') as file:
    experiment_definitions = yaml.safe_load(file)
experiment_definitions['original'] = {'codename': 'original', 'contig': '[A5-331/0 B1-651/0]', }


def make_fake_trb(contig):
    trb = {'con_ref_idx0': [], 'complex_con_ref_idx0': [], 'inpaint_seq': [], }
    pdb_idx = None
    idx0 = 0
    previus_aft = 9999
    for subcontig in contig[1:-1].split():
        for span in subcontig[:-2].split('/'):
            if span[0].isdecimal():
                continue
            fore, aft = span[1:].split('-')
            if span[0] == 'A':
                for pdb_idx in range(previus_aft + 1, int(fore)):
                    trb['inpaint_seq'].append(False)
                    idx0 += 1
                previus_aft = int(aft)
            # extend inpainted chain A+B
            for pdb_idx in range(int(fore), int(aft) + 1):
                trb['complex_con_ref_idx0'].append(idx0)
                if span[0] == 'A':
                    trb['con_ref_idx0'].append(idx0)
                trb['inpaint_seq'].append(True)
                idx0 += 1
    trb['complex_con_hal_idx0'] = trb['complex_con_ref_idx0']
    trb['con_hal_idx0'] = trb['con_ref_idx0']
    return trb


from Bio.Data.IUPACData import protein_letters_3to1


def extract_sequence_from_block(pdb_block):
    sequence = []
    for line in pdb_block.splitlines():
        if line.startswith('ATOM') and line[21] == 'A':
            residue_name = line[17:20].strip()
            residue_id = line[22:26].strip()
            atom_name = line[12:16].strip()
            if atom_name == 'CA':
                sequence.append(protein_letters_3to1.get(residue_name.strip().capitalize(), 'X'))
    return ''.join(sequence)


tunings = []
template_block = open('pentakaihemimer_renumbered2.pdb').read()
bb_names = ['N', 'CA', 'C', 'O']
shaven_block = '\n'.join(
    [line for line in template_block.split('\n') if 'ATOM' in line and line[13:16].strip() in bb_names])

for group, defs in experiment_definitions.items():
    contig = defs['contig']
    name = defs['codename']
    trb = make_fake_trb(contig)
    with open(f'output/ref/ref_{name}.trb', 'wb') as fh:
        pickle.dump(trb, fh)
    lines = shaven_block.split('\n')
    pacified_block = ''
    for idx0, is_frozen in enumerate(trb['inpaint_seq']):
        for o in range(4):
            line = lines[idx0 * 4 + o]
            if not is_frozen:
                line = line[:17] + 'GLY' + line[20:]
            pacified_block += line + '\n'
    with open(f'output/ref/ref_{name}.pdb', 'w') as fh:
        fh.write(pacified_block)

    hallucinated_seq = extract_sequence_from_block(pacified_block)
    template_seq = extract_sequence_from_block(template_block)
    with open(f'output/ref/seqs/ref_{name}.fa', 'w') as fh:
        fh.write(
            f">ref_{name}, score=0.0, global_score=0.0, fixed_chains=['B'], designed_chains=['A'], model_name=v_48_020, git_hash=unknown, seed=888\n")
        fh.write(f"{hallucinated_seq}\n")

        fh.write(f">T=0.1, sample=1, score=0.0, global_score=0.0, seq_recovery=0.0\n")
        fh.write(f"{template_seq}\n")

    tunings.append(dict(target_folder=Path('/opt/xchem-fragalysis-2/mferla/crysalin-redux/output/ref'),
                        parent_filename='/opt/xchem-fragalysis-2/mferla/crysalin-redux/pentakaihemimer_renumbered2.pdb',
                        target_name=f'ref_{name}Ø',
                        target_sequence=hallucinated_seq,
                        hallucination_name=f'ref_{name}',
                        metadata={},
                        ))
    tunings.append(dict(target_folder=Path('/opt/xchem-fragalysis-2/mferla/crysalin-redux/output/ref'),
                        parent_filename='/opt/xchem-fragalysis-2/mferla/crysalin-redux/pentakaihemimer_renumbered2.pdb',
                        target_name=f'ref_{name}A',
                        target_sequence=template_seq,
                        hallucination_name=f'ref_{name}',
                        metadata={},
                        ))

if False == True:
    from pathlib import Path
    import sys

    sys.path.append('/opt/xchem-fragalysis-2/mferla/crysalin-redux/repo/code')
    from repo.code.thread_tune_v2 import run_process as tune_single
    from repo.code.thread_tune_v2 import main as tune_enamass
    from repo.code.thread_tune_v2 import init_pyrosetta, SETTINGS

    SETTINGS['exception_to_catch'] = Exception
    SETTINGS['clash_dist_cutoff'] = 0.1
    SETTINGS['bond_dist_cutoff'] = 2.0

    init_pyrosetta()

    for tuning in tunings:
        tune_single(**tuning)
