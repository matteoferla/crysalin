import re
import json
import os
import pandas as pd
from pathlib import Path
# this is a refactored version of the original MPNN helpers
from functional_proteinMPNN_helper import (parse_PDBblock, define_fixed_chains,
                                           define_unfixed_positions, define_unfixed_positions)
from typing import Dict, List, Any, Sequence

# -------------------------------------------
work_path = Path(os.environ.get('WORKPATH', 'output'))
#input_folder = 'output_redux/filtered'
input_folder = f'{work_path}/superposed_skeletons'
chains_definitions_path = f'{work_path}/chains_definitions.jsonl'
fixed_chains_path = f'{work_path}/fixed_chains.json'
fixed_positions_path = f'{work_path}/fixed_positions.json'

# -------------------------------------------

no_fix = []
definitions = []
global_fixed_chains: Dict[str, List[List[str]]] = {}
global_fixed_positions = {}

def get_ATOM_only(pdbblock: str) -> str:
    """
    This gets all ATOM, regardless of name and chain
    """
    return '\n'.join([line for line in pdbblock.splitlines() if line.startswith('ATOM')])

three_to_one = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',
    'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
    'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
    'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}

def get_chainA_sequence(pdbblock: str) -> str:
    sequence = ''
    residues_seen = set()
    for line in pdbblock.splitlines():
        if line.startswith("ATOM") and " CA " in line and " A " in line:
            res_info = line[17:26]  # Residue name and number for uniqueness
            if res_info not in residues_seen:
                residues_seen.add(res_info)
                res_name = line[17:20].strip()
                sequence += three_to_one.get(res_name, '?')
    return sequence

if __name__ == '__main__':
    n = -1 # <-- positive number: for testing
    paths = [path for path in Path(input_folder).glob('*.pdb')]
    woA_block = Path('woA_pentakaihemimer.relax.pdb').read_text()
    for path in paths:
        if 'complex' in path.stem:
            continue
        # break after n if enabled
        n -= 1
        if n == 0:
            print('break')
            break
        # skip if seqs exists
        if (work_path / f'seqs/{path.stem}.fa').exists():
            continue
        # chain and seq def
        name = path.stem
        print(f'Prepping {name}', flush=True)
        complex_pdbblock = get_ATOM_only(path.read_text())
        definition = parse_PDBblock(complex_pdbblock, name)
        definitions.append(definition)
        fixed_chains = define_fixed_chains(definition, 'A')
        global_fixed_chains[name] = fixed_chains

        momomer_pdbblock = get_ATOM_only(path.read_text())
        complex_pdbblock = momomer_pdbblock + '\n' + woA_block
        sequence = get_chainA_sequence(complex_pdbblock)

        # fixed chain


        # fixed pos
        masked = re.sub(r'(G{3,})', lambda match: '_' * len(match.group(1)), sequence)
        fixed_list = [[i for i, r in enumerate(masked) if r == '_']]
        fixed_positions = define_unfixed_positions(definition, ['A'], fixed_list)
        global_fixed_positions[name] = fixed_positions

    # write out definitions, global_fixed_chains, global_fixed_positions
    with open(chains_definitions_path, 'w') as fh:
        for definition in definitions:
            fh.write(json.dumps(definition) + '\n')

    with open(fixed_chains_path, 'w') as fh:
        json.dump(global_fixed_chains, fh)
        fh.write('\n')

    with open(fixed_positions_path, 'w') as fh:
        json.dump(global_fixed_positions, fh)
        fh.write('\n')

