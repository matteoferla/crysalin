import re
import json
import os
import pandas as pd
from pathlib import Path
# this is a refactored version of the original MPNN helpers
from functional_proteinMPNN_helper import parse_PDBblock, define_fixed_chains, define_unfixed_positions, define_unfixed_positions
from typing import Dict, List, Any, Sequence

# -------------------------------------------
data_path = 'scored_complexes3.pkl.gz'
input_folder = 'output_redux/filtered'
chains_definitions_path = 'output_MPNN/chains_definitions.jsonl'
fixed_chains_path = 'output_MPNN/fixed_chains.json'
fixed_positions_path = 'output_MPNN/fixed_positions.json'
woA_block = Path('woA_pentakaihemimer.pdb').read_text()

# -------------------------------------------

no_fix = []
definitions = []
global_fixed_chains: Dict[str, List[List[str]]] = {}
global_fixed_positions = {}

def get_ATOM_only(pdbblock: str) -> str:
    return '\n'.join([line for line in pdbblock.splitlines() if line.startswith('ATOM')])

three_to_one = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',
    'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
    'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
    'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}

def get_sequence(pdbblock: str) -> str:
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


for path in Path(input_folder).glob('*.pdb'):
    if 'complex' in path.stem:
        continue
    # chain and seq def
    name = path.stem
    momomer_pdbblock = get_ATOM_only(path.read_text())
    complex_pdbblock = momomer_pdbblock + '\n' + woA_block
    sequence = get_sequence(complex_pdbblock)
    definition = parse_PDBblock(complex_pdbblock, name)
    definitions.append(definition)
    # fixed chain
    fixed_chains = define_fixed_chains(definition, 'A')
    global_fixed_chains[name] = fixed_chains
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

