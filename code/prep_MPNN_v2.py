"""
version 2.

This script generates the files of ProteinMPNN,
based on seq alignment to a reference sequence,
not number of glycines in a row.
"""

import re
import json
import os
import sys
import pickle
import pandas as pd
from pathlib import Path
# this is a refactored version of the original MPNN helpers
from functional_proteinMPNN_helper import (parse_PDBblock, define_fixed_chains,
                                           define_unfixed_positions, define_unfixed_positions)
from typing import Dict, List, Any, Sequence
from Bio.Align import substitution_matrices, PairwiseAligner, Alignment

# -------------------------------------------
work_path = Path(os.environ.get('WORKPATH', 'output'))

# ## chain definitions
# this is the only real JSONL file
# passed to --jsonl_path
chains_definitions_path = work_path / 'chains_definitions.jsonl'
definitions = [] # appended, no need to re-read...

# ## global_fixed_chains
# passed to --chain_id_jsonl
fixed_chains_path = work_path / 'fixed_chains.json'
global_fixed_chains: Dict[str, List[List[str]]] = {}
if fixed_chains_path.exists():
    global_fixed_chains = json.loads(fixed_chains_path.read_text())

# global_fixed_positions
# passed to --fixed_positions_jsonl
fixed_positions_path = work_path / 'fixed_positions.json'
global_fixed_positions = {}
if fixed_positions_path.exists():
    global_fixed_positions = json.loads(fixed_positions_path.read_text())

# ## metadata
# this is for my own reference
meta_path = work_path / 'meta.pkl'
metadata = {}
if meta_path.exists():
    metadata = pickle.loads(meta_path.read_bytes())

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


def get_unfixed_by_alignment(ref_seq: str, pose_seq: str) -> List[int]:
    """
    This does a custom alignment that only allows mutations to Gly
    The values return is a list of 1-based indices of the mutatated positions

    :param ref_seq:
    :param pose_seq:
    :return:
    """
    # ## Define custom substitution_matrix
    # custom weights where only mutations to G are fine
    # this is not going to be really universal?
    # natively there are 3x Gly in a row. This fails if 1 mutation
    # Test:
    # aln = aligner.align(ref_seq, pose_seq)[0]
    # print(aln.score)
    # print(aln.format())
    blosum = substitution_matrices.load("BLOSUM62")
    extend_gap_score = 0.05
    custom = substitution_matrices.Array(alphabet=blosum.alphabet, dims=2)
    for a in custom.alphabet:
        for b in custom.alphabet:
            if a == 'G' and b == 'G':
                custom[a, b] = extend_gap_score * 2
            elif a == b:
                custom[a, b] = 1
            else:
                custom[a, b] = -10
    # ## align
    aligner = PairwiseAligner()
    aligner.substitution_matrix = custom
    aligner.extend_gap_score = extend_gap_score
    aln: Alignment = aligner.align(ref_seq, pose_seq)[0]
    # ## return all different
    return [int(m)+1 for r, m in zip(aln.indices[0], aln.indices[1]) if m != -1 and ref_seq[r] != pose_seq[m]]


# -------------------------------------------

if __name__ == '__main__':
    # input_folder = 'output_redux/filtered'
    input_folder = sys.argv[1]

    n = -1  # <-- positive number: for testing
    paths = [path for path in Path(input_folder).glob('*.pdb')]

    ref_seq = get_chainA_sequence(Path('pentakaihemimer_renumbered.pdb').read_text())
    for path in paths:
        if 'complex' in path.stem:
            continue
        n -= 1
        if n == 0:
            print('break')
            break
        if (work_path / f'seqs/{path.stem}.fa').exists():
            continue
        # chain and seq def
        name = path.stem
        print(f'Prepping {name}', flush=True)
        complex_pdbblock = path.read_text()
        sequence = get_chainA_sequence(complex_pdbblock)
        definition = parse_PDBblock(complex_pdbblock, name)
        definitions.append(definition)
        # fixed chain
        fixed_chains = define_fixed_chains(definition, 'A')
        global_fixed_chains[name] = fixed_chains
        # fixed pos
        # masked = re.sub(rprep_MPNN_v2.py'(G{3,})', lambda match: '_' * len(match.group(1)), sequence)
        # unfixed_list = [[i for i, r in enumerate(masked) if r == '_']]
        unfixed_list = get_unfixed_by_alignment(ref_seq, sequence)
        fixed_positions = define_unfixed_positions(definition, ['A'], [unfixed_list])
        global_fixed_positions[name] = fixed_positions
        metadata[name] = {'name': name, 'sequence': sequence, 'unfixed': unfixed_list}

    # write out definitions, global_fixed_chains, global_fixed_positions
    with open(chains_definitions_path, 'a') as fh:
        for definition in definitions:
            # jsonl
            fh.write(json.dumps(definition) + '\n')

    with open(fixed_chains_path, 'w') as fh:
        json.dump(global_fixed_chains, fh)
        fh.write('\n') # why?

    with open(fixed_positions_path, 'w') as fh:
        json.dump(global_fixed_positions, fh)
        fh.write('\n') # why?

    with open(meta_path, 'wb') as fh:
        pickle.dump(metadata, fh)

