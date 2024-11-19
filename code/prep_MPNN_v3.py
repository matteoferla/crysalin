"""
version 3.

This script generates the files of ProteinMPNN,
based on the trb file, which I did not realise had the sequence.
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
class DefinitionStore:
    def __init__(self, work_path = None):
        # work path
        if work_path is None:
            work_path = Path(os.environ.get('WORKPATH', 'output'))
        self.work_path = work_path

        # ## chain definitions
        # this is the only real JSONL file
        # passed to --jsonl_path
        self.chains_definitions_path = work_path / 'chains_definitions.jsonl'
        self.definitions = []  # appended, no need to re-read...

        # ## global_fixed_chains
        # passed to --chain_id_jsonl
        self.fixed_chains_path = work_path / 'fixed_chains.json'
        self.global_fixed_chains: Dict[str, List[List[str]]] = {}

        # global_fixed_positions
        # passed to --fixed_positions_jsonl
        self.fixed_positions_path = work_path / 'fixed_positions.json'
        self.global_fixed_positions = {}

    def read(self, including_definitions=False):
        if including_definitions and self.chains_definitions_path.exists():
            self.definitions = [json.loads(line) for line in self.chains_definitions_path.read_text().splitlines()]
        if self.fixed_chains_path.exists():
            self.global_fixed_chains = json.loads(self.fixed_chains_path.read_text())
        if self.fixed_positions_path.exists():
            self.global_fixed_positions = json.loads(self.fixed_positions_path.read_text())
        return self

    def write(self):

        with open(self.chains_definitions_path, 'a') as fh:
            for definition in self.definitions:
                # jsonl
                fh.write(json.dumps(definition) + '\n')

        with open(self.fixed_chains_path, 'w') as fh:
            json.dump(self.global_fixed_chains, fh)
            fh.write('\n')  # why?

        with open(self.fixed_positions_path, 'w') as fh:
            json.dump(self.global_fixed_positions, fh)
            fh.write('\n')  # why?
        return self



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


def get_munged_fixed_def_by_trb(path: Path) -> Dict[str, List[int]]:
    """
    reads the trb file and returns the list of 1-based indices of the mutatated positions

    :param ref_seq:
    :param pose_seq:
    :return:
    """
    trb_path = path.parent / (path.stem + '.trb')
    assert trb_path.exists(), f'{trb_path} does not exist'
    trb: Dict[str, Any] = pickle.load(trb_path.open('rb'))
    definitions = {}
    for chain, resi in trb['con_hal_pdb_idx']: # ('A', 1)
        if chain not in definitions:
            definitions[chain] = []
        definitions[chain].append(int(resi))
    # ## return all different
    return definitions

def check_done(root_folder='output', ongoing = ()):
    """
    This is not called by the script.
    It is a snippet to check what worked.
    Prints the stats of each folder in the output folder.

    :param root_folder:
    :return:
    """
    import json

    repeaters = []

    for path in Path(root_folder).glob('*'):
        if not path.is_dir() or path.name == '.ipynb_checkpoints':
            continue
        is_ongoing = path.stem in ongoing
        if (path / 'fixed_chains.json').exists():
            n_expected = len(json.load((path / 'fixed_chains.json').open('r')))
        else:
            n_expected = -1
        pdb_names = [p.stem for p in path.glob('*.pdb')]
        if not (path / 'seqs').exists():
            seqs_paths = []
            diff = pdb_names
        else:
            pdb_names = [p.stem for p in path.glob('*.pdb')]
            seqs_paths = [p.stem for p in (path / 'seqs').glob('*.fa')]
            diff = set(pdb_names) - set(seqs_paths)
        if len(diff) > 0 and not is_ongoing:
            repeaters.append(path)
        print(f'{path} - missing: {len(diff)}, pdbs: {len(pdb_names)}, seqs: {len(seqs_paths)}, json: {n_expected}. Is ongoing {is_ongoing}')
    return repeaters

# -------------------------------------------

if __name__ == '__main__':
    # input_folder = 'output_redux/filtered'
    input_folder = sys.argv[1]

    work_path = Path(os.environ.get('WORKPATH', 'output'))
    defstore = DefinitionStore(work_path)
    defstore.read(including_definitions=False)

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
        definition = parse_PDBblock(complex_pdbblock, name)
        defstore.definitions.append(definition)
        # fixed chain
        fixed_chains = define_fixed_chains(definition, designed_chain_list='A')
        defstore.global_fixed_chains[name] = fixed_chains
        # fixed pos
        defstore.global_fixed_positions[name]: Dict[str, List[int]] = get_munged_fixed_def_by_trb(path)
        # write out definitions, global_fixed_chains, global_fixed_positions etc.
        defstore.write()
    print(f'Processed {len(defstore.definitions)} models')



