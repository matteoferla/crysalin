import re
import json
import os
import pandas as pd
# this is a refactored version of the original MPNN helpers
from functional_proteinMPNN_helper import parse_PDBblock, define_fixed_chains, define_unfixed_positions, define_unfixed_positions
from typing import Dict, List, Any, Sequence

# -------------------------------------------
data_path = 'scored_complexes3.pkl.gz'
chains_definitions_path = 'output_MPNN/chains_definitions.jsonl'
fixed_chains_path = 'output_MPNN/fixed_chains.json'
fixed_positions_path = 'output_MPNN/fixed_positions.json'

# -------------------------------------------
# read data

# df1 = pd.read_pickle('scored_complexes.pkl.gz').set_index('name')
# df2 = pd.read_pickle('scored_complexes2.pkl.gz') #.set_index('name')
#df = pd.concat([df2, df1])
df = pd.read_pickle(data_path)
df = df.reset_index().drop_duplicates('name', keep='first').set_index('name')
print(f'{len(df)} entries')

# -------------------------------------------

no_fix = []
definitions = []
global_fixed_chains: Dict[str, List[List[str]]] = {}
global_fixed_positions = {}

for name, row in df.iterrows():
    if str(row.minicomplex) == 'nan':
        continue
    elif row.dG_bind > 5e3:
        continue
    elif row.N_close_residues_4 <= 2:
        continue
    # chain and seq def
    definition = parse_PDBblock(row.minicomplex, row.name)
    definitions.append(definition)
    # fixed chain
    fixed_chains = define_fixed_chains(definition, 'A')
    global_fixed_chains[name] = fixed_chains
    # fixed pos
    masked = re.sub(r'(G{3,})', lambda match: '_' * len(match.group(1)), row.sequence)
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

