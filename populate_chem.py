from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd
import numpy as np
import sys
import os

if sys.stdin.isatty():
    df = pd.read_csv(sys.stdin)
else:
    df = pd.read_csv('./merger/merge.csv')

smiles = df['Smiles'].values.copy()

for idx, smile in enumerate(smiles):
    smiles[idx] = Chem.MolFromSmiles(smile)
    

[AllChem.ComputeGasteigerCharges(y) for y in smiles]
    
cycles = [len(x.GetRingInfo().BondRings()) for x in smiles]
count_valence = [sum([x.GetExplicitValence() for x in y.GetAtoms()]) for y in smiles]
gesteiger_charges = [[float(x.GetProp("_GasteigerCharge")) for x in y.GetAtoms()] for y in smiles]
positive_gesteiger_sum = [np.sum([y for y in x if y >= 0]) for x in gesteiger_charges]
negative_gesteiger_sum = [np.sum([y for y in x if y < 0]) for x in gesteiger_charges]
hydrogen_count = [sum([x.GetExplicitValence() for x in y.GetAtoms()]) for y in smiles]

df['cycles'] = cycles
df['atom_valence'] = count_valence
df['negative_gesteiger'] = negative_gesteiger_sum
df['positive_gesteiger'] = positive_gesteiger_sum

rings = [x.GetRingInfo().BondRings() for x in smiles]
binarystr = ['0000000000000000'] * len(rings)

for idx, ring in enumerate(rings):
    for cycle in ring:
        binarystr[idx] = binarystr[idx][:16-(len(cycle)-2)] + '{:x}'.format(int(binarystr[idx][16-(len(cycle)-2)], 16) + 1) + binarystr[idx][16-(len(cycle)-3):]

converted = []
for binstr in binarystr:
    converted.append(int(binstr, 16))

df['cycle_type_counts'] = converted

df.to_csv('./chem/new_attrib.csv', index=False)

if not os.name == 'nt':
    df.to_csv(sys.stdout, index=False)
