from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt

df = pd.read_csv('./data/new_merged.csv')
smiles = df['Canonical SMILES'].values

for idx, smile in enumerate(smiles):
    smiles[idx] = Chem.MolFromSmiles(smile)
    

[AllChem.ComputeGasteigerCharges(y) for y in smiles];
    
cycles = [len(x.GetRingInfo().BondRings()) for x in smiles]
count_valence = [sum([x.GetExplicitValence() for x in y.GetAtoms()]) for y in smiles]
mean_gesteiger_charge = [np.mean([float(x.GetProp("_GasteigerCharge")) for x in y.GetAtoms()]) for y in smiles]

df['cycles'] = cycles
df['atom_valence'] = count_valence
df['gesteiger'] = mean_gesteiger_charge

rings = [x.GetRingInfo().BondRings() for x in smiles]
ring_counts = [0]*len(rings)
min_cycle = 3

for idx, ring in enumerate(rings):
    ring_counts[idx] = [0] * 12
    for cycle in ring:
        ring_counts[idx][len(cycle) - min_cycle] += 1
        
df = pd.concat([df, pd.DataFrame(ring_counts)], axis=1)

df.to_csv('./data/new_attrib.csv', index=False)