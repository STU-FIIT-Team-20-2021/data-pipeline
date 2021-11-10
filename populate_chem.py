from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd
import numpy as np
from math import log
import os


def populate():
    df = pd.read_csv('data/merger_output/merged.csv')

    smiles = df['Smiles'].values.copy()

    for idx, smile in enumerate(smiles):
        smiles[idx] = Chem.MolFromSmiles(smile)

    [AllChem.ComputeGasteigerCharges(y) for y in smiles]
    
    cycles = [len(x.GetRingInfo().BondRings()) for x in smiles]
    count_valence = [sum([x.GetExplicitValence() for x in y.GetAtoms()]) for y in smiles]
    gesteiger_charges = [[float(x.GetProp("_GasteigerCharge")) for x in y.GetAtoms()] for y in smiles]
    positive_gesteiger_sum = [np.sum([y for y in x if y >= 0]) for x in gesteiger_charges]
    negative_gesteiger_sum = [np.sum([y for y in x if y < 0]) for x in gesteiger_charges]
    hydrogen_count = [sum([x.GetAtomicNum() for x in y.GetAtoms() if x.GetAtomicNum() == 1]) for y in smiles]

    df['cycles'] = cycles
    df['atom_valence'] = count_valence
    df['negative_gesteiger'] = negative_gesteiger_sum
    df['positive_gesteiger'] = positive_gesteiger_sum
    df['sum_gesteiger'] = (df['negative_gesteiger'] + df['positive_gesteiger']) * 10 ** 16
    df['hydrogen_count'] = hydrogen_count

    rings = [x.GetRingInfo().BondRings() for x in smiles]
    binarystr = ['0000000000000000'] * len(rings)

    for idx, ring in enumerate(rings):
        for cycle in ring:
            binarystr[idx] = binarystr[idx][:16-(len(cycle)-2)] + '{:x}'.format(int(binarystr[idx][16-(len(cycle)-2)], 16) + 1) + binarystr[idx][16-(len(cycle)-3):]

    converted = []
    for binstr in binarystr:
        converted.append(int(binstr, 16))

    converted = [0 if x == 0 else log(x, 16) for x in converted]
    df['cycle_type_counts'] = converted

    if os.path.isfile('data/final/phototox.csv'):
        df_old = pd.read_csv('data/final/phototox.csv')
        df = pd.concat([df_old, df], axis=0)
        df.to_csv('data/chem_output/chem_populated.csv', index=False)
        df_old.to_csv('data/final/phototox_old.csv')
    else:
        df.to_csv('data/chem_output/chem_populated.csv', index=False)


if __name__ == "__main__":
    populate()
