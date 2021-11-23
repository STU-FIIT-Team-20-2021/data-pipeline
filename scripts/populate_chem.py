"""Calculating features using rdkit."""
from math import log
import os

from rdkit import Chem
from rdkit.Chem import AllChem
import rdkit.Chem.rdMolDescriptors as Mol
import rdkit.Chem.Crippen as Crip
import rdkit.Chem.Descriptors as Desc
import numpy as np
import pandas as pd


def populate(file=None) -> pd.DataFrame:
    """
    Takes smiles from the output .csv of merger, and populates the csv with chemical descriptors based on them.

    The default input path is: "data/merger_output/merged.csv"
    the default output path is: "data/chem_output/chem_populated.csv

    The calculated descriptors currently are:
    cycles
    atom_valence
    gesteiger sum and mean of positives, negatives, then total mean and total sum
    hydrogen_count
    csp3
    atomatic_heavy_atoms
    balban_j
    hallkier
    labute_asa
    h_donors and h_acceptors
    carbon_count
    tpsa
    druglikeness
    vsa_logP
    cycle_type_counts

    :param file: Input file for the function for testing, otherwise leave empty for pipeline automation
    :return The modified dataframe, but it's also edited inplace
    """

    df = pd.read_csv(file if file else 'data/merger_output/merged.csv')

    smiles = df['Smiles'].values.copy()

    for idx, smile in enumerate(smiles):
        smiles[idx] = Chem.AddHs(Chem.MolFromSmiles(smile))

    [AllChem.ComputeGasteigerCharges(y) for y in smiles]

    cycles = [len(x.GetRingInfo().BondRings()) for x in smiles]
    count_valence = [sum([x.GetExplicitValence() for x in y.GetAtoms()]) for y in smiles]
    gesteiger_charges = [[float(x.GetProp("_GasteigerCharge")) for x in y.GetAtoms()]
                         for y in smiles]
    positive_gesteiger_sum = [np.sum([y for y in x if y >= 0]) for x in gesteiger_charges]
    negative_gesteiger_sum = [np.sum([y for y in x if y < 0]) for x in gesteiger_charges]
    positive_gesteiger_mean = [np.mean([y for y in x if y >= 0]) for x in gesteiger_charges]
    negative_gesteiger_mean = [np.mean([y for y in x if y < 0]) for x in gesteiger_charges]
    total_gesteiger_sum = [np.sum([y for y in x]) * 10 ** 15 for x in gesteiger_charges]
    total_gesteiger_mean = [np.mean([y for y in x]) * 10 ** 15 for x in gesteiger_charges]
    hydrogen_count = [sum([x.GetAtomicNum() for x in y.GetAtoms() if x.GetAtomicNum() == 1])
                      for y in smiles]
    carbon_count = [sum([x.GetAtomicNum() for x in y.GetAtoms() if x.GetAtomicNum() == 6])
                    for y in smiles]
    csp3 = [Mol.CalcFractionCSP3(y) for y in smiles]
    aha = [len([x for x in y.GetAtoms() if x.GetAtomicNum() != 1 and x.GetIsAromatic()])
           for y in smiles]
    drug = [sum([0 if -0.4 <= Crip.MolLogP(y, True) <= 5.6 else 1,
                0 if 160 <= Desc.ExactMolWt(y) <= 480 else 1,
                0 if 20 <= len(y.GetAtoms()) <= 70 else 1,
                0 if 40 <= Crip.MolMR(y) <= 130 else 1,
                0 if Mol.CalcNumRotatableBonds(y) <= 10 else 1]) for y in smiles]
    kappa = [Mol.CalcKappa1(y) for y in smiles]
    chi = [Mol.CalcChi0n(y) for y in smiles]
    hallkier = [Mol.CalcHallKierAlpha(y) for y in smiles]
    labute_asa = [Mol.CalcLabuteASA(y) for y in smiles]
    phi = [Mol.CalcPhi(y) for y in smiles]
    h_donors = [Mol.CalcNumHBD(y) for y in smiles]
    h_acceptors = [Mol.CalcNumHBA(y) for y in smiles]
    tpsa = [Mol.CalcTPSA(y) for y in smiles]
    vsa = np.array([Mol.SlogP_VSA_(y) for y in smiles])

    df['cycles'] = cycles
    df['atom_valence'] = count_valence
    df['negative_gesteiger_sum'] = negative_gesteiger_sum
    df['positive_gesteiger_sum'] = positive_gesteiger_sum
    df['negative_gesteiger_mean'] = negative_gesteiger_mean
    df['positive_gesteiger_mean'] = positive_gesteiger_mean
    df['sum_gesteiger_fromsum'] = (df['negative_gesteiger_sum']
                                   + df['positive_gesteiger_sum']) * 10 ** 15
    df['sum_gesteiger_frommean'] = (df['negative_gesteiger_mean'] + df['positive_gesteiger_mean'])
    df['sum_gesteiger_totalsum'] = total_gesteiger_sum
    df['sum_gesteiger_totalmean'] = total_gesteiger_mean
    df['hydrogen_count'] = hydrogen_count
    df['csp3'] = csp3
    df['aromatic_heavy_atoms'] = aha
    df['chi0n'] = chi
    df['hallkier'] = hallkier
    df['kappa_1'] = kappa
    df['labute_asa'] = labute_asa
    df['phi'] = phi
    df['h_donors'] = h_donors
    df['h_acceptors'] = h_acceptors
    df['carbon_count'] = carbon_count
    df['tpsa'] = tpsa
    df['druglikeliness'] = drug

    for i, vsa_vals in enumerate(vsa.transpose()):
        df[f'vsa_logP_{i}'] = vsa_vals

    rings = [x.GetRingInfo().BondRings() for x in smiles]
    binarystr = ["0000000000000000"] * len(rings)

    for idx, ring in enumerate(rings):
        for cycle in ring:
            binarystr[idx] = (binarystr[idx][:16-(len(cycle)-2)]
                              + "{:x}".format(int(binarystr[idx][16-(len(cycle)-2)], 16) + 1)
                              + binarystr[idx][16-(len(cycle)-3):])

    converted = []
    for binstr in binarystr:
        converted.append(int(binstr, 16))

    converted = [0 if x == 0 else log(x, 16) for x in converted]
    df['cycle_type_counts'] = converted

    if os.path.isfile("data/final/phototox.csv"):
        df_old = pd.read_csv("data/final/phototox.csv")
        df = pd.concat([df_old, df], axis=0)
        df.to_csv("data/chem_output/chem_populated.csv", index=False)
        df_old.to_csv("data/final/phototox_old.csv")
    else:
        df.to_csv("data/chem_output/chem_populated.csv", index=False)

    return df

if __name__ == "__main__":
    populate("data/merger_output/merged.csv")
