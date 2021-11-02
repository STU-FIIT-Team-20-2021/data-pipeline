from unittest.mock import inplace

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os.path

import logging
import networkx as nx

logging.basicConfig()
logging.getLogger().setLevel(logging.INFO)


# poands discribe
# by model analysis vuci phototox

# *********  Other ********************************************************

def drop_duplicits(df: pd.DataFrame, subset, keep='first'):
    logging.warning(f"Found duplicates:{df[df.duplicated(subset=subset, keep=keep)][subset].to_list()} deleted.")
    df.drop_duplicates(subset=subset, keep=keep, inplace=True)


# ********* PUBCHEM ********************************************************

def load_multiple_csv(csv_names: list):
    all_df = [pd.read_csv(name) for name in csv_names]
    return pd.concat(all_df)


def clean_pubchem(pubchem):
    pubchem = pubchem.drop(columns='Unnamed: 0')
    drop_duplicits(pubchem, subset="Smiles", keep='first')

    pubchem.loc[:, 'Topological Polar Surface Area'] = pubchem.loc[:, 'Topological Polar Surface Area'].astype(
        str).str.replace(' Å²', '')
    pubchem.replace({'Compound Is Canonicalized': {'Yes': True, 'No': False}})

    dtypes_columns = {bool: ['Phototoxic', 'Compound Is Canonicalized'], str: ['Name', ],
                      float: ['Molecular Weight', 'Exact Mass', 'Monoisotopic Mass', 'Topological Polar Surface Area',
                              'Complexity', 'XLogP3', 'XLogP3-AA'],
                      int: ['Hydrogen Bond Donor Count', 'Hydrogen Bond Acceptor Count', 'Rotatable Bond Count',
                            'Heavy Atom Count', 'Formal Charge', 'Isotope Atom Count',
                            'Defined Atom Stereocenter Count',
                            'Defined Bond Stereocenter Count', 'Undefined Atom Stereocenter Count',
                            'Covalently-Bonded Unit Count']}
    dtypes_columns = {col: dt for dt, columns in dtypes_columns.items() for col in columns}
    pubchem = pubchem.astype(dtypes_columns)

    return pubchem


def get_pubchem(input_csv_files, output_csv_file='./data/pubchem.csv', replace=True):
    if os.path.isfile(output_csv_file) and not replace:
        return pd.read_csv(output_csv_file)

    pubchem = load_multiple_csv(input_csv_files)
    pubchem = clean_pubchem(pubchem)

    pubchem.to_csv(output_csv_file)
    return pubchem


# ********* SWISS ********************************************************

def clean_swiss(swiss):
    drop_duplicits(swiss, subset="Smiles", keep='first')

    replace_dict = {'Low': False, 'High': True, 'Yes': True, 'No': False}

    for repl, val in replace_dict.items():
        swiss.replace(repl, val, inplace=True)

    # its only cathegory based on previous values
    swiss.drop(columns={'ESOL Class', 'Ali Class', 'Silicos-IT class'}, inplace=True)
    # only different unit
    swiss.drop(columns={'ESOL Solubility (mol/l)', 'Ali Solubility (mol/l)', 'Silicos-IT Solubility (mol/l)'},
               inplace=True)
    # todo check their logs on correlation

    dtypes_columns = {str: ['Smiles', 'Formula', 'Name', 'Formula', ],
                      bool: ['BBB permeant', 'Pgp substrate', 'CYP1A2 inhibitor', 'CYP2C19 inhibitor',
                             'CYP2C9 inhibitor', 'CYP2D6 inhibitor', 'CYP3A4 inhibitor', 'GI absorption']}

    dtypes_columns = {col: dt for dt, columns in dtypes_columns.items() for col in columns}
    swiss = swiss.astype(dtypes_columns)

    return swiss


def get_swiss(input_csv_files, output_csv_file='./data/swiss.csv', replace=True):
    if os.path.isfile(output_csv_file) and not replace:
        return pd.read_csv(output_csv_file)

    swiss = load_multiple_csv(input_csv_files)
    swiss = clean_swiss(swiss)

    swiss.to_csv(output_csv_file)
    return swiss


def log_unique_ds_count(left, right, on='Name'):
    left_count = len(left[left[on].isin(right[on])])
    right_count = len(right[right[on].isin(left[on])])
    total = left_count + right_count
    if total > 0:
        logging.warning(f"Number of molecules not included in both datasets: {left_count + right_count}")


def check_column_consistency(df, mapper: dict, percent_error=None):
    for col1, col2 in mapper.items():
        if percent_error is None:
            is_equal = df[col1] != df[col2]
            if not all(is_equal):
                logging.error(
                    f'Pubchem and swiss columns are not equal. '
                    f'Columns "{col1}" and "{col2}" contain diffrent values by {df[~is_equal][col1].count()} rows')
        else:
            percent_difference = df[col2] / (df[col1] / 100)
            difference = percent_difference < percent_error
            if not all(difference):
                logging.error(
                    f'Pubchem and swiss columns are not equal. '
                    f'Columns "{col1}" and "{col2}" differ more than {percent_error}% by '
                    f'{df[~difference][col1].count()} rows')


def merge_equal_columns(df):
    int_pubchem_swiss_mapping = {'Heavy Atom Count': '#Heavy atoms', 'Hydrogen Bond Donor Count': '#H-bond donors',
                                 'Hydrogen Bond Acceptor Count': '#H-bond acceptors',
                                 'Rotatable Bond Count': '#Rotatable bonds'}
    float_pubchem_swiss_mapping = {'Molecular Weight': 'MW'}

    check_column_consistency(df, int_pubchem_swiss_mapping)
    check_column_consistency(df, float_pubchem_swiss_mapping, 5)

    df.drop(columns=[*int_pubchem_swiss_mapping.keys()], inplace=True)
    df.drop(columns=[*float_pubchem_swiss_mapping.keys()], inplace=True)

    return df


def merge_pubchem_swiss(pubchem, swiss):
    log_unique_ds_count(pubchem, swiss)
    df = pd.merge(pubchem, swiss, how='outer', on='Name')
    df = merge_equal_columns(df)
    return df


def unique_count(column):
    return column.name, len(column.unique())


def remove_useless_columns(df):
    column_unique_count = df.apply(unique_count, axis=0).transpose()
    columns_to_remove = column_unique_count[column_unique_count[1] <= 1][0].to_list()

    if len(columns_to_remove):
        logging.info(f'Removing columns {columns_to_remove} with only one unique value')
        df.drop(columns=columns_to_remove, inplace=True)

    return df


def get_correlated_descriptors(df, threshold):
    correlations = df.corr().abs()

    # half of table only
    mask = np.triu(np.ones_like(correlations, dtype=bool))

    # high correlations only, other values nan
    correlations.mask(mask, inplace=True)
    correlations = correlations.applymap(lambda x: x if x > threshold else np.nan)
    correlations = correlations.loc[~(correlations.isna()).all(axis=1)]
    correlations = correlations.loc[:, ~(correlations.isna()).all(axis=0)]

    fig, ax = plt.subplots(figsize=(8.27, 8.27))
    sns.heatmap(correlations, cmap='gist_heat', square=True, ax=ax)

    fig.suptitle('Strong correlations values', fontsize=16)

    plt.tight_layout()
    plt.show()
    plt.close()

    fig.savefig('./correlations_heatmap.jpeg')

    return correlations


def create_correlation_graphs(correlations):
    MG = nx.MultiGraph()
    columns = correlations.columns

    # get edges between nodes
    for index, row in correlations.iterrows():
        new_edges = [(index, col, row[col]) for col in columns if not np.isnan(row[col])]
        MG.add_weighted_edges_from(new_edges)

    subgraphs = [MG.subgraph(c) for c in nx.connected_components(MG)]

    fig, axes = plt.subplots(len(subgraphs), figsize=(8.27, 11.69))
    axes = axes.flatten()

    for i, sg in enumerate(subgraphs):
        nodes = sg.nodes
        nx.draw(sg, with_labels=True, font_weight='bold' , font_size=8 , pos=nx.circular_layout(sg), ax=axes[i])
        logging.info(f'Subgraph {i}: {nodes}')

    fig.suptitle('Strong correlations graph', fontsize=16)

    plt.tight_layout()
    plt.show()
    plt.close()

    fig.savefig('./correlations_grapgs.jpeg')

    return subgraphs


def get_correlated_descriptors_to_delete(graphs) -> list:
    return [node for sg in graphs for i, node in enumerate(sg.nodes) if i != 0]


def remove_correlated_column(df, trashhold=0.9) -> pd.DataFrame:
    correlations = get_correlated_descriptors(df, trashhold)
    graphs = create_correlation_graphs(correlations)
    to_delete = get_correlated_descriptors_to_delete(graphs)
    df.drop(to_delete, axis=1, inplace=True)
    logging.info(f'Removing correlated columns: {to_delete}')
    return df


def main():
    pubchem = get_pubchem(["./data/pubchem_raw.csv", ], replace=True)
    swiss = get_swiss(["./data/swissadme_raw.csv", ], replace=True)

    df = merge_pubchem_swiss(pubchem, swiss)
    df = remove_useless_columns(df)
    df = remove_correlated_column(df)

    df.to_csv('./data/cleaned.csv', index=False)


if __name__ == "__main__":
    main()