from libxml2 import treeError
from unittest.mock import inplace

import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
import matplotlib.pyplot as plt
import os.path

import logging
import networkx as nx

logging.basicConfig()
logging.getLogger().setLevel(logging.INFO)


# poands discribe
# by model analysis vuci phototox


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


#  *************************** NEW ************************************
def load_data(csv_file_name):
    return pd.read_csv(csv_file_name)


def set_proper_data_type(df):
    replace_dict = {'Low': False, 'High': True, 'Yes': True, 'No': False}

    for repl, val in replace_dict.items():
        df.replace(repl, val, inplace=True)

    dtypes_columns = {
        bool: ['Phototoxic', 'Compound Is Canonicalized', 'BBB permeant', 'Pgp substrate', 'CYP1A2 inhibitor',
               'CYP2C9 inhibitor', 'CYP2D6 inhibitor', 'CYP3A4 inhibitor', 'GI absorption'],
        str: ['Name', 'Smiles', 'Formula', ],
        float: ['Molecular Weight', 'Exact Mass', 'Monoisotopic Mass', 'Topological Polar Surface Area',
                'Complexity', 'XLogP3', 'XLogP3-AA'],
        int: ['Hydrogen Bond Donor Count', 'Hydrogen Bond Acceptor Count', 'Rotatable Bond Count',
              'Heavy Atom Count', 'Formal Charge', 'Isotope Atom Count',
              'Defined Bond Stereocenter Count', 'Undefined Atom Stereocenter Count',
              'Covalently-Bonded Unit Count']
    }

    all_columns = df.columns
    dtypes_columns = {col: dt for dt, columns in dtypes_columns.items() for col in columns if col in all_columns}
    df = df.astype(dtypes_columns)

    return df


def unique_count(column):
    return column.name, len(column.unique())


def remove_useless_columns(df):
    # its only cathegory based on previous values
    df.drop(columns={'ESOL Class', 'Ali Class', 'Silicos-IT class'}, inplace=True)
    # only different unit
    df.drop(columns={'ESOL Solubility (mol/l)', 'Ali Solubility (mol/l)', 'Silicos-IT Solubility (mol/l)'},
            inplace=True)

    # remove columns with one value only
    column_unique_count = df.apply(unique_count, axis=0).transpose()
    columns_to_remove = column_unique_count[column_unique_count[1] <= 1][0].to_list()

    if len(columns_to_remove):
        logging.info(f'Removing columns {columns_to_remove} with only one unique value')
        df.drop(columns=columns_to_remove, inplace=True)

    # remove duplicit columns in pubchem - they contain values based on bad smiles code
    punches_duplicity_columns = ['Heavy Atom Count', 'Hydrogen Bond Donor Count', 'Hydrogen Bond Acceptor Count',
                                 'Rotatable Bond Count', 'Topological Polar Surface Area', 'Molecular Weight']
    # todo check na tretnuti
    logp_duplicity_columns = ['XLogP3-AA', 'XLogP3']
    df.drop(columns=punches_duplicity_columns + logp_duplicity_columns, inplace=True)

    # remove duplicit columns caused by merge
    df.drop(columns=['Smiles.1', 'Name.1'], inplace=True)
    return df


def get_correlated_descriptors_to_delete(graphs) -> list:
    return [node for sg in graphs for i, node in enumerate(sg.nodes) if i != 0]


def get_correlated_descriptors(df, threshold, fig_location):
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

    fig.suptitle(f'Strong correlations values ({threshold})', fontsize=16)

    plt.tight_layout()
    plt.show()
    plt.close()

    fig.savefig(fig_location)

    return correlations


def create_correlation_graphs(correlations, fig_location, trashold):
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
        nx.draw(sg, with_labels=True, font_weight='bold', font_size=8, pos=nx.circular_layout(sg), ax=axes[i])
        logging.info(f'Subgraph {i}: {nodes}')

    fig.suptitle(f'Strong correlations graph {trashold}', fontsize=16)

    plt.tight_layout()
    plt.show()
    plt.close()

    fig.savefig(fig_location)

    return subgraphs


def check_correlated_column(df, trashhold=0.9, remove=False, graph_location='./plot/correlations_grapgs.jpeg',
                            heatmap_location='./plot/correlations_heatmap.jpeg') -> pd.DataFrame:
    correlations = get_correlated_descriptors(df, trashhold, heatmap_location)
    graphs = create_correlation_graphs(correlations, graph_location, trashhold)
    to_delete = get_correlated_descriptors_to_delete(graphs)
    if remove:
        df.drop(to_delete, axis=1, inplace=True)
    logging.info(f'Removing correlated columns: {to_delete}')
    return df


def remove_duplicits(df: pd.DataFrame, subset, keep='first'):
    # todo check phototox
    duplicits = df[df.duplicated(subset=subset, keep=keep)]

    logging.warning(f"Found duplicates:{duplicits[subset].to_list()} deleted.")
    df.drop_duplicates(subset=subset, keep=keep, inplace=True)

    return df


def check_outliers(df, remove=False):
    numerics = ['int16', 'int32', 'int64', 'float16', 'float32', 'float64']
    df_numeric = df.select_dtypes(include=numerics)

    outliers = df[~(np.abs(stats.zscore(df_numeric)) < 4.2).all(axis=1)]['Name'].unique()

    logging.info(f'Outliners: {outliers}')

    if remove:
        df = df[~df['Name'].isin(outliers)]
        logging.info(f'Outliners {outliers} deleted.')

    return df


def main():
    input = 'data/chem_output/chem_populated.csv'
    output = 'data/final/phototox.csv'

    df = load_data(input)
    df = remove_useless_columns(df)
    df = set_proper_data_type(df)

    df = check_correlated_column(df)
    df = remove_duplicits(df, subset="Smiles", keep='first')
    df = check_outliers(df)

    df.to_csv(output, index=False)


if __name__ == "__main__":
    main()
