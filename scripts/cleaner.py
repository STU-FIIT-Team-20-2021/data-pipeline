import os
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import configparser
import os.path
import logging
import networkx as nx

logging.basicConfig()
logging.getLogger().setLevel(logging.INFO)


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

    df.rename(columns={'Ali Solubility (mg/ml)': 'Ali Solubility', 'ESOL Solubility (mg/ml)': 'ESOL Solubility'})

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
                                 'Rotatable Bond Count', 'Topological Polar Surface Area', 'Molecular Weight',
                                 'XLogP3-AA', 'XLogP3']
    df.drop(columns=punches_duplicity_columns, inplace=True)

    # remove duplicit columns from swiss and computed columns
    swiss_duplicity_columns = ['TPSA', '#H-bond acceptors', '#H-bond donors', '#Aromatic heavy atoms', 'Fraction Csp3',
                               'Ali Solubility (mg/ml)', 'ESOL Solubility (mg/ml)']
    df.drop(columns=swiss_duplicity_columns, inplace=True)

    # remove duplicit columns caused by merge
    df.drop(columns=['Smiles.1', 'Name.1'], inplace=True)
    return df


def get_duplicit_correlated_descriptors(graphs, descriptors) -> list:
    duplicit_correlated = []
    for subgraph in graphs:
        nodes = [node for node in subgraph.nodes]
        nodes_without_preserve = [node for node in nodes if node not in descriptors]

        if len(nodes) == len(nodes_without_preserve):
            duplicit_correlated = duplicit_correlated + nodes[:-1]
        else:
            duplicit_correlated = duplicit_correlated + nodes_without_preserve

    return duplicit_correlated


def get_correlated_descriptors(df, threshold, fig_location):
    correlations = df.corr().abs()

    # half of table only
    mask = np.triu(np.ones_like(correlations, dtype=bool))

    # high correlations only, other values nan
    correlations.mask(mask, inplace=True)
    correlations = correlations.applymap(lambda x: x if x > threshold else np.nan)
    correlations = correlations.loc[~(correlations.isna()).all(axis=1)]
    correlations = correlations.loc[:, ~(correlations.isna()).all(axis=0)]

    # plot
    fig, ax = plt.subplots(1, 1, figsize=(8.27, 8.27))

    img = ax.imshow(correlations, cmap="gist_heat", extent=[-1, 1, -1, 1])
    fig.colorbar(img)

    tick_len = (1 / (len(correlations.columns)))
    ticks = [(i * tick_len + tick_len / 2) * 2 - 1 for i in range(len(correlations.columns))]

    ax.set_xticks(ticks)
    ax.set_yticks(ticks)

    ax.set_xticklabels(correlations.columns.to_list())
    ax.set_yticklabels(correlations.index.to_list()[::-1])
    ax.tick_params(axis='x', labelrotation=90)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)

    fig.suptitle(f'Strong correlations values ({threshold})', fontsize=16)

    plt.tight_layout()
    plt.show()
    plt.close()

    fig.savefig(fig_location)

    return correlations


def create_correlation_graphs(correlations, fig_location, threshold):
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

    fig.suptitle(f'Strong correlations graph {threshold}', fontsize=16)

    plt.tight_layout()
    plt.show()
    plt.close()

    fig.savefig(fig_location)

    return subgraphs


def save_mask(df, mask_name, columns: [None, list] = None, rows: [None, list] = None, data_dir='data'):
    if rows is None and columns is not None:
        mask = (df == df) & (~df.columns.isin(columns))
    elif rows is not None and columns is None:
        mask = pd.DataFrame().reindex_like(df)
        for col in df.columns:
            mask[col] = rows
    else:
        raise KeyError('Columns or mask must be specified')

    mask_path = os.path.join(data_dir, f'masks/{mask_name}.csv')
    mask.to_csv(mask_path, index=False)


def check_correlated_column(df, threshold=0.9, remove=False, preserve_columns=None, plot_dir='plot',
                            data_dir='data') -> pd.DataFrame:
    if preserve_columns is None:
        preserve_columns = []
    graph_location = os.path.join(plot_dir, f'correlations_grapgs_{threshold}.jpeg')
    heatmap_location = os.path.join(plot_dir, f'correlations_heatmap_{threshold}.jpeg')

    correlations = get_correlated_descriptors(df, threshold, heatmap_location.format(threshold))
    graphs = create_correlation_graphs(correlations, graph_location.format(threshold), threshold)
    duplicit_correlated = get_duplicit_correlated_descriptors(graphs, preserve_columns)

    if remove:
        df.drop(duplicit_correlated, axis=1, inplace=True)
        logging.info(f'Removing correlated columns: {duplicit_correlated}')
    else:
        save_mask(df, mask_name='correlations_' + str(threshold), columns=duplicit_correlated, data_dir=data_dir)
    return df


def remove_duplicits(df: pd.DataFrame, subset, keep='first'):
    duplicits = df[df.duplicated(subset=subset, keep=keep)]
    duplicits_smiles = duplicits['Smiles'].to_list()

    for smile in duplicits_smiles:
        phototox = df.loc[df["Smiles"] == smile, 'Phototoxic']
        if phototox.all() != phototox.any():
            logging.error(
                f'There are two duplicit rows with same smiles, but diffrent phototox value. Removing one of them.'
                f'Names: {df.loc[df["Smiles"] == smile, "Name"]}')

    logging.warning(f"Found duplicates:{duplicits[subset].to_list()} deleted.")
    df.drop_duplicates(subset=subset, keep=keep, inplace=True)

    return df


def check_outliers(df, threshold=4.2, remove=False):
    numerics = ['int16', 'int32', 'int64', 'float16', 'float32', 'float64']
    df_numeric = df.select_dtypes(include=numerics)

    outliers = df[~(np.abs(stats.zscore(df_numeric)) < threshold).all(axis=1)]['Name'].unique()
    logging.info(f'Outliners: {outliers}')

    rows = ~df['Name'].isin(outliers)

    if remove:
        df = df[rows]
        logging.info(f'Outliners {outliers} deleted.')
    else:
        save_mask(df, mask_name='outliers_' + str(threshold), rows=rows)

    return df


def main(input_file='data/chem_output/chem_populated.csv', output_file='data/final/phototox.csv',
         project_dir='.', correlation_threshold=0.95, outliers_threshold=4.2, preserve_columns=None,
         remove=False):
    if preserve_columns is None:
        preserve_columns = []

    df = load_data(input_file)
    df = remove_useless_columns(df)
    df = set_proper_data_type(df)

    plot_dir = os.path.join(project_dir, 'plot')
    data_dir = os.path.join(project_dir, 'data')

    df = check_correlated_column(df, threshold=correlation_threshold, remove=remove, preserve_columns=preserve_columns,
                                 plot_dir=plot_dir, data_dir=data_dir)
    df = remove_duplicits(df, subset="Smiles", keep='first')
    df = check_outliers(df, threshold=outliers_threshold, remove=remove)

    df.to_csv(output_file, index=False)


def process_config(config_file='conf/cleaner.ini') -> dict:
    """
    Process input config and extract
    :param config_file: config file name location
    :return: confuguration
    """
    if not os.path.exists(config_file):
        logging.warning(f'Config file {config_file} does not exist. Running with default setting.')
        return dict()

    config = configparser.ConfigParser()
    config.read(config_file)
    default = config['DEFAULT']

    method = default['method']
    method_dict = {'mask': False, 'remove': True}
    try:
        method = method_dict[method]
    except KeyError:
        raise ValueError('Bad value in ini file. Method can be eather "mask" of "remove".')

    return {'correlation_threshold': float(default['correlation_threshold']),
            'outliers_threshold': float(default['outliers_threshold']),
            'preserve_columns': default['preserve_columns'],
            'remove': method}


if __name__ == "__main__":
    script_dir = os.path.dirname(os.path.realpath(__file__))
    project_dir = os.path.dirname(script_dir)

    args = process_config(os.path.join(project_dir, 'conf/cleaner.ini'))
    input_data = os.path.join(project_dir, 'data/chem_output/chem_populated.csv')
    output_data = os.path.join(project_dir, 'data/chem_output/phototox.csv')
    main(input_data, output_data, project_dir=project_dir, **args)
