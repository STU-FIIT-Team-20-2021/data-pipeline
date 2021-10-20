import pandas as pd
import numpy as np

with open('./data/columns.dat', 'r') as columnfile:
    pubchem_columns = columnfile.readline().strip().split(',')
    swiss_columns = columnfile.readline().strip().split(',')

df_merge = pd.read_csv('./data/merged.csv')
df_swiss = pd.read_csv('./data/new_swiss.csv')
df_pubchem = pd.read_csv('./data/new_pubchem.csv', index_col=0)

df_swiss = df_swiss[swiss_columns]
df_pubchem = df_pubchem[pubchem_columns]

df_concat = pd.concat([df_swiss, df_pubchem], axis=1)

df_merge = pd.concat([df_merge, df_concat], axis=0)

df_merge.to_csv('./data/new_merged.csv', index=False)