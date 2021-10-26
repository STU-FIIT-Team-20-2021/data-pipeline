import pandas as pd
import argparse
import sys
import os

parser = argparse.ArgumentParser(description='Accepting pubchem and swiss files.')
parser.add_argument("-pub", type=str, help="Pubchem file to be merged.")
parser.add_argument("-swiss", type=str, help="Swiss file to be merged.")
args = parser.parse_args()

file_pubchem = './data/new_pubchem.csv'
file_swiss = './data/new_swiss.csv'

if args.pub:
    file_pubchem = args.pub

if args.swiss:
    file_swiss = args.swiss

with open(file_pubchem, 'r') as pubfile, open('./data/new_pubchem_clean.csv', 'w') as cleanpubfile:
    lines = pubfile.readlines()
    
    for line in lines:
        new_line = ''.join(letter for letter in line if ord(letter) in range(127))
        cleanpubfile.writelines(new_line)


df_swiss = pd.read_csv(file_swiss)
df_pubchem = pd.read_csv('./data/new_pubchem_clean.csv', index_col=0)

df_concat = pd.concat([df_pubchem, df_swiss], axis=1)

df_concat.to_csv('./merger/merge.csv', index=False, header=True)

if not os.name == 'nt':
    df_concat.to_csv(sys.stdout, index=False, header=True)
