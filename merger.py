import pandas as pd

with open('./data/new_pubchem.csv', 'r') as pubfile, open('./data/new_pubchem_clean.csv', 'w') as cleanpubfile:
    lines = pubfile.readlines()
    
    for line in lines:
        new_line = ''.join(letter for letter in line if ord(letter) in range(127))
        cleanpubfile.writelines(new_line)


df_swiss = pd.read_csv('./data/new_swiss.csv')
df_pubchem = pd.read_csv('./data/new_pubchem_clean.csv', index_col=0)

df_concat = pd.concat([df_swiss, df_pubchem], axis=1)

df_concat.to_csv('./merger/merge.csv', index=False, header=True)