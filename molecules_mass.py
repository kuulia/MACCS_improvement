#Author: Linus Lind Jan. 2024
#GNU General Public License v3.0

import pandas as pd
import numpy as np
import csv

df = pd.read_csv('all_smiles_mols.csv') #generated using openbabel report
df.columns = ['column']
df = df[ df[ 'column' ].str.contains('MASS: ')==True ] #We just need the molecule formula
df = df.iloc[::2, :]
df = df.replace('MASS: ', '', regex=True)
df = df.reset_index()
df = df.drop(['index'],axis='columns')
df = df.astype(float)
df = df
print(df)
print(min(df['column']))
print(max(df['column']))
print(np.mean(df['column']))
smiles = pd.read_csv('all_smiles.csv', header=None)
smiles.columns = ['SMILES']
#print(smiles)
output = pd.DataFrame()
output.insert(0, 'SMILES', smiles['SMILES'])
output.insert(1, 'molar_mass', df['column'])
output['log_p_sat'] = pd.read_clipboard(header=None)
print(output)
output.to_csv('all_smiles_molar_mass.csv', index = False)