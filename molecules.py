#Author: Linus Lind Jan. 2024
#GNU General Public License v3.0

import pandas as pd
import numpy as np
import csv

df = pd.read_csv('all_smiles_mols.csv') #generated using openbabel report
df.columns = ['column']
df = df[ df[ 'column' ].str.contains('FORMULA: ')==True ] #We just need the molecule formula
df = df.replace('FORMULA: ', '', regex=True)
df = df.reset_index()
df = df.drop(['index'],axis='columns')
print(df)

smiles = pd.read_csv('all_smiles.csv', header=None)
smiles.columns = ['SMILES']
print(smiles)
output = pd.DataFrame()
output.insert(0, 'compound', df['column'])
output.insert(1, 'SMILES', smiles['SMILES'])
output.SMILES = output.SMILES.astype(str)
output.compound = output.compound.astype(str)
print(output)
output.to_csv('all_smiles_molecules.csv', index = False, quotechar='"',quoting=csv.QUOTE_NONNUMERIC)