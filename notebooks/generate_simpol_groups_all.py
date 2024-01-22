#Author: Linus Lind Jan. 2024
#GNU General Public License v3.0
#uses APRL-SSP by Satoshi Takahama https://doi.org/10.5281/zenodo.34975

import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'scripts'))
from aprl_ssp.util import searchgroups
from aprl_ssp.substructure_search import count_groups
from openbabel.pybel import readstring
from aprl_ssp.userdef import *

groups = pd.read_csv('../maccs_improvement/scripts/aprl_ssp/SMARTSpatterns/SIMPOLgroups_noring_nonnitrophenol.csv').set_index('substructure')

################################################################################
#Code from: simpol-groups_demo-test.ipynb by Hilda Sandström
data = pd.read_csv('../maccs_improvement/all_smiles_molecules.csv')
output_file = '../maccs_improvement/data/output/all_smiles_simpol_all_groups.csv'
#Define which groups for which stats should be exported to output_file. If none everything will be exported. 
export = None
#Rearrange input to fit the program
inp = data.set_index('compound')
###_* --- Apply search function
search = searchgroups(groups.pattern, export) 
output = count_groups(inp, search)
################################################################################

#add oxygen count
output['oxygen count'] = data['SMILES'].str.count('O').values
# Add ring groups using functions in aprl_ssp userdef
data['mols'] = data['SMILES'].apply(lambda x: readstring('smiles', x))
output['aromatic ring'] = data['mols'].apply(lambda x: count_aromatic_rings(x)).values
output['non-aromatic ring']  = data['mols'].map(count_nonaromatic_rings).values
output['nitrophenol']  = data['mols'].map(lambda x: count_nitrophenols(x, '[OX2H][cX3]:[c]', '[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8] ')).values
###_* --- Export to output
output.to_csv(output_file, index_label='compound')
print(output)