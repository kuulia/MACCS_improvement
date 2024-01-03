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

################################################################################
#Code from: simpol-groups_demo-test.ipynb by Hilda Sandström
groups = pd.read_csv('../maccs_improvement/scripts/aprl_ssp/SMARTSpatterns/SIMPOLgroups_SANE.csv').set_index('substructure')
data = pd.read_csv('../maccs_improvement/all_smiles_molecules.csv')
output_file = '../maccs_improvement/data/output/all_smiles_simpol_groups.csv'
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
#Export
output.to_csv(output_file, index_label='compound')
print(output)