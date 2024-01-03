#Author: Linus Lind Jan. 2024
#GNU General Public License v3.0

import sys
import os 
import numpy as np
import pandas as pd
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'scripts'))
from fps.fingerprints import MACCSAnalysis

################################################################################
#Code from: maccs_stats_explanations.ipynb by Hilda Sandström
data = pd.read_csv('../maccs_improvement/geckoq_smiles_molecules.csv')
maccs_data = MACCSAnalysis(data)
percentages  = np.linspace(0,1,5)
feature_counts, nz_features = maccs_data.get_counts()
shared_features = maccs_data.shared_features(nz_features, percentages)
maccs_data.plot_feature_distribution(shared_features, nz_features, percentages, '../maccs_improvement/reports/figures/test')
maccs_data.plot_feature_fp(feature_counts, '../maccs_improvement/reports/figures/test')
nz_keys = pd.DataFrame(nz_features, columns=['in X percent of dataset'])
nz_keys.index.name = 'key'
key_def = pd.read_csv('../maccs_improvement/scripts/fps/explained_patterns.csv')
nz_keys['Explanation'] = key_def['Question'].iloc[nz_keys.index]
nz_keys.to_csv('../maccs_improvement/reports/nz_keys_geckoq.csv')
################################################################################

print(nz_keys)
all_keys = list(np.linspace(0, 166, 167, dtype=int))
used_keys = list(nz_keys.index)

def Diff(A, B): #calculates the set difference A - B
    li_dif = [i for i in A + B if i not in A or i not in B]
    return li_dif

unused_keys = Diff(all_keys, used_keys)
unused_keys = pd.DataFrame(key_def['Question'].loc[unused_keys], columns=['Question'])
unused_keys.index.name = 'key'

print(unused_keys)
print(f'Unused keys: {len(unused_keys)}')
print(f'Used keys: {len(nz_keys)}')
unused_keys.to_csv('../maccs_improvement/reports/unused_keys_geckoq.csv')

feature_percentages = feature_counts / 3414
feature_percentages = feature_percentages[feature_percentages <= 0.001]
low_usage_keys = pd.DataFrame()
low_usage_keys['Explanation'] = key_def['Question'].iloc[feature_percentages.index]
low_usage_keys.to_csv('../maccs_improvement/reports/low_usage_keys_geckoq.csv')

feature_counts.to_csv('../maccs_improvement/reports/feature_counts_geckoq.csv')

print(low_usage_keys)
print(f'Unused or low usage keys: {len(low_usage_keys)}')