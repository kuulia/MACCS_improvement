#Author: Linus Lind Jan. 2024
#GNU General Public License v3.0

import pandas as pd
import numpy as np

raw_data = pd.read_csv('../maccs_improvement/data/output/all_smiles_simpol_all_groups.csv')

df = raw_data.drop(columns=['compound']) #drop compound column to calculate max and min

maxes = df.max(axis=0)
idx_maxes = df.idxmax(axis=0)
compound_maxes = raw_data.iloc[:,0]
compound_maxes = compound_maxes.iloc[idx_maxes]

means = df.mean(axis=0)

mins = df.min(axis=0)
idx_mins = df.idxmin(axis=0)
compound_mins = raw_data.iloc[:,0]
compound_mins = compound_mins.iloc[idx_mins]

data = pd.DataFrame()
data['Min'] = mins
data['Mean'] = means
data['Max'] = maxes
data['Min_index'] = idx_mins
data['Max_index'] = idx_maxes

def get_counts(
        df: pd.DataFrame, 
        comp: str) -> int:
    df_pruned = df[ df[comp] != 0]
    return df_pruned.shape[0]

def counts_of_groups(df: pd.DataFrame) -> list:
    group_counts = list()
    groups = df.columns
    for group in groups:
        group_counts.append(get_counts(df, group))
    return group_counts

def weighted_mean(df: pd.DataFrame) -> list:
    group_weighted_means = list()
    groups = df.columns
    for group in groups:
        counts = get_counts(df, group)
        if (counts != 0):
            sums = df[group].sum()
            group_weighted_means.append(sums / counts)
        else: group_weighted_means.append(0)
    return group_weighted_means

#def weighted_mean(
#        df: pd.DataFrame, 
#        comp: str) -> int:
#    df_pruned = df[ df[comp] != 0]
#    sum = df_pruned.sum()
#    counts = get_counts(df, comp)
#    return sum / counts

#pats = pd.read_csv('../maccs_improvement/scripts/aprl_ssp/SMARTSpatterns/SIMPOLgroups_noring_nonnitrophenol.csv')
#print(len(data))
#print(len(pats))
#data['SMARTSpattern'] = pats['pattern'].astype(str)
data['Counts'] = counts_of_groups(df)
data['Weighted_mean'] = weighted_mean(df)
data.index.name = 'Compound'
data.to_csv('../maccs_improvement/data/output/mins_and_maxes.csv')
print(data)

pruned_data = data[ data['Max'] != 0 ]

print(pruned_data)
print(f'Pruning results in {len(pruned_data)} potentially useful keys ' \
      + f'out of the total {len(raw_data.columns)} simpol groups')

pruned_data.to_csv('../maccs_improvement/data/output/potential_simpol_groups.csv')

further_pruned_data = data[ data['Max'] > 1 ]
further_pruned_data.to_csv('../maccs_improvement/data/output/multiple_simpol_groups.csv')
print(f'\n{further_pruned_data}')

print(f'There are {len(further_pruned_data)} groups that appear more than once in the data set ' \
      + f'out of the pruned data')

more_than_four = data[ data['Max'] > 4 ]
print(more_than_four)
print(f'There are {len(more_than_four)} groups that can appear more than four times in the data set')
more_than_four.to_csv('../maccs_improvement/data/output/more_than_four_simpol_groups.csv')