#Author: Linus Lind Jan. 2024
#GNU General Public License v3.0
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

data = pd.read_csv('data/validation/all_smiles_simpol_fp.csv', index_col=[0])

#print(data)

print(np.sum(data['(alkane CH)_4plus']))
print(np.sum(data['(carbonyl)_4plus']))
print(np.sum(data['(ketone)_4plus']))

sums = np.sum(data / len(data))
sums = sums

plt.figure(figsize=(6, 10))
plt.barh(sums.index, sums.values, height=0.6)
plt.subplots_adjust(left=0.5)
plt.grid(True, axis='x')
plt.yticks(fontsize=8)
#fig, ax = plt.subplots()
#ax.bar(sums.values, sums.index, width=0.3)
#ax.invert_yaxis()
plt.show()