#Author: Linus Lind Jan. 2024
#GNU General Public License v3.0
import pandas as pd
import numpy as np

#unit conversions for comparing SIMPOL.1 to ML model
data = pd.read_csv('Dataframe.csv')
psats = data['pSat_Pa']
psats_kpa = psats / 1_000
log_psats_kpa = np.log10(psats_kpa)
print(log_psats_kpa)
log_psats_kpa.to_csv('geckoq_log_p_sat.txt', index=False, header=False)