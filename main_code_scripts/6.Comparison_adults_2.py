##############################################################################
# ---- Adding LD block number to reference dataset metabQTLs

import pandas as pd
import numpy as np
import re

regions = pd.read_csv('runCMS-master/others/regions.bed', sep='\t')
df = pd.read_csv('Results_all_papers_urine_combined.txt', sep ='\t')
reg = pd.DataFrame(np.zeros(df.shape[0], dtype = int))

for r in range(0, 1706):
    left = regions[' start '][r]    
    right = regions[' stop'][r]    
    chrnum = [int(i) for i in re.findall('\d+', regions['chr '][r])]   
    ind = np.where((df['CHR'] == chrnum) & (df['POS'] >= str(left) ) & (df['POS'] <= str(right) ))    
    reg.iloc[ind] = r + 1

df['Region'] = reg
df.to_csv(r'./urine_all_studies_with_regions.csv',sep='\t', index=False, header=False)
