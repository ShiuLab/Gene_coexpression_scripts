'''
export PATH=/mnt/home/azodichr/miniconda3/bin:$PATH
### should ssh dev-intel16-k80

## when submitting to queue, need the working directory where your file is
## change the scratch address to your scratch for output
'''
import os,sys
import pandas as pd
import numpy as np
import random
from scipy.stats.stats import pearsonr
from scipy import stats
from scipy.stats import chisqprob

file = sys.argv[1]
SAVE = '/mnt/scratch/john3784/PCC/PCC_' + file

df = pd.read_csv(file, sep='\t', index_col = 0, header = 0)  
pcc = df.T.corr(method='pearson').round(3)
pcc.to_csv(SAVE, index=True, header=True,sep="\t")


