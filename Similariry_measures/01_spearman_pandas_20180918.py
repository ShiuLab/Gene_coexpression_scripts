'''
export PATH=/mnt/home/azodichr/miniconda3/bin:$PATH
### should ssh dev-intel16-k80
'''
import os,sys
import pandas as pd
import numpy as np
import random
from scipy.stats.stats import pearsonr
from scipy.stats.stats import pearsonr
from scipy import stats
from scipy.stats import chisqprob

file = sys.argv[1]
SAVE = '/mnt/scratch/peipeiw/Spearman_EC/Spearman_' + file

df = pd.read_csv('/mnt/home/peipeiw/Documents/Pathway_prediction/20180827_all_EC_pathway/Expression_for_EC_genes/'+file, sep='\t', index_col = 0, header = 0)  
pcc = df.T.corr(method='spearman').round(3)
pcc.to_csv(SAVE, index=True, header=True,sep="\t")


