'''
export PATH=/mnt/home/azodichr/miniconda3/bin:$PATH
### should ssh dev-intel16-k80

input1: expression data, Fold change or FPKM
input2: start
input3: stop
'''
import os,sys
import pandas as pd
import numpy as np
import random
from scipy.stats.stats import pearsonr
#from scipy.stats import spearmanr
from scipy import stats
#from scipy.stats import chisqprob
import itertools
import math
from sklearn.metrics.cluster import normalized_mutual_info_score

file = sys.argv[1]
#pair = open(sys.argv[2],'r').readlines()
start = int(sys.argv[2])
stop = int(sys.argv[3])
out = open('/mnt/scratch/peipeiw/MI_EC/MI_%s_%s_%s'%(file,start,stop),'w')

df = pd.read_csv('/mnt/home/peipeiw/Documents/Pathway_prediction/20180827_all_EC_pathway/Expression_for_EC_genes/'+file, sep='\t', index_col = 0, header = 0)
D = {} ###
rowname = df.index.tolist()
title = 'gene'
for name in rowname:
	title = title + '\t' + name
out.write(title + '\n')
out.flush()

x = start -1
while x < stop:
	gene1 = rowname[x]
	result = gene1
	for gene2 in rowname:
		MI = float(normalized_mutual_info_score(df.loc[gene1,:],df.loc[gene2,:]))
		result = result + '\t%s'%MI
	out.write(result + '\n')
	out.flush()
	x += 1

out.close()

