#use R cor2pcor for all pathways
#input directory pathways
#expression dataset

import os,sys
import numpy as np
import random
#from scipy.stats import scoreatpercentile

def get_pcor(file, inf, path):
    tmpR = open("%s_tmp.R" % inf,"w")
    tmpR.write("library('corpcor')\n")
    tmpR.write("setwd(%r)\n" % ("%s"  % path))
    tmpR.write("examp = read.table('%s', row.names=1)\n" % (file))
    tmpR.write("pcr1 = pcor.shrink(t(examp))\n")
    tmpR.write("write.table(pcr1, 'PC_%s', sep=',')"%inf)
    tmpR.close()
    os.system("R CMD BATCH %s_tmp.R" % file)

current = "."
expression = sys.argv[1] #whole path
expression_open = open(expression, "r")
line1 = expression_open.readline()
line1 = expression_open.readline()
dict_e = {}
gen_list = []
while line1:
    info = line1.strip().split()
    gen_list.append(info[0])
    dict_e[info[0]] = [float(i) for i in info[1:]]
    line1 = expression_open.readline()
#creating random pathways
pc =[]
for i in range(1, 1000):
    output = open("random_genes_%i" %i, "w")
    for j in range(1, 4):
        a = random.choice(gen_list)
        x = dict_e[a]
        output.write("%s\t" % a)
#       y = "\t".join(x)
        for y in x:
            output.write("%f\t" % y)
        output.write("\n")
    output.close()
    file = current+"/random_genes_%i" % i
    info = "random_genes_%i" % i
    mat1= get_pcor(file, info, current)
    if "PC_%s" %info in os.listdir(current):
        op1 = open("PC_%s" %info , "r")
        pc_temp = []
        line = op1.readline()
        line = op1.readline()
        while line:
            toparse = line.strip().split(",")
            for u in toparse[1:]:
                if u!="1":
                    if float(u) not in pc_temp:
                        pc_temp.append(float(u))
            line = op1.readline()
        for o in pc_temp:
            pc.append(o)
        op1.close()
output = open("random_PCs", "w")
for p in pc:
    output.write("%f\n" %p)
output.close()
"""
directory_pathways = sys.argv[2]
liste = []
for file in os.listdir(directory_pathways):
    if file.startswith("selected"):
        info = file.split("_")
        toput = directory_pathways+file
        mat = get_pcor(directory_pathways+file, info,current)
        
        
        os.system("rm %s_tmp.R" % info)
                       
"""        
        
        




