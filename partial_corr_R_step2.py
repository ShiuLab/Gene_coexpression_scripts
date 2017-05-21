#use R cor2pcor for all pathways
#input directory pathways
#expression dataset

import os,sys
import numpy as np
import random, math
#from scipy.stats import scoreatpercentile

def round_sigfigs(num, sig_figs):
    """Round to specified number of sigfigs.

    >>> round_sigfigs(0, sig_figs=4)
    0
    >>> int(round_sigfigs(12345, sig_figs=2))
    12000
    >>> int(round_sigfigs(-12345, sig_figs=2))
    -12000
    >>> int(round_sigfigs(1, sig_figs=2))
    1
    >>> '{0:.3}'.format(round_sigfigs(3.1415, sig_figs=2))
    '3.1'
    >>> '{0:.3}'.format(round_sigfigs(-3.1415, sig_figs=2))
    '-3.1'
    >>> '{0:.5}'.format(round_sigfigs(0.00098765, sig_figs=2))
    '0.00099'
    >>> '{0:.6}'.format(round_sigfigs(0.00098765, sig_figs=3))
    '0.000988'
    """
    if num != 0:
        return round(num, -int(math.floor(math.log10(abs(num))) - (sig_figs - 1)))
    else:
        return 0  # Can't take the log of 0
def get_pcor(file, inf, path):
    tmpR = open("%s_tmp.R" % inf,"w")
    tmpR.write("library('corpcor')\n")
    tmpR.write("setwd(%r)\n" % ("%s"  % path))
    tmpR.write("examp = read.table('%s', row.names=1, header=T)\n" % (file))
    tmpR.write("pcr1 = pcor.shrink(t(examp))\n")
    tmpR.write("write.table(pcr1, 'PC_%s', sep=',')"%inf)
    tmpR.close()
    os.system("R CMD BATCH %s_tmp.R" % inf)

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
threshold = 0.3205512

directory_pathways = sys.argv[2]
liste = []
output_son = open("ecs_R_PC", "w")
for file in os.listdir(directory_pathways):
    if file.startswith("selected"):
        info = file.split("_")[3]
        toput = directory_pathways+file
        mat = get_pcor(directory_pathways+file, info,current)
        if "PC_%s" %info in os.listdir(current):
            op1 = open("PC_%s" %info , "r")
            pc_temp = []
            line = op1.readline()
            numbers= line.strip().split(",")
            ngenes = len(numbers)
            if ngenes>2:
                line = op1.readline()
                while line:
                    toparse = line.strip().split(",")
                    for u in toparse[1:]:
                        if u!="1":
                            uu = round_sigfigs(float(u), 3)
                            if uu not in pc_temp:
                                pc_temp.append(uu)
                    line = op1.readline()
                count = 0
                for v in pc_temp:
                    if v > threshold:
                        count = count+1
                ec = float(count)/float((ngenes*(ngenes-1)/2))
                output_son.write("%s\t%f\n" % (info, ec))
                if info=="FASYN-ELONG-PWY" or info=="PWY-5271":
                    print pc_temp
                    print info, count, ngenes
            else:
                op1.close
output_son.close()
        
        
        #os.system("rm %s_tmp.R" % info)
                       
        
        
        




