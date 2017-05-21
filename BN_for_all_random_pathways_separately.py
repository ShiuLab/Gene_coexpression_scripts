import random
import sys, os
#from scipy.stats import pearsonr
from scipy.stats import spearmanr
from scipy.stats import scoreatpercentile
from itertools import combinations

#aracyc_arrange = convert aracyc file into dictionary; dict[pathway] = [gene1, gene2, gene3....]
#def aracyc_arrange():
aracyc = sys.argv[1] #directory
expression = sys.argv[2]
output_direc = sys.argv[3]
dict = {}
t = []
for file in os.listdir(aracyc):
    if file.startswith("genelist_exp"):
        inf = file.strip().split("_")
        p = inf[2]
        file1 = aracyc+file
        a = open(file1, "r")
        line = a.readline()
        count = 0
        while line:
            count = count+1
            inf = line.strip().split()
            if p not in dict:
                dict[p] = [inf[0]]
            else:
                dict[p].append(inf[0])
            line = a.readline()
        t.append(count)
        a.close()
#print dict

expression_open = open(expression, 'r')
line1 = expression_open.readline()
line1 = expression_open.readline()
dict_e = {}
gen_list = []
while line1:
    info = line1.strip().split()
    gen_list.append(info[0])
    dict_e[info[0]] = [float(i) for i in info[1:]]
    line1 = expression_open.readline()
#outpcc = open('random_pcc_distribution_%s' % type , 'w')    
#random genes, calculate pcc
#output_pcc = open("pcc_random_%s" % type, "w")
#output = open('5_EC_sp_random', 'w')
#output.write("Pathway\tEC95pcc_%f\tEC95_random\t#ofgenes\t#ofgenesinexpress\n" %(threshold2))
#pcc for pathways
#dict_pcc_path ={}
for p in dict:
    for sayi in range(0,100):
        rand_gen_list= dict[p]
        if shr <3:
            pass
        else:
            output = open("random_genes_pluspathways_%i_%i" %(shr, sayi),"w")
            for fcv in range(1, shr+1):
                rand1 = random.choice(gen_list)
                if rand1 in rand_gen_list:
                    rand1 = random.choice(gen_list)
                    output.write("%s\n" %rand1)
                else:
                    rand_gen_list.append(rand1)
                    output.write("%s\n" %rand1)
            output.close()
            