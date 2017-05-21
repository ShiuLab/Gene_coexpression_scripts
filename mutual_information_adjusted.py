import random, math
import sys, os
from sklearn.metrics.cluster import adjusted_mutual_info_score
from scipy.stats import scoreatpercentile
from itertools import combinations

#aracyc_arrange = convert aracyc file into dictionary; dict[pathway] = [gene1, gene2, gene3....]
#def aracyc_arrange():
aracyc = sys.argv[1] #directory
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
expression = sys.argv[2]
type = sys.argv[3]
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
mi = []
for i in range(1, 1000):
    datamatrix = []
    a = random.choice(gen_list)
    b = random.choice(gen_list)
    x = dict_e[a]
    y = dict_e[b]
    z = adjusted_mutual_info_score(x, y)
    mi.append(z)
    print z
    
threshold2 = scoreatpercentile(mi, 95)
print threshold2
output = open('EC_mi_adj_%s' % type, 'w')

output.write("Pathway\tEC95pcc_%f\t#ofgenes\t#ofgenesinexpress\n" %(threshold2))


dict_ec = {}
dict_ec2 ={}

for bio in dict:
    ec =0
    pcc_path_list = []
    new_list =[]
    genez = dict[bio]
    rand_gen_list =[]
    pcc_random =[]
    for el in genez:
        if el in dict_e:
            new_list.append(el)
    leng_genez = float(len(new_list))
    for combo in combinations(new_list, 2):
        gen1=list(combo)[0]
        gen2=list(combo)[1]
        gen1exp = dict_e[gen1]
        gen2exp = dict_e[gen2]
        pcc_path = adjusted_mutual_info_score(gen1exp, gen2exp)
        pcc_path_list.append(pcc_path)
    count2 = 0
    for p in pcc_path_list:
        if p > threshold2:
            count2 = count2 + 1
    if leng_genez > 2:
        ec2 = float(count2)/float((leng_genez*(leng_genez-1)/2))
    else: 
        ec2 = "NA"
    output.write("%s\t%s\t%i\t%i\n" % (bio, str(ec2), len(genez), len(new_list)))
    

expression_open.close()
output.close()




