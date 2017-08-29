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
    #print z
    
threshold2 = scoreatpercentile(mi, 95)
print threshold2
output = open('EC_mi_adj_random_%s' % type, 'w')

output.write("Pathway\tEC95pcc_%f\t#ofgenes\t#ofgenesinexpress\n" %(threshold2))


dict_ec = {}
dict_ec2 ={}

#output = open('%s_EC_pc_random' % type, 'w')
#output.write("Pathway\tEC95pcc_%f\tEC95_random\t#ofgenes\t#ofgenesinexpress\n" %(threshold2))
#pcc for pathways
#dict_pcc_path ={}

#dict_ec = {}
dict_ec2 =[]
for sayi in range(0,100):
    for shr in t:
        rand_gen_list=[]
        pcc_random = []
        count2_random = 0
        if shr <3:
            pass
        else:
            for fcv in range(1, shr+1):
                rand1 = random.choice(gen_list)
                if rand1 in rand_gen_list:
                    rand1 = random.choice(gen_list)
                    rand_gen_list.append(rand1)
                else:
                    rand_gen_list.append(rand1)
            for combo2 in combinations(rand_gen_list, 2):
                gen_rand1=list(combo2)[0]
                gen_rand2=list(combo2)[1]
                gen1exprand = dict_e[gen_rand1]
                gen2exprand = dict_e[gen_rand2]
                pcc_rand = adjusted_mutual_info_score(gen1exprand, gen2exprand)
                pcc_random.append(pcc_rand)
            for p2 in pcc_random:
                if p2 > threshold2:
                    count2_random = count2_random + 1
            leng_genez = len(rand_gen_list)
            ec2_random = float(count2_random)/float((leng_genez*(leng_genez-1)/2))
            print ec2_random
            output.write("%i\t%f\n" % (shr, ec2_random))
output.close()



