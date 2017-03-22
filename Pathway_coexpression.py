# calculate all pairwise PCCs for all pathways
# module load SciPy
import random
import sys, os
import os.path
from scipy.stats import pearsonr
# module load scipy ! before running this code
# for 3 pathways it took 1min to run this script
from scipy.stats import scoreatpercentile
from numpy import median
from itertools import combinations

print('''
inp1 = Expression dataset, rows are genes and columns are samples, with header (File)
inp2 = Percentile to use for random background, Measure:PCC (Integer:95 for 95th percentile)
inp3 = Pathway annotation first column pathway, second column gene with header(File)
inp4 = Directory to output files (Path)
inp5 = Number of gene pairs to pick randomly for background distribution  (Integer: 500000)
''')

expression = sys.argv[1] #"/mnt/home/uygunsah/projects/1_Expression_Database/expression_matrix/5_Stress_FC"
perc = int(sys.argv[2]) #95
aracyc = sys.argv[3] #aracyc, go, or manual annotation
direc = sys.argv[4] #"/mnt/home/uygunsah/projects/4_Plastid/_June2016_resubmission/PCC_distr/"
pairs = int(sys.argv[5])


#all pairwise-gene pccs
def calculate_pairwise_PCC(list_of_genes):
    pcc_path_list = []
    for combo in combinations(list_of_genes, 2):
        gen1=list(combo)[0]
        gen2=list(combo)[1]
        gen1exp = dict_e[gen1]
        gen2exp = dict_e[gen2]
        pcc_path = pearsonr(gen1exp, gen2exp)[0]
        pcc_path_list.append(pcc_path)
    return pcc_path_list

#ec function
def calculate_EC(pcc_list, threshold, ngenes):
    count = 0
    for p2 in pcc_list:
        if p2 > threshold:
            count = count + 1
    ec2 = float(count)/float((ngenes*(ngenes-1)/2))
    return ec2
    
# expression dictionary
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


# get random pcc distribution and print to file, also print pcc95
oup_random = open("random_pairs_PCC_%s" %pairs, "w")
pcc_list =[]
for i in range(0, pairs):
    a = random.choice(gen_list)
    b = random.choice(gen_list)
    x = dict_e[a]
    y = dict_e[b]
    pcc = pearsonr(x, y)[0]
    oup_random.write("%f\n" %pcc)
    pcc_list.append(pcc)
pcc95 = scoreatpercentile(pcc_list, perc)
#pcc95 = 0.404774
threshold2 = pcc95


# pathway dictionary
aracyc_op = open(aracyc, "r")
line = aracyc_op.readline()
line = aracyc_op.readline()
dict = {}
while line:
    info = line.strip().split("\t")
    pathway = info[0]
    gene = info[1]
    if pathway not in dict:
        dict[pathway] = [gene]
    else:
        if gene not in dict[pathway]:
            dict[pathway].append(gene)
    line = aracyc_op.readline()

# pathway PCCs and EC calculations
output = open("EC_pathways", "w")
all_pathway_genes = []
for bio in dict:
    pcc_path_list = []
    new_list = []
    genez = dict[bio] # genes in pathway
    for el in genez: 
        if el in dict_e: # check if genes in pathway is present in expression dataset
            new_list.append(el)
        if el not in all_pathway_genes:
            all_pathway_genes.append(el)
    leng_genez = float(len(new_list))
    if leng_genez > 2: # only consider pathways with more than 2 genes
        if os.path.isfile(direc+"PCC_distr_%s" % bio): #if the file already present, pass
            print bio
            pass
        else:
            pcc_path_list = calculate_pairwise_PCC(new_list
            #print bio, "writing"
            out = open("PCC_distr_%s" % bio, "w")
            out.write("%s\n%s\n" % (bio, str(med))) # name of pathway and median PCC as header
            for num in pcc_path_list:
                out.write("%s\n" % str(num))
            out.close()                
            ec = calculate_EC(pcc_path_list, threshold2,len(new_list))                
            output.write("%s\t%f\n" % (bio, ec))
output.close()
               

# random ECs, randomly form pathways 50 times for pathway size of 3-100 genes
output = open("EC_random", "w")
ec2 = []
for a in range(1, 50):
    for sayi in range(3,100):
        rand_gen_list=[]
        pcc_random = []
        count2_random = 0
        for fcv in range(1, sayi+1):
            rand1 = random.choice(all_pathway_genes)
            if rand1 in rand_gen_list:
                rand1 = random.choice(all_pathway_genes)
                rand_gen_list.append(rand1)
            else:
                rand_gen_list.append(rand1)
        for combo2 in combinations(rand_gen_list, 2):
            gen_rand1=list(combo2)[0]
            gen_rand2=list(combo2)[1]
            gen1exprand = dict_e[gen_rand1]
            gen2exprand = dict_e[gen_rand2]
            pcc_rand = pearsonr(gen1exprand, gen2exprand)[0]
            pcc_random.append(pcc_rand)
        ec_random = calculate_EC(pcc_random, threshold2,len(rand_gen_list)) 
        output.write("%i\t%f\n" % (sayi, ec_random))
output.close()
