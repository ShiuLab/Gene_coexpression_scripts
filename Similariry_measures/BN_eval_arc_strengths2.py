import os, sys
from itertools import combinations
direc = sys.argv[1]

output = open("bayesian_random_ec", "w")
for file in os.listdir(direc):
    if file.startswith("arc_strength_pvalueMI"):
        pathway1 = file.strip().split("_")[3]
        pathway2 = file.strip().split("_")[4]
        dict = {}
        liste = []
        #pathway2 = pathway.replace("-", "_")
        if pathway2 != "100":
            file2 = "selected_gene_expression_%s_%s" %(pathway1, pathway2)
        else:
            file2 = "selected_gene_expression_%s" %(pathway1)
        file2_ = open(direc+file2, "r")
        line = file2_.readline()
        line = file2_.readline()
        while line:
            info = line.strip().split("\t")
            gene = info[0]
            liste.append(gene)
            line = file2_.readline()
        leng_genez = len(liste)
        file2_.close()
        file1_ = open(direc+file, "r")
        fop = file1_.readline()
        fop = file1_.readline()
        while fop:
            inf = fop.strip().split("\t")
            g1 = inf[1].strip('"')
            g2 = inf[2].strip('"')
            g3 = g1+"_"+g2
            p = float(inf[3].strip('"'))
            #print g3
            print p
            dict[g3] = p
            fop = file1_.readline()
        file1_.close()
        count = 0
        for combo in combinations(liste, 2):
            gen1=list(combo)[0]
            gen2=list(combo)[1]
            possible1 = gen1+"_"+gen2
            possible2 = gen2+"_"+gen1
            if possible1 in dict:
                #if p <0.05:
                count =count+1
            elif possible2 in dict:
                #if p <0.05:
                count =count+1
            else:
                print "nothing"
        ec2 = float(count)/float((leng_genez*(leng_genez-1)/2))
        output.write("%s_%s\t%f\n" % (pathway1, pathway2, ec2))
        
        
        
        
    
#arc_strength_pvalueMI
#selected
