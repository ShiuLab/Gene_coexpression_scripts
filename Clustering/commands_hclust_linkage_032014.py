#commands_kmeans
#quick_kmeans_for various k
#kmeans.R file should be in the same directory as wd

import sys, os

expression_data = sys.argv[1]
#.rstrip("/")
#startk 
#a = int(sys.argv[2])
#endk 
#b = int(sys.argv[3])
#increment 
#c = int(sys.argv[4])
#distance = sys.argv[2]
#linkage = sys.argv[3]
output_file = sys.argv[2]
path = sys.argv[3]
#run = int(sys.argv[5])
type = sys.argv[4]

oup=open("commands_hclust_%s" % type, "w" )

#b = [100]
b = [250, 300, 400, 500, 750, 1000]
list= ["complete", "average" , "ward"]
for i in b:
    for j in list:
        k = int(i)
        output_file2 = output_file+"h"+str(k)+ "_%s_%s" % (type, j)
        oup.write("R --vanilla --slave --args %s %i %s %s %s < hclust_031314.R\n" %(expression_data, k, j, path, output_file2))
oup.close()

'''
import sys, os

files_path = sys.argv[1]
#.rstrip("/")
code = sys.argv[2]
go = sys.argv[3]
a = sys.argv[4]

oup=open("commands_enrichment_%s" % a ,"w" )

for file in os.listdir(files_path):
    if file.startswith("genelist"):
        file_name = file.strip().split("_")
        kmeans = file_name[1][0]
        type = file_name[2]
        run = int(file_name[3][3])
        file1 = files_path+file
        oup.write("python\t%s\t%s\t%s\t%i\t%s\t%s\n" % (code, file1, go, run, kmeans, type)) 

#python ../goenrichment_cluster.py genelist_k1000_stressint_run2 ~/1_Expression_Database/GO/updated_go_gene_012013 2 k strint

oup.close()
'''
