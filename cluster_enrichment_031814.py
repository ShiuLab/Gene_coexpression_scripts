import sys, os
cluster_result = sys.argv[1]
go_annot = sys.argv[2]

#cluster
cluster_result_open = open(cluster_result, "r")
line = cluster_result_open.readline() #header
line = cluster_result_open.readline()
expre_gen =[]
dict_cluster = {}
while line:
    info = line.strip().split()
    gene = info[0].replace('"', '')
    cluster = info[1]
    expre_gen.append(gene)
    if cluster in dict_cluster:
        dict_cluster[cluster].append(gene)
    else:
        dict_cluster[cluster] = [gene]
    line = cluster_result_open.readline()
#print len(expre_gen) #this will be the total genes used in contingency table
print gene
#AraCyc
go_gene_open = open(go_annot, "r")
dict = {}
line1 = go_gene_open.readline()
line1 = go_gene_open.readline()
genler_list = []
while line1:
    info = line1.strip().split("\t")
    pathway = info[0]
    for i in info:
        if i.startswith("AT"):
            #gene = i.split(".")
            genes = i.split(".")
            gene = genes[0].split("-")
            gene = gene[0]
            if gene.startswith("AT"):
                gene1 = gene
                if pathway in dict:
                    if gene1 not in dict[pathway]:
                        dict[pathway].append(gene1)
                else:
                    dict[pathway] = [gene1]
    line1 = go_gene_open.readline()

output_table = open("tableforEnrichment", "w")
for item in dict:
    gene_list_for_go = dict[item]
    for i in dict_cluster.keys():
        cluster_list = dict_cluster[i]
        count1 = 0
        count2 = 0
        count3 = 0
        count4 = 0
        for j in cluster_list:
            if j in gene_list_for_go:
                count1 = count1 + 1
            else:
                count2 = count2 + 1
        count3 = len(gene_list_for_go) - count1
        count4 = 20060 - (count1 + count2 + count3)
        output_table.write("%s_%s\t%i\t%i\t%i\t%i\n" % (item, i, count1, count2, count3, count4))
output_table.close()