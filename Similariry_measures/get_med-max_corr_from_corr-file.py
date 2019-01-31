'''This script will get the max and median values of each gene to a list of given genes given 
an expression correlation matrix'''

import os, sys
import numpy as np

corrmatrix = open(sys.argv[1], 'r') #expression correlation matrix- gene X gene correlation
genepath_file = open(sys.argv[2], 'r') #file with gene:path or gene:class
name= str(sys.argv[1]).strip().split('/')[1]
output= open(str(name)+'-'+str(sys.argv[2])+"med-max.txt","w")

def get_path_genes(inp1, D):
    gene_path_list=[]
    head= inp1.readline()
    for line in inp1:
        L = line.strip().split('\t')
        gene = L[0]
        path_obj = L[1]
        if path_obj not in D:
            D[path_obj]= [gene]
        else:
            D[path_obj].append(gene)
        if gene not in gene_path_list:
            gene_path_list.append(gene)
        else:
            pass
    return(D, gene_path_list)
    
D={}
path_dict, gene_path_list = get_path_genes(genepath_file, D)
genepath_file.close()
#print(path_dict)

def get_gene_pcc(inp2, D2):
    genes_header= corrmatrix.readline()
    genes_header1= genes_header.strip().split('\t')
    for line in inp2:
        L = line.strip().split('\t')
        if len(L) > 1:
            gene = L[0]
            data = L[1:]
            if gene not in D2:
                D2[gene]=data
            else:
                print(gene, "duplicate- check")
        else:
            pass
    return(D2, genes_header1)

D2={}
print("Getting pcc dictionary of genes")
pcc_dict, genes_header= get_gene_pcc(corrmatrix, D2)
corrmatrix.close()

def write_path_pcc(path_dict, pcc_dict, genes_header):
    path_list= sorted(path_dict.keys())
    newpath_list=[]
    for p in path_list:
        pmed= str(p)+'_med'
        pmax= str(p)+'_max'
        newpath_list.append(pmed)
        newpath_list.append(pmax)
    pathstr= "\t".join(newpath_list)
    output.write("gene\t%s\n" % pathstr)
    for gene in pcc_dict:
        output.write("%s\t" % gene)
        data= pcc_dict[gene]
        for path in path_list:
            genes= path_dict[path]
            data_list= []
            for gene1 in genes:
                if gene1 == gene:
                    pass
                else:
                    if gene1 in genes_header:
                        ind= genes_header.index(gene1)
                        #print(ind)
                        #print(len(data))
                        x=data[ind]
                        if x == '': #check if completely blank
                            pass
                        else:
                            if x != "NA" and (not x.isspace()): #don't get whitespace characters
                                if float(x) != float(1.0):
                                    data_list.append(float(x))
                                else:
                                    pass
                            else:
                                pass
                    else:
                        pass
            if not data_list:
                path_med= 'NA'
                path_max= 'NA'
            else:
                #print(data_list)
                path_med= np.median(data_list)
                path_max= max(data_list)
                print(path_max)
            output.write("%s\t%s\t" % (path_med, path_max))
        output.write("\n")
    output.close()

print("writing median and maximum for each gene for each pathway")    
write_path_pcc(path_dict, pcc_dict, genes_header)
          