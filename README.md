# Gene_coexpression_scripts
Once you obtain your gene expression matrix
1) Calculate random background gene co-expression


2) Calculate pathway EC and random pathway ECs

#EC=expression coherence=(# of gene pairs with PCC > PCC95)/total # of gene pairs

#For a pathway with n genes, total number of gene pairs is taken as: [n*(n-1)]/2, without the self-pairs 

3) Clustering

#kmeans, hclust (ward, complete and average linkages), cmeans, akkmeans, WGCNA
