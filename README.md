# Gene_coexpression_scripts
Once you obtain your gene expression matrix
1) Calculate random background gene co-expression


2) Calculate pathway EC and random pathway ECs

      EC=expression coherence=(# of gene pairs with PCC > PCC95)/total # of gene pairs

      For a pathway with n genes, total number of gene pairs is taken as: [n*(n-1)]/2, without the self-pairs 

3) Clustering

      kmeans, hclust (ward, complete and average linkages), cmeans, akkmeans, WGCNA
      
4) Visualize clusters

      1. Normalize expression matrix: all values are normalized from 0 to 1 per gene
      
                  python normalization.py <expression matrix> <row or col>
  
      2. Combine the expression matrix with each cluster
      
                  python combine_exressionmatrix.py <cluster file> <normalized expression file>
  
      3. Get visualized expression cluster
      
          input is the output from step 2. This script is meant to run locally on your computer.
          
                  coexpression_profile_from_cluster_plot_loop.R
      
      

