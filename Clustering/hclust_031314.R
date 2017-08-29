args = commandArgs(TRUE)
exprs = args[1]
#distance = args[2]
n = args[2]
linkage = args[3]
path = args[4]
path_output = args[5]

library('amap')
setwd(path)

exprs = as.matrix(exprs)
exprs = read.table(exprs, header=T,row.names=1,sep="\t")

d <- Dist(as.matrix(exprs), method="euclidean")
tree.euclidian <- hclust(d, method=linkage)
clusters.euclidian50 <- cutree(tree.euclidian,k=n)
list1 <- as.list(clusters.euclidian50)
write.table(as.matrix(list1), path_output, sep="\t") 









#input <- read.table("/mnt/home/uygunsah/PlastidNetworkProject/Effec_Input/stress_expression_combined_c1_P0.05_FC", header=T, row.names=1, sep="\t")
#d <- dist(input, method="euclidian")
#tree.euclidian <- hclust(d)
#clusters.euclidian50 <- cutree(tree.euclidian,k=100)
#list1 <- as.list(clusters.euclidian50)
#write.table(as.matrix(list1), "/mnt/home/uygunsah/PlastidNetworkProject/Effec_method/hclust200_2", sep=",") 



#dendro<-hclust(as.dist(combo2),method="average")



#cut2<-cutreeDynamicTree(dendro,maxTreeHeight=1,deepSplit=TRUE,minModuleSize=1)

#I then used the following to visualise the data:

#cut2colour<-labels2colors(cut2)
#plotDendroAndColors(dendro,cut2colour,"Dynamic Tree Cut", dendroLabels=FALSE,hang=0.03,addGuide=TRUE,guideHang=0.05)

#hc <- hclust(dist(USArrests))
#labels =  cutree(hc, k=1:5)
#require("WGCNA")
#plotDendroAndColors(hc, labels)
