args = commandArgs(TRUE)
exprs = args[1]
#distance = args[2]
#n = args[3]
#linkage = args[4]
path = args[2]
path_output = args[3]

library('amap')
library('dynamicTreeCut')
setwd(path)

exprs = as.matrix(exprs)
exprs = read.table(exprs, header=T,row.names=1,sep="\t")

d <- Dist(as.matrix(exprs), method="euclidean")
tree.euclidian <- hclust(d)
cut2 <-cutreeDynamic(tree.euclidian, distM = as.matrix(d), deepSplit=2) 
#clusters.euclidian50 <- cutree(tree.euclidian,k=n)
list1 <- as.list(cut2)
write.table(as.matrix(list1), path_output, sep="\t") 
