# simple kmeans with 100 iterations and 50 random nstart

setwd(path)
exprs = read.table(expression_data, header=T,row.names=1,sep="\t")
c1 <- kmeans(exprs, n, iter.max = 100, nstart = 50)
list <- as.list(c1$cluster)
write.table(as.matrix(list), path_output, sep="\t")


# from class notes:
# K-means
#

# Now let's generate cluster with kmeans. Do this on all genes
K50 <- kmeans(M,50,iter.max=500,nstart=10)

# Take a look at the k-means results
# First we need to install the fpc package
install.packages("fpc")
library(fpc)

# See how well the cluster separate from one another
plotcluster(M,K50$cluster)

# Check out cluster sizes
K50$size

# Get the cluster assignments
K50_assign <- K50$cluster
K50_assign

# Check out which genes in the cluster assignment belong to cluster 1
C1_has <- which(K50_assign==1)
C1_has

# Check size of cluster 1
length(C1_has)

# Get expression matrix for genes in cluster 1
MC1 <- M[C1_has,]
dim(MC1)

# Show expression pattern of genes in cluster 1
# Get the max and min values
MAX <- max(MC1)
MIN <- min(MC1)

# Get the sample names for the plot
CN <- colnames(MC1)
CN

# Number of samples
CNN <- length(CN)

# Plot the first gene
plot(t(MC1[1,]),type="l",ylim=c(MIN,MAX),col="gray",xaxt="n",xlab="")
axis(1,at=seq(0,CNN,length.out=ncol(MC1)),labels=CN,las=2)

# Plot the rest of the genes
for (i in seq(2,dim(MC1)+1)){
  lines(t(MC1[i,]),col="gray")
}

# Plot the median
MED <- sapply(MC1,median)
lines(MED,col="blue")

# Q. Now plot the expression patterns of genes in cluster 2.
