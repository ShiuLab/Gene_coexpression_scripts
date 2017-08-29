args = commandArgs(TRUE)
exprs = args[1]
n = args[2]
path = args[3]
path_output = args[4]

library('e1071')

setwd(path)

exprs = as.matrix(exprs)
exprs = read.table(exprs, header=T,row.names=1,sep="\t")

n = as.numeric(n)

c1 <- cmeans(exprs, n, iter.max = 100)
list <- as.list(c1$cluster)
vector <- as.vector(c1$cluster)

write.table(as.matrix(list), path_output, sep="\t")