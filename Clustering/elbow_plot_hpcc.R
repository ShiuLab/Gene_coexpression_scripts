#elbow plot 
#submit on hpcc:
#R --vanilla --slave --args	[start_directory]	[expression_matrix_file]	<	/mnt/home/john3784/3-Solanaceae_project/expr_cluster/elbow_plot_hpcc.R
args = commandArgs(TRUE) #make arguments
setwd(args[1]) #working directory
nc= args[2] #file name (matrix)
#read in data
p1 = read.table(nc, header=T, sep="\t", row.names=1) #will remove NAs later
#dataframe
mydata <- data.frame(p1)
#omit NAs
mydata <- na.omit(mydata)
#get within-group sum of squares
wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
for (i in 5:2000) wss[i] <- sum(kmeans(mydata,
                                     centers=i)$withinss)
#make plot
p.b <-plot(1:2000, wss, type="b", xlab="Number of Clusters",
           ylab="Within groups sum of squares")
#can try plotting smaller range
# wss1 = wss[1:200]
# p.b <-plot(1:200, wss1, type="b", xlab="Number of Clusters",
#            ylab="Within groups sum of squares")
# write to pdf
nb=paste(c(basename(nc), "_elbowplot.pdf"), collapse='')
pdf(file = nb)
p.b
dev.off()
