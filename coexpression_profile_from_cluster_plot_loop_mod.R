#plot gene coexpression for each cluster and different treatment using log2FC normalization to 0 and 1.

#change lines 19 and 50 based on the number of columns in dataset

setwd("working directory path") #set wd

#input file name with gene:cluster:normalized expression
nc="Col_SH_pinwound_0.25hr_cluster_up_pos_only.txt-wound-stress-hormone_k100_nocbt_run7-sum_contrasts_horm_wound_abiotic_biotic_logFC_P0.05_gene_matrix_nodups.txt_norm.txt"

#read in data
p1 = read.table(nc, header=T, row.names=1, sep="\t") 

#input the AucRoc and F1 values using Random Forest for each cluster.
#p2 = read.table("Slyc_FC-FPKM_k10_M82_up_homologexpr_run7-Results_median_FC_for_all_Slyc_20180125.txt_norm.txt", header=T, row.names=1, sep="\t")
#head(p2)

#get clusters
cluster = unique(p1[,1])

#output name
nb=paste(c(basename(nc), "_expr_profile.pdf"), collapse='')

#pdf
pdf(file =nb)
par(mfrow=c(3,2))
e = c()
e1 = c()
e2 = c()
#loop for each cluster to visulalize in pdf
for(i in cluster){
	e = subset(p1,p1[,1]==i)
	e1 = e[,2:93] #number of colomns
	e2 = t(e1)
	#AucRoc = p2[i,2]
	#F1 = p2[i,3]
	plot(e2[,1], ylim=c(0,1), main = paste("cluster",i,", No.",nrow(e),"\n",sep=" "), 
	xlab = "", ylab="log2FC", type="l", lwd=1, col="gray", las = 1, xaxt = 'n',mar=c(0,0,0,0))
	for (h in (1:dim(e2)[2])) {
		lines(e2[,h], ylim=c(0,1), xlab = "", ylab="log2FC_nor", type="l", lwd=1, col="gray", xaxt = 'n',mar=c(0,0,0,0))
	}
	
	##get median
	m=c()
	for (j in 1:dim(e2)[1]){
		med = median(e2[j,(1:dim(e2)[2])])
		m = c(m, med)
	}
	lines(m, ylim=c(0,1), xlab = "", ylab="log2FC",type="l", lwd=1, col="red", xaxt = 'n',mar=c(0,0,0,0))

	##get quartile
	q=c()
	q2=c()
	for (k in 1:dim(e2)[1]){
		quant= quantile(e2[k,(1:dim(e2)[2])])
		up= quant[2]
		dwn= quant[4]
		q = c(q,up)
		q2 = c(q2, dwn)
	}
	lines(q, ylim=c(0,1), xlab = "", ylab="log2FC_nor", type="l", lwd=1, col="blue", xaxt = 'n',mar=c(0,0,0,0))
	lines(q2, ylim=c(0,1), xlab = "", ylab="log2FC_nor", type="l", lwd=1, col="blue", xaxt = 'n',mar=c(0,0,0,0))
	x_labe1 <- c(row.names(e2))
	##change based on the number of columns-1
	axis(1,at=1:92, labels=x_labe1, col.axis="blue", las=0, cex.lab=1,cex.axis=1)}

dev.off() #turn off pdf

