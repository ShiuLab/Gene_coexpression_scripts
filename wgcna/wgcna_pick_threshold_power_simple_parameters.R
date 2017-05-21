args = commandArgs(TRUE)

inputFileName = args[1]
outputPrefix = args[2]
path = args[3]

setwd(path)

library(WGCNA);
options(stringsAsFactors = FALSE);

inputFile = paste(inputFileName, sep="")
lnames = load(file = inputFile);
#The variable lnames contains the names of loaded variables.
lnames

nGenes = ncol(riceData)
nSamples = nrow(riceData)


# Choose a soft-thresholding power value.
powers = c(c(1:60), seq(from = 12, to=70, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(riceData, powerVector = powers, verbose = 5)
# Plot the results:
#sizeGrWindow(9, 5)
pdfName = paste(outputPrefix, "_pick_soft_threshold_", ".pdf", sep="")
pdf(file=pdfName, width=9, height=5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
main = paste("Scale independence COV ", outputPrefix, sep=""));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

# Save our progress.
save(riceData, gsg, sampleTree, powers, nGenes, nSamples, cex1, sft,
     file = paste(outputPrefix, "-threshold_data.RData", sep=""))

q(save="no")
