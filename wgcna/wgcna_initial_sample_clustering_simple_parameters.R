args = commandArgs(TRUE)

inputFileName = args[1]
outputPrefix = args[2]
thisPath = args[3]

path = thisPath
setwd(path)

library(WGCNA);
options(stringsAsFactors = FALSE);

inputFile = paste(inputFileName, sep="")
riceData0 = read.table(inputFile, header = TRUE, sep = "\t")

riceData = as.data.frame(t(riceData0[,-c(1)]));

#riceData

names(riceData) = t(riceData0[,1])

rownames(riceData) = names(riceData0)[-c(1)];

# Genes with all 0.0 expression values or too many missing or 0.0 values will be discarded.
gsg = goodSamplesGenes(riceData, verbose = 3);

if (!gsg$allOK)
{
# Optionally, print the gene and sample names that were removed:
if (sum(!gsg$goodGenes)>0)
printFlush(paste("Removing genes:", paste(names(riceData)[!gsg$goodGenes], collapse = ", ")));
if (sum(!gsg$goodSamples)>0)
printFlush(paste("Removing samples:", paste(rownames(riceData)[!gsg$goodSamples], collapse = ", ")));
# Remove the offending genes and samples from the data:
riceData = riceData[gsg$goodSamples, gsg$goodGenes]
}

sampleTree = flashClust(dist(riceData), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
sizeGrWindow(10,9)
# pdf() paired with dev.off will write the results to a pdf file.
# Use either sizeGrWindow() or pdf().  Both will not work at the same time.
chipTreeOutput = paste(outputPrefix, "_chip_tree.pdf", sep="")
pdf(file=chipTreeOutput, width=10, height=9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plotTitle = paste(outputPrefix, " Sample clustering to detect outliers - average", sep="")
plot(sampleTree, main = plotTitle, sub="", xlab="", cex.lab = 1.5,
cex.axis = 1.5, cex.main = 2)
dev.off()

outputFile = paste(outputPrefix, "-dataInput.RData", sep="")
save(riceData, gsg, sampleTree, file = outputFile)
q(save="no")
