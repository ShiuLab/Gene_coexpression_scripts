args = commandArgs(TRUE)

inputFile = args[1]
outputPrefix = args[2]
selectedPower = as.numeric(args[3])
cutHeight = as.numeric(args[4])
sign = args[5]
path = args[6]

setwd(path)

library(WGCNA);
options(stringsAsFactors = FALSE);

# To read the data back into memory.
inputFileName = paste(inputFile, sep="")
lnames = load(file = inputFileName);
#The variable lnames contains the names of loaded variables.
lnames

selectedPower
cutHeight

# The Scale Independence graph and the Mean connectivity graph show that a thresholding power of 9 or 10 are the values closest to an R^2 value of 0.9.

# maxBlockSize is the maximum number of probes/genes that will be analyzed at once.  It is best that this value is greater than the total number of probes, but the true best value is dependent on the memory in the machine.  Supposedly, 20000 probes will require about 16 Gb.  
# So far, only 4.5 Gb of RAM is being used.  So, maybe the block size can be set to 40000 without fear of using too much RAM.

tomFile = paste(outputPrefix, "_TOM", sep="")

net = blockwiseModules(riceData, power = selectedPower, minModuleSize = 30,
maxBlockSize = 25000, networkType = sign,
reassignThreshold = 0, mergeCutHeight = 0.25,
detectCutHeight = cutHeight,
numericLabels = TRUE, pamRespectsDendro = FALSE,
saveTOMs = TRUE,
saveTOMFileBase = tomFile,
verbose = 3)

table(net$colors)

# 30 modules created (the max, by default).
# open a graphics window or save to a pdf
#sizeGrWindow(12, 9)
#networkTreeFile = paste("network_tree_w_modules_", experimentName, "_cov_", cov, "_beta_", selectedPower, "_cut_", cutHeight, ".pdf", sep="")
#pdf(file=networkTreeFile, width=9, height=5)

# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)

# Plot the dendrogram and the module colors underneath
#plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
#"Module colors",
#dendroLabels = FALSE, hang = 0.03,
#addGuide = TRUE, guideHang = 0.05)
#dev.off()

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];

saveFileName = paste(outputPrefix, "-networkConstruction-auto.RData", sep="")
save(riceData, gsg, sampleTree, powers, nGenes, nSamples, cex1, sft, net, geneTree, 
     mergedColors, MEs, moduleColors, moduleLabels, selectedPower,
     file = saveFileName)
#save(sampleTree, nGenes, nSamples, geneTree, gsg, 
#mergedColors, MEs, moduleColors, moduleLabels, net, powers, riceData, 
#sampleTree, file = "rice_GSE6893_flower_seed_x-networkConstruction-auto_no_cex_sft.RData")


# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
saveDissTOM = paste(outputPrefix, "-all_data_up_to_dissTOM-auto.RData", sep="")
dissTOM = 1-TOMsimilarityFromExpr(riceData, power = selectedPower);
save(riceData, gsg, sampleTree, powers, nGenes, nSamples, cex1, sft, net, geneTree, 
     mergedColors, MEs, moduleColors, moduleLabels, dissTOM, selectedPower,
     file = saveDissTOM)

# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
#plotTOM = dissTOM^7;
# Set diagonal to NA for a nicer plot
# diag(plotTOM) = NA;  # This function caused trouble in TOMplot().

#nSelect = 400
# For reproducibility, we set the random seed
#set.seed(10);
#select = sample(nGenes, size = nSelect);
#selectTOM = dissTOM[select, select];
# There is no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
#selectTree = flashClust(as.dist(selectTOM), method = "average")
#selectColors = moduleColors[select];
# Open a graphical window
#sizeGrWindow(9,9)
#smallHeatMap = paste("all_by_all_network_heatmap_400_random_genes_", experimentName, "_cov_", cov, "_beta_", selectedPower, "_cut_", cutHeight, ".pdf", sep="")
#title = paste("Network heatmap plot ", experimentName, " cov ", cov, " beta ", selectedPower, " cut ", cutHeight, ", selected genes", sep="")
#pdf(file=smallHeatMap, width=9, height=9)
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing
# the color palette; setting the diagonal to NA also improves the clarity of the plot
#plotDiss = selectTOM^7;
#diag(plotDiss) = NA;  # This function caused trouble in TOMplot().
#TOMplot(plotDiss, selectTree, selectColors, main = title)
#dev.off()


nGenes = ncol(riceData)
nSamples = nrow(riceData)
nSelect = 1000
if (nGenes < 1000) {
  nSelect = nGenes
}
# For reproducibility, we set the random seed
set.seed(10);
select = sample(nGenes, size = nSelect);
selectTOM = dissTOM[select, select];
# There isno simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree = flashClust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select];
# Open a graphical window
#sizeGrWindow(9,9)
largeHeatMap = paste(outputPrefix, "all_by_all_network_heatmap_", nSelect, "_random_genes", ".pdf", sep="")
title = paste("Network heatmap plot ", outputPrefix, sep="")
pdf(file=largeHeatMap, width=9, height=9)
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing
# the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDiss = selectTOM^7;
diag(plotDiss) = NA;  # This function caused trouble in TOMplot().
TOMplot(plotDiss, selectTree, selectColors, main = title)
dev.off()

