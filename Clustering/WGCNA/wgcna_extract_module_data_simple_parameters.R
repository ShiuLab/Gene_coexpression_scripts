args = commandArgs(TRUE)

inputFile = args[1]
newDirectoryName = args[2]
path = args[3]

setwd(path)

library(WGCNA);
options(stringsAsFactors = FALSE);
inputFileName = paste(inputFile, sep="")
lnames = load(file = inputFileName)
unique(moduleColors)

# Select module probes
probes = names(riceData)

newDir = paste(path, newDirectoryName, sep="")
dir.create(newDir)
setwd(newDir)

print_module_colors <- function(moduleColors) {
    colors = unique(moduleColors)
    for (i in 1:length(colors)) {
        inModule = (moduleColors==colors[i]);
	modProbes = probes[inModule];
	output_name = paste(colors[i], "_module.txt", sep = "");
	write(modProbes, file = output_name, ncolumns = 1, append = FALSE)
    }
}

# Write a list of tu ids for each module.
print_module_colors(moduleColors)

nGenes
