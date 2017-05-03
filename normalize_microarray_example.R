source("http://bioconductor.org/biocLite.R")
biocLite("affy")
library(AnnotationDbi)
library(limma)
library(affy)

M <- ReadAffy() # Get the CEL files in the directory (all cel files that you want to normalize and get differential expression should be in the working directory)
M.rma   <- rma(M, background = T, normalize = T) # Background correction and normalization.
M.rma.e <- exprs(M.rma)
list.celfiles()
dim(M.rma.e)    # Check the dimension of the dataset
write.table(M.rma.e , file="normalized_intensity_0324", quote=F, sep='\t') #output normalized intensity

##QC arrays
source("http://bioconductor.org/biocLite.R")
biocLite("simpleaffy")
library(simpleaffy)  	# Load the simpleaffy package
M.qc <- qc(M)			# Call the qc function
plot(M.qc)			# Plot the qc results

## check for degradation
RNAdeg <- AffyRNAdeg(M)
plotAffyRNAdeg(RNAdeg, col=c(1,2,3,4))
summaryAffyRNAdeg(RNAdeg)

#for differential expression :
design7641 <- cbind(
  epi_lat.c          = c(1,1,1, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0),
  epi_lat.s          = c(0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 1,1,1, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0),
  col.c              = c(0,0,0, 1,1,1, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0),
  col.s              = c(0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 1,1,1, 0,0,0, 0,0,0, 0,0,0, 0,0,0),
  cortex.c           = c(0,0,0, 0,0,0, 1,1,1, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0),
  cortex.s           = c(0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 1,1,1, 0,0,0, 0,0,0, 0,0,0),
  endo_qui.c         = c(0,0,0, 0,0,0, 0,0,0, 1,1,1, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0),
  endo_qui.s         = c(0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 1,1,1, 0,0,0, 0,0,0),
  stele.c            = c(0,0,0, 0,0,0, 0,0,0, 0,0,0, 1,1,1, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0),
  stele.s            = c(0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 1,1,1, 0,0,0),
  protophl.c         = c(0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 1,1,1, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0),
  protophl.s         = c(0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 1,1,1)) #design is designating which CEL files are control and which ones are treatment

#design hormone

design <-cbind(
  mock30 = c(1,1, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 
             0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 
             0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0,
             0,0, 0,0, 0,0, 0,0, 0,0, 0,0),
  IAA30 = c(0,0, 1,1, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 
             0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 
             0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0,
             0,0, 0,0, 0,0, 0,0, 0,0, 0,0),
  zeatin30 = c(0,0, 0,0, 1,1, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 
            0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 
            0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0,
            0,0, 0,0, 0,0, 0,0, 0,0, 0,0),
  GA30 = c(0,0, 0,0, 0,0, 1,1, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 
            0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 
            0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0,
            0,0, 0,0, 0,0, 0,0, 0,0, 0,0),
  ABA30 = c(0,0, 0,0, 0,0, 0,0, 1,1, 0,0, 0,0, 0,0, 0,0, 0,0, 
            0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 
            0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0,
            0,0, 0,0, 0,0, 0,0, 0,0, 0,0),
  MJ30 = c(0,0, 0,0, 0,0, 0,0, 0,0, 1,1, 0,0, 0,0, 0,0, 0,0, 
            0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 
            0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0,
            0,0, 0,0, 0,0, 0,0, 0,0, 0,0),
  ACC30 = c(0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 1,1, 0,0, 0,0, 0,0, 
            0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 
            0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0,
            0,0, 0,0, 0,0, 0,0, 0,0, 0,0),
  BL30 = c(0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 1,1, 0,0, 0,0, 
            0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 
            0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0,
            0,0, 0,0, 0,0, 0,0, 0,0, 0,0),
  mock1 = c(0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 1,1, 0,0, 
            0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 
            0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0,
            0,0, 0,0, 0,0, 0,0, 0,0, 0,0),
  IAA1 = c(0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 1,1, 
            0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 
            0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0,
            0,0, 0,0, 0,0, 0,0, 0,0, 0,0),
  zeatin1 = c(0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 
            1,1, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 
            0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0,
            0,0, 0,0, 0,0, 0,0, 0,0, 0,0),
  GA1 = c(0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 
            0,0, 1,1, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 
            0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0,
            0,0, 0,0, 0,0, 0,0, 0,0, 0,0),
  ABA1 = c(0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 
            0,0, 0,0, 1,1, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 
            0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0,
            0,0, 0,0, 0,0, 0,0, 0,0, 0,0),
  MJ1 = c(0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 
            0,0, 0,0, 0,0, 1,1, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 
            0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0,
            0,0, 0,0, 0,0, 0,0, 0,0, 0,0),
  ACC11 = c(0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 
            0,0, 0,0, 0,0, 0,0, 1,1, 0,0, 0,0, 0,0, 0,0, 0,0, 
            0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0,
            0,0, 0,0, 0,0, 0,0, 0,0, 0,0),
  BL1 = c(0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 
            0,0, 0,0, 0,0, 0,0, 0,0, 1,1, 0,0, 0,0, 0,0, 0,0, 
            0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0,
            0,0, 0,0, 0,0, 0,0, 0,0, 0,0),
  mock3 = c(0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 
            0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 1,1, 0,0, 0,0, 0,0, 
            0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0,
            0,0, 0,0, 0,0, 0,0, 0,0, 0,0),
  IAA3 = c(0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 
            0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 1,1, 0,0, 0,0, 
            0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0,
            0,0, 0,0, 0,0, 0,0, 0,0, 0,0),
  zeatin3 = c(0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 
            0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 1,1, 0,0, 
            0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0,
            0,0, 0,0, 0,0, 0,0, 0,0, 0,0),
  GA3 = c(0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 
            0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 1,1, 
            0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0,
            0,0, 0,0, 0,0, 0,0, 0,0, 0,0),
  ABA3 = c(0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 
            0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 
            1,1, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0,
            0,0, 0,0, 0,0, 0,0, 0,0, 0,0),
  MJ3 = c(0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 
            0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 
            0,0, 1,1, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0,
            0,0, 0,0, 0,0, 0,0, 0,0, 0,0),
  ACC3 = c(0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 
            0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 
            0,0, 0,0, 1,1, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0,
            0,0, 0,0, 0,0, 0,0, 0,0, 0,0),
  BL3 = c(0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 
            0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 
            0,0, 0,0, 0,0, 1,1, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0,
            0,0, 0,0, 0,0, 0,0, 0,0, 0,0),
  GA15m30 = c(0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 
               0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 
               0,0, 0,0, 0,0, 0,0, 1,1, 0,0, 0,0, 0,0, 0,0, 0,0,
               0,0, 0,0, 0,0, 0,0, 0,0, 0,0),
  GA15G30 = c(0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 
               0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 
               0,0, 0,0, 0,0, 0,0, 0,0, 1,1, 0,0, 0,0, 0,0, 0,0,
               0,0, 0,0, 0,0, 0,0, 0,0, 0,0),
  GA15m1 = c(0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 
               0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 
               0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 1,1, 0,0, 0,0, 0,0,
               0,0, 0,0, 0,0, 0,0, 0,0, 0,0),
  GA15G1 = c(0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 
               0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 
               0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 1,1, 0,0, 0,0,
               0,0, 0,0, 0,0, 0,0, 0,0, 0,0),
  GA15m3 = c(0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 
               0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 
               0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 1,1, 0,0,
               0,0, 0,0, 0,0, 0,0, 0,0, 0,0),
  GA15G3 = c(0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 
               0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 
               0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 1,1,
               0,0, 0,0, 0,0, 0,0, 0,0, 0,0),
  det2m30 = c(0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 
               0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 
               0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0,
               1,1, 0,0, 0,0, 0,0, 0,0, 0,0),
  det2BL30 = c(0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 
               0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 
               0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0,
               0,0, 1,1, 0,0, 0,0, 0,0, 0,0),
  det2m1 = c(0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 
               0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 
               0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0,
               0,0, 0,0, 1,1, 0,0, 0,0, 0,0),
  det2BL1 = c(0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 
               0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 
               0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0,
               0,0, 0,0, 0,0, 1,1, 0,0, 0,0),
  det2m3 = c(0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 
               0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 
               0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0,
               0,0, 0,0, 0,0, 0,0, 1,1, 0,0),
  det2BL3 = c(0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 
               0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 
               0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0, 0,0,
               0,0, 0,0, 0,0, 0,0, 0,0, 1,1))
  
design_aphid = cbind(cont = c(1, 1, 1, 0, 0), treat = c(0, 0, 0, 1, 1))
contrast = makeContrasts(treat-cont, levels = design_aphid)
  
contrast <- makeContrasts(IAA30-mock30, zeatin30-mock30, GA30-mock30, ABA30-mock30, MJ30-mock30, ACC30-mock30, BL30-mock30, 
                          IAA1-mock1, zeatin1-mock1, GA1-mock1, ABA1-mock1, MJ1-mock1, ACC11-mock1, BL1-mock1,
                         IAA3-mock3, zeatin3-mock3, GA3-mock3, ABA3-mock3, MJ3-mock3, ACC3-mock3, BL3-mock3,
                         GA15m30-mock30, det2m30-mock30, GA15G30- mock30, GA15m3-mock3, det2m3-mock3, GA15G3-mock3,
                         GA15m1-mock1, det2m1-mock1, GA15G1-mock1, det2BL30-mock30, det2BL3-mock3, det2BL1-mock1,
                         levels = design)
  
#contrast <- makeContrasts(epi_lat.s-epi_lat.c, col.s-col.c, cortex.s-cortex.c, levels = design) 
#make contrast between treatment(salt in this case) - control according to your design #only 3 celltypes as examples
              
design.df <- as.data.frame(design_aphid)
fit <- lmFit(M.rma.e, design_aphid) #fit data to the design
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)
n <- length(featureNames(M.rma))

#n
#top.epi<- topTable(fit2, coef ="epi_lat.s-epi_lat.c", n=n, adjust = 'fdr')
#top.col<- topTable(fit2, coef ="col.s-col.c", n=n, adjust = 'fdr')
top_aphid <- topTable(fit2, coef ="treat - cont", n=n, adjust = 'fdr')

top_IAA30 <- topTable(fit2, coef ="IAA30 - mock30", n=n, adjust = 'fdr')
top_zeatin30<- topTable(fit2, coef ="zeatin30 - mock30", n=n, adjust = 'fdr')
top_GA30<- topTable(fit2, coef ="GA30 - mock30", n=n, adjust = 'fdr')
top_ABA30<- topTable(fit2, coef ="ABA30 - mock30", n=n, adjust = 'fdr')
top_MJ30<- topTable(fit2, coef ="MJ30 - mock30", n=n, adjust = 'fdr')
top_ACC30<- topTable(fit2, coef ="ACC30 - mock30", n=n, adjust = 'fdr')
top_BL30<- topTable(fit2, coef ="BL30 - mock30", n=n, adjust = 'fdr')
top_GA15m30<- topTable(fit2, coef ="GA15G30 - mock30", n=n, adjust = 'fdr')
top_det2m30<- topTable(fit2, coef ="det2m30 - mock30", n=n, adjust = 'fdr')
top_GA15G30<- topTable(fit2, coef ="GA15G30 - mock30", n=n, adjust = 'fdr')
top_det2BL30<- topTable(fit2, coef ="det2BL30 - mock30", n=n, adjust = 'fdr')
top_IAA1<- topTable(fit2, coef ="IAA1 - mock1", n=n, adjust = 'fdr')
top_zeatin1<- topTable(fit2, coef ="zeatin1 - mock1", n=n, adjust = 'fdr')
top_GA1<- topTable(fit2, coef ="GA1 - mock1", n=n, adjust = 'fdr')
top_ABA1<- topTable(fit2, coef ="ABA1 - mock1", n=n, adjust = 'fdr')
top_MJ1<- topTable(fit2, coef ="MJ1 - mock1", n=n, adjust = 'fdr')
top_ACC1<- topTable(fit2, coef ="ACC11 - mock1", n=n, adjust = 'fdr')
top_BL1<- topTable(fit2, coef ="BL1 - mock1", n=n, adjust = 'fdr')
top_GA15m1<- topTable(fit2, coef ="GA15G1 - mock1", n=n, adjust = 'fdr')
top_det2m1<- topTable(fit2, coef ="det2m1 - mock1", n=n, adjust = 'fdr')
top_GA15G1<- topTable(fit2, coef ="GA15G1 - mock1", n=n, adjust = 'fdr')
top_det2BL1<- topTable(fit2, coef ="det2BL1 - mock1", n=n, adjust = 'fdr')
top_IAA3<- topTable(fit2, coef ="IAA3 - mock3", n=n, adjust = 'fdr')
top_zeatin3<- topTable(fit2, coef ="zeatin3 - mock3", n=n, adjust = 'fdr')
top_GA3<- topTable(fit2, coef ="GA3 - mock3", n=n, adjust = 'fdr')
top_ABA3<- topTable(fit2, coef ="ABA3 - mock3", n=n, adjust = 'fdr')
top_MJ3<- topTable(fit2, coef ="MJ3 - mock3", n=n, adjust = 'fdr')
top_ACC3<- topTable(fit2, coef ="ACC3 - mock3", n=n, adjust = 'fdr')
top_BL3<- topTable(fit2, coef ="BL3 - mock3", n=n, adjust = 'fdr')
top_GA15m3<- topTable(fit2, coef ="GA15G3 - mock3", n=n, adjust = 'fdr')
top_det2m3<- topTable(fit2, coef ="det2m3 - mock3", n=n, adjust = 'fdr')
top_GA15G3<- topTable(fit2, coef ="GA15G3 - mock3", n=n, adjust = 'fdr')
top_det2BL3<- topTable(fit2, coef ="det2BL3 - mock3", n=n, adjust = 'fdr')



write.table(top_IAA30, file = "contrast_IAA30", quote = F, sep = '\t', row.names = T)
write.table(top_zeatin30, file = "contrast_zeatin30", quote = F, sep = '\t', row.names = T)
write.table(top_GA30, file = "contrast_GA30", quote = F, sep = '\t', row.names = T)
write.table(top_ABA30, file = "contrast_ABA30", quote = F, sep = '\t', row.names = T)
write.table(top_MJ30, file = "contrast_MJ30", quote = F, sep = '\t', row.names = T)
write.table(top_ACC30, file = "contrast_ACC30", quote = F, sep = '\t', row.names = T)
write.table(top_BL30, file = "contrast_BL30", quote = F, sep = '\t', row.names = T)
write.table(top_GA15m30, file = "contrast_GA15m30", quote = F, sep = '\t', row.names = T)
write.table(top_det2m30, file = "contrast_det2m30", quote = F, sep = '\t', row.names = T)
write.table(top_GA15G30, file = "contrast_GA15G30", quote = F, sep = '\t', row.names = T)
write.table(top_det2BL30, file = "contrast_det2BL30", quote = F, sep = '\t', row.names = T)
write.table(top_IAA1, file = "contrast_IAA1", quote = F, sep = '\t', row.names = T)
write.table(top_zeatin1, file = "contrast_zeatin1", quote = F, sep = '\t', row.names = T)
write.table(top_GA1, file = "contrast_GA1", quote = F, sep = '\t', row.names = T)
write.table(top_ABA1, file = "contrast_ABA1", quote = F, sep = '\t', row.names = T)
write.table(top_MJ1, file = "contrast_MJ1", quote = F, sep = '\t', row.names = T)
write.table(top_ACC1, file = "contrast_ACC1", quote = F, sep = '\t', row.names = T)
write.table(top_BL1, file = "contrast_BL1", quote = F, sep = '\t', row.names = T)
write.table(top_GA15m1, file = "contrast_GA15m1", quote = F, sep = '\t', row.names = T)
write.table(top_det2m1, file = "contrast_det2m1", quote = F, sep = '\t', row.names = T)
write.table(top_GA15G1, file = "contrast_GA15G1", quote = F, sep = '\t', row.names = T)
write.table(top_det2BL1, file = "contrast_det2BL1", quote = F, sep = '\t', row.names = T)
write.table(top_IAA3, file = "contrast_IAA3", quote = F, sep = '\t', row.names = T)
write.table(top_zeatin3, file = "contrast_zeatin3", quote = F, sep = '\t', row.names = T)
write.table(top_GA3, file = "contrast_GA3", quote = F, sep = '\t', row.names = T)
write.table(top_ABA3, file = "contrast_ABA3", quote = F, sep = '\t', row.names = T)
write.table(top_MJ3, file = "contrast_MJ3", quote = F, sep = '\t', row.names = T)
write.table(top_ACC3, file = "contrast_ACC3", quote = F, sep = '\t', row.names = T)
write.table(top_BL3, file = "contrast_BL3", quote = F, sep = '\t', row.names = T)
write.table(top_GA15m3, file = "contrast_GA15m3", quote = F, sep = '\t', row.names = T)
write.table(top_det2m3, file = "contrast_det2m3", quote = F, sep = '\t', row.names = T)
write.table(top_GA15G3, file = "contrast_GA15G3", quote = F, sep = '\t', row.names = T)
write.table(top_det2BL3, file = "contrast_det2BL3", quote = F, sep = '\t', row.names = T)


write.table(top_aphid, file = "contrast_aphid", quote = F, sep = '\t', row.names = T)
write.table(top.col, file = "contrast_col", quote = F, sep = '\t', row.names = T)