#module load R/3.2.0
#install.packages("corpcor", repos="http://cran.r-project.org")
library(corpcor)
cat("
args:
inp1 = expression dataset
")
args = commandArgs(TRUE)
files = args[1]
dat <- read.table(files,head=T,sep='\t',row.names=1,stringsAsFactors=F)
dat2 <- t(dat)
co1 <- pcor.shrink(dat2,lambda=0, verbose=TRUE)
co2 <- round(co1,3)
write.table(co2,paste('/mnt/scratch/peipeiw/Pathway_prediction/corpcor/corpcor_',substr(files,73,nchar(files)),sep=''),quote=F,sep='\t')




