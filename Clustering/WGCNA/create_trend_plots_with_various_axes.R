args = commandArgs(TRUE)

path = args[1]
moduleColor = args[2]
numberTreatments = as.integer(args[3])
graphSize = as.integer(args[4])
lineColor = args[5]
experimentName = args[6]
zLimit = as.integer(args[7])
negZLimit = as.integer(args[8])
treatmentLabels = strsplit(args[9],";")

setwd(path)

inputFile =  paste(moduleColor, "_module_expression_matrix_z_scores.txt", sep="")

expression0 = read.table(inputFile, header = TRUE, sep = "	")
expression = (t(expression0[,-c(1)]));

output1 = paste(experimentName, "_", moduleColor, "_module_A.pdf", sep="")
output2 = paste(experimentName, "_", moduleColor, "_module_B.pdf", sep="")
output3 = paste(experimentName, "_", moduleColor, "_module_C.pdf", sep="")
output4 = paste(experimentName, "_", moduleColor, "_module_D.pdf", sep="")

# output1 gets both x and y axis labels.
pdf(output1, height=5, width=graphSize, paper="special", pagecentre=TRUE) 
par(mar=c(12,8,4,2)+0.1)
split.screen(c(1,1));
plot(expression[,1], axes=F, ylim=c(negZLimit,zLimit), xlab="", ylab="Z-score", type="l", lwd=0.5, cex.lab=1, col=lineColor);
for(i in 2:length(expression[1,])) {screen(1, new=FALSE); lines(expression[,i], lwd=0.5, col=lineColor); };
axis(1, at=1:numberTreatments, labels=FALSE, tick=TRUE)
text(1:numberTreatments, par("usr")[3] - 0.75, srt=45, adj=1, labels=treatmentLabels[[1]], xpd=TRUE, cex=1)
axis(2, las=1, at=1*zLimit:negZLimit, cex.axis=1)
box()
dev.off()

# output2 only gets x axis labels.
pdf(output2, height=5, width=graphSize, paper="special", pagecentre=TRUE) 
par(mar=c(12,8,4,2)+0.1)
split.screen(c(1,1));
plot(expression[,1], axes=F, ylim=c(negZLimit,zLimit), xlab="", ylab="", type="l", lwd=0.5, col=lineColor);
for(i in 2:length(expression[1,])) {screen(1, new=FALSE); lines(expression[,i], lwd=0.5, col=lineColor); };
axis(1, at=1:numberTreatments, labels=FALSE, tick=TRUE)
text(1:numberTreatments, par("usr")[3] - 0.75, srt=45, adj=1, labels=treatmentLabels[[1]], xpd=TRUE, cex=1)
box()
dev.off() 

# output3 only gets y axis labels.
pdf(output3, height=5, width=graphSize, paper="special", pagecentre=TRUE) 
par(mar=c(12,8,4,2)+0.1)
split.screen(c(1,1));
plot(expression[,1], axes=F, ylim=c(negZLimit,zLimit), xlab="", ylab="Z-score", type="l", lwd=0.5, cex.lab=1, col=lineColor);
for(i in 2:length(expression[1,])) {screen(1, new=FALSE); lines(expression[,i], lwd=0.5, col=lineColor); };
axis(2, las=1, at=1*zLimit:negZLimit, cex.axis=1)
box()
dev.off() 

# output4 gets no axis labels
pdf(output4, height=5, width=graphSize, paper="special", pagecentre=TRUE) 
par(mar=c(12,8,4,2)+0.1)
split.screen(c(1,1));
plot(expression[,1], axes=F, ylim=c(negZLimit,zLimit), xlab="", ylab="", type="l", lwd=0.5, col=lineColor);
for(i in 2:length(expression[1,])) {screen(1, new=FALSE); lines(expression[,i], lwd=0.5, col=lineColor); };
box()
dev.off()
