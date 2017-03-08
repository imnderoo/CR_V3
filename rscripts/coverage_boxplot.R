args=(commandArgs(TRUE))
print(args)

covBase = basename(args[1])
outFolder = args[2]
depthThreshold = (args[3])

covFileName = sub("^([^.]*).*", "\\1", covBase) 
#Read in coverage data from file
covData = read.table(args[1], header=T, sep="\t", row.names=NULL)
geneName = covData[,1]
exonName = covData[,3]
startPos = covData[,5]
endPos = covData[,6]
q1=covData[,12]
q2=covData[,13]
q3=covData[,14]
covStats=covData[,11:15]
yMax=max(covData[,15])
yMax=1.1*yMax

#The numver of Observations (numObs) can ideally estimated by total_coverage / average_coverage --> yielding number of bases (observations)
#However, the exact number is not needed since it is only used for confidence interval to display outliers (which we don't)
numObs=endPos-startPos

fileName = paste(outFolder, "/", covFileName, ".png", sep = "")

#conf is calculated as: median +/- 1.58 * (Q3-Q1) / sqrt(n)
confMin = q2 - 1.58 * (q3 - q1) / sqrt(numObs)
confMax = q2 + 1.58 * (q3 - q1) / sqrt(numObs)

#This list is the format required by the bxp command
confSum = list(stats=t(covStats), n=numObs, conf=t(matrix(c(confMin,confMax), nrow=nrow(covData), ncol=2)), out=numeric(0), group=numeric(0), names=exonName)

pdf.options()

#sampleName = "sample"
#geneName = "gene"

a=strsplit(covFileName, "-")
	
# If we don't want the _S2, then change to sampleName=a[[1]][1]
#sampleName=paste(a[[1]][1:2], collapse="_")

# Need to use grepl because covFileName is the raw sampleID, sampleID contains the  "_S##" suffix added by MiSeq		
sampleName = sub("^(.+)-aln-pe-sorted-merged-dedup-recal-covQ20.*", "\\1", covFileName) 
		
# Need this fix for MiSeq cases
if (grepl("_S", covFileName))
{
	sampleName = sub("^(.+)_S.*-aln-pe-sorted-merged-dedup-recal-covQ20.*", "\\1", covFileName) 
}

#geneName=a[[1]][length(a[[1]])]		

title=paste("Distribution of Coverage for",geneName[1],"in Sample ",sampleName,sep=" ")

#Plot and save to PDF
png(filename=fileName, width=864, height=605, pointsize=12)
bxp(confSum, main=title, ylim=c(0,yMax), yaxs="i", ylab="Depth of Coverage", las=2, cex.lab=1, cex.axis=0.7) +
abline(h=depthThreshold, lty=2, col="red")
dev.off()
