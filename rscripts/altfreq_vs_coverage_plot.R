if (!require("ggplot2",character.only = TRUE))
{
      install.packages("ggplot2",dep=TRUE,repos='http://probability.ca/cran/')
}

if (!require("grid",character.only = TRUE))
{
      install.packages("grid",dep=TRUE,repos='http://probability.ca/cran/')
}


if (!require("scales",character.only = TRUE))
{
      install.packages("scales",dep=TRUE,repos='http://probability.ca/cran/')
}

library("ggplot2")
library("grid")
library("scales")

# Remove objects from previous sessions
rm(list = ls())

args=(commandArgs(TRUE))

sampleID = args[1]
sampleAnalysisFile = args[2]
dirPath = args[3]
outFolder = args[4]

print(args)

#Read in coverage data from file

fileName = paste(outFolder, "/", sampleID, "_altfreq_vs_coverage.png", sep = "")

#Plot and save to PDF
png(filename=fileName, width=864, height=605, pointsize=12)

	files <- list.files(path=dirPath, pattern="*_filter_av.txt", full.names=T, recursive=TRUE)
	
	sampleTable <- read.table(sampleAnalysisFile, header=TRUE, sep=",", row.names=NULL)
	
	for (file in files) {

		fileName = basename(file)
		
		# Need to use grepl because covFileName is the raw sampleID, sampleID contains the  "_S##" suffix added by MiSeq
		covFileName = sub("^(.+)_filter_av.txt", "\\1", fileName) 
		print (fileName)
		print (covFileName)
		
		# If covFileName is in the samplesheet, then process
		# Later, can pass in analysis type and do the same thing so only coverage of samples from same type is compared
		if(any(sampleTable[c("sample")] == covFileName))
		{
			newCovTable <- read.table(file, header=TRUE, sep=",", row.names=NULL, stringsAsFactors=FALSE)	

			# Select only DP and AD.ALT
			covTable <- subset(newCovTable, select=c(DP, AD.ALT))

			# Replace values "-" with NA
			covTable$AD.ALT[covTable$AD.ALT == "-"] <- NA
			covTable$DP[covTable$DP == "-"] <- NA

			# Drops row with NA
			covTable <- covTable[complete.cases(covTable),]
			
			if (nrow(covTable) > 0)
			{
				covTable$DP <- as.numeric(covTable$DP)
				covTable$AD.FREQ <- as.numeric(covTable$AD.ALT) / covTable$DP
			}
			else
			{
				covTable <- data.frame(DP=c(0), AD.FREQ=c(0))
			}

			print(sampleID)
			print(covFileName)
			
			if(grepl(covFileName, sampleID)) {
					mainTable <- covTable
			}
			else {
				if(exists("myPlot"))
				{
					myPlot <- myPlot + 
					geom_point(data = covTable, aes(x = DP, y = AD.FREQ, colour = "Other Samples"), alpha = 0.5, size = 3)
				}
				else
				{
					myPlot <- ggplot(covTable, aes(x = DP, y = AD.FREQ, colour="Other Samples")) + 
					geom_point(alpha = 0.5, size = 3) 
				
				}	
			}
		}

	}

	# Plot
	
	myPlot <- myPlot + geom_point(data = mainTable, aes(x = DP, y = AD.FREQ, colour = sampleID), alpha = 1.0, size = 3) +
	scale_colour_manual("", breaks = c(sampleID, "Other Samples"), values = c("red", "black")) +
	theme_classic() +
	labs(list(y = "Allelic Fraction", x = "Depth of Coverage")) +
	labs(title = paste("Allelic Fraction of ", sampleID, " against other samples in run", sep = "")) +
	labs(colour = "Legend") + 
	expand_limits(y = 0) +
	expand_limits(y = 1.01) +
	expand_limits(x = 800) +
	scale_y_continuous(expand = c(0, 0)) +
	scale_x_continuous(expand = c(0, 0)) +
	theme(axis.text = element_text(size = 12)) + 
	theme(axis.title.x = element_text(size = 16, vjust = -0.5)) + 
	theme(axis.title.y = element_text(size = 16, vjust = -0.05)) + 
	theme(plot.title = element_text(size = 22, vjust = 2)) +
	theme(legend.text = element_text(size=14)) + 
	theme(legend.key.size = unit(0.8, "cm")) + 
	theme(plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm"))
	
	myPlot

dev.off()

warnings()
