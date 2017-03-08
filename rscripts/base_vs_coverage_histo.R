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

rm(list = ls())

args=(commandArgs(TRUE))

sampleID = args[1]
sampleAnalysisFile = args[2]
dirPath = args[3]
outFolder = args[4]

print(args)

#Read in coverage data from file

fileName = paste(outFolder, "/", sampleID, "_bases_vs_coverage.png", sep = "")

#Plot and save to PDF
png(filename=fileName, width=864, height=605, pointsize=12)

	files <- list.files(path=dirPath, pattern="*all_only.hist", full.names=T, recursive=FALSE)
	sampleTable <- read.table(sampleAnalysisFile, header=TRUE, sep=",", row.names=NULL)

	maxY = 0

	for (file in files) 
	{
		# Coverage File Name - try to parse out the sampleID name
		
		fileName = basename(file)
		
		# Need to use grepl because covFileName is the raw sampleID, sampleID contains the  "_S##" suffix added by MiSeq		
		covFileName = sub("^(.+)-aln-pe-sorted-merged-dedup-recal-covQ20_all_only.hist", "\\1", fileName) 
		
		# Need this fix for MiSeq cases
		if (grepl("_S", fileName))
		{
			covFileName = sub("^(.+)_S.*-aln-pe-sorted-merged-dedup-recal-covQ20_all_only.hist", "\\1", fileName) 
		}

		# If covFileName is in the samplesheet, then process
		# Later, can pass in analysis type and do the same thing so only coverage of samples from same type is compared
		if(any(sampleTable[c("sample")] == covFileName))
		{
			covTable <- read.table(file, header=FALSE, sep="\t", row.names=NULL)	
			colnames(covTable) <- c("region", "coverage", "numBases", "totalBases", "percentTotal")
		
			covTableRow = 1
	
			covVector <- c()

			i <- 1
	
			while(i <= covTable[1,4])
			{	
				for(j in 1:covTable[covTableRow,3])
				{
					covVector[i] <- covTable[covTableRow,2]
					i = i + 1
				}
				covTableRow = covTableRow + 1
			}
			
			#colnames(covVector) <- c("depth")

			df <- data.frame(
				depth=covVector
			)

			if(grepl(covFileName, sampleID))
			{
				mainTable <- df
			}

			
			if(exists("myplot"))
			{
				myplot <- myplot + geom_histogram(data=df, aes(fill='black', x=depth, y=..density..), binwidth=50, alpha=0.05)
			}
			else
			{
				myplot <- ggplot(df, aes(x=depth)) + geom_histogram(aes(fill='black',y=..density..), binwidth=50, alpha=0.05)
			}

		}

	}

	myplot <- myplot + geom_histogram(data=mainTable, aes(fill='red',x=depth, y=..density..), binwidth=50, color="black", alpha=0.2, show_guide=FALSE) +
		theme_classic() + 
		labs(list(y = "% Bases in ROI", x = "Depth Of Coverage")) + 
		labs(title = paste("Coverage of ", sampleID, " against other samples in run", sep = "")) +
		expand_limits(y = 0) +
		scale_fill_manual(values = c('red' = 'red', 'black' = 'black'), guide = 'legend', labels = c('Other', sampleID), name = '') + 
                scale_x_continuous(expand = c(0, 0), limits=c(0, 1500)) +
                theme(axis.text = element_text(size = 12)) +
                theme(axis.title.x = element_text(size = 16, vjust = -1.5)) +
                theme(axis.title.y = element_text(size = 16, vjust = -0.05)) +
                theme(plot.title = element_text(size = 24, vjust = 3)) +
                theme(legend.text = element_text(size=14)) +
                theme(legend.key.size = unit(0.8, "cm")) +
                theme(plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm"))


	myplot

if (FALSE) 
{

	#maxY = maxY * 1.1

	# Plot
	if(exists("myPlot")) {
		myPlot <- myPlot + geom_line(data = mainTable, aes(x = coverage, y = percentTotal * 100, colour = sampleID), alpha = 1.0) +
		scale_colour_manual("", breaks = c(sampleID, "Other Samples"), values = c("red", "black")) +
		theme_classic() +
		labs(list(y = "% Bases in ROI", x = "Depth of Coverage")) +
		labs(colour = "Legend") +
		scale_y_continuous(expand = c(0, 0), limits=c(0, maxY)) +
		scale_x_continuous(expand = c(0, 0), limits=c(0, 1000)) +
		theme(axis.text = element_text(size = 12)) + 
		theme(axis.title.x = element_text(size = 16, vjust = -1.5)) + 
		theme(axis.title.y = element_text(size = 16, vjust = -0.05)) + 
		theme(plot.title = element_text(size = 24, vjust = 3)) +
		theme(legend.text = element_text(size=14)) + 
		theme(legend.key.size = unit(0.8, "cm")) + 
		theme(plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm"))
	}

	myPlot
}

dev.off()
