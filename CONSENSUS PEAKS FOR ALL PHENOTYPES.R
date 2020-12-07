#####Load packages#####
library(DiffBind)
library(ChIPpeakAnno)
library(org.Hs.eg.db) 

SampleList <- read.csv("E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/macs2/SampleList1.csv", header=T) #Read sample list into R#
myDBA <- dba(sampleSheet="E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/macs2/SampleList1.csv") #create dba (an object within which data can be stored)#
NLLTNCTC <- dba(myDBA) #select for only liver samples#
NLLTNCTC


####Consensus peak and occupancy analysis####
myDBA_consensus_NLLTNCTC <- dba.peakset(NLLTNCTC, consensus=c(DBA_TISSUE,DBA_CONDITION), minOverlap=2) #computes consensus peaks - peaks that appear in at least 2 replicates#
myDBA_consensus_NLLTNCTC <- dba(myDBA_consensus_NLLTNCTC, mask=myDBA_consensus_NLLTNCTC$masks$Consensus, minOverlap=1)
myDBA_consensus_NLLTNCTC #returns number of consensus peaks for Liver Tumour across all LT replicates and Normal Liver across all NL replicates# 
consensus_peaks_NLLTNCTC <- dba.peakset(myDBA_consensus_NLLTNCTC, bRetrieve=TRUE) 
consensus_peaks_NLLTNCTC
myDBA_consensus_NLLTNCTC
dba.plotVenn(myDBA_consensus_NLLTNCTC, myDBA_consensus_NLLTNCTC$masks$Consensus)
NLLTNCTC.count <- dba.count(NLLTNCTC, consensus_peaks_NLLTNCTC, summits = 250) #counts number of peaks for each sample. Only counts consensus peaks - peaks that appear in at least 2 samples more likely to be "real". Summits makes the width of all the peaks the same (250 bp up and downstream of peak summit)
NLLTNCTC.count #displays number of consensus peaks across samples: Samples, 8682 sites in matrix#
NLLTNCTC_normed <- dba.normalize(NLLTNCTC.count) #normalises number of peaks based on lib method - normalises by the mean number of reads across all samples being analysed#
NLLTNCTC_normed

dba.plotHeatmap(NLLTNCTC_normed, correlations = TRUE, colScheme = "Greens")#plot heatmaps of normalised data#
dba.plotHeatmap(NLLTNCTC_normed, correlations = FALSE, colScheme = "Greens")#plot heatmaps of normalised data#
