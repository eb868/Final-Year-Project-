
#####Load packages#####
library(DiffBind)
library(ChIPpeakAnno)
library(org.Hs.eg.db) 

SampleList <- read.csv("E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/macs2/SampleList1.csv", header=T) #Read sample list into R#
myDBA <- dba(sampleSheet="E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/macs2/SampleList1.csv") #create dba (an object within which data can be stored)#
dba.plotHeatmap(myDBA, correlations = TRUE)#plot heatmaps #
hmap <- colorRampPalette(c("white", "grey", "black"))(n = 13)
dba.plotHeatmap(myDBA, correlations = FALSE, colScheme = hmap )#plot heatmaps #

LTTC <- dba(myDBA, mask=myDBA$masks$Tumour) #select for only liver samples#

LTTC #returns number of enriched peaks called for the Liver samples by MACS2#

####Consensus peak and occupancy analysis####
myDBA_consensus_LTTC <- dba.peakset(LTTC, consensus=c(DBA_TISSUE,DBA_CONDITION), minOverlap=2) #computes consensus peaks - peaks that appear in at least 2 replicates#
myDBA_consensus_LTTC <- dba(myDBA_consensus_LTTC, mask=myDBA_consensus_LTTC$masks$Consensus, minOverlap=1)###I DO NOT UNDERSTAND THIS LINE - PLEASE EXPLAIN!###
myDBA_consensus_LTTC #returns number of consensus peaks for Liver Tumour across all LT replicates and Normal Liver across all NL replicates# 
consensus_peaks_LTTC <- dba.peakset(myDBA_consensus_LTTC, bRetrieve=TRUE) 
consensus_peaks_LTTC
dba.plotVenn(myDBA_consensus_LTTC, myDBA_consensus_LTTC$masks$Consensus)


####Differential binding analysis####
LTTC.count <- dba.count(LTTC, consensus_peaks_LTTC, summits = 250) #counts number of peaks for each sample. Only counts consensus peaks - peaks that appear in at least 2 samples more likely to be "real". Summits makes the width of all the peaks the same (250 bp up and downstream of peak summit)
LTTC.count #displays number of consensus peaks across samples: Samples, 8682 sites in matrix#
LTTC_normed <- dba.normalize(LTTC.count) #normalises number of peaks based on lib method - normalises by the mean number of reads across all samples being analysed#
LTTC_normed
LTTC_contrast <- dba.contrast(LTTC_normed, categories=DBA_CONDITION) #Tells diffbind which comparison we are interested in i.e. we want to compare the 2 conditions for Liver: normal and tumour#
LTTC_contrast #displays which group each sample is in for the differential analysis#
LTTC_analyzed <- dba.analyze(LTTC_contrast, bBlacklist = FALSE, bGreylist = FALSE) #Runs default differential binding analysis (statistical sinificance given by p-value) based on negative binomial distribution = DESeq2. Ourdata has already been filtered against a "greylist#
LTTC.peaks <- dba.report(LTTC_analyzed, bCalled = TRUE, th=1) #Shows the differential binding analysis results. Shows which peaks were called in group1 and group2. th=1 will return all peak results (not filtered by a p value)#
data("TSS.human.GRCh38") #Load human genome dataset from org.Hs.eg.db package#
LTTC.peaks <- annotatePeakInBatch(LTTC.peaks, AnnotationData=TSS.human.GRCh38) # match the peaks found with genes in the human genome dataset TSS#
LTTC.peaks <- addGeneIDs(LTTC.peaks, "org.Hs.eg.db", IDs2Add = c("entrez_id", 'symbol')) #annotate the peaks found with gene names#
LTTC.peaks

library(tidyverse)
out <- as.data.frame(LTTC.peaks)
write.csv(out, file="E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/macs2/Results_NLLT/Normal.Liver_vs_Liver.Tumour_desq21.csv")

Normal.Liver_enrich <- out %>%
  filter(FDR < 0.1 & Fold > 0) %>% 
  select(seqnames, start, end, feature, FDR) 
write.csv(Normal.Liver_enrich, file="E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/macs2/Results_NLLT/Normal.Liver_enriched_LTTC.csv")

Liver.Tumour_enrich<- out %>%
  filter(FDR < 0.1 & Fold < 0) %>% 
  select(seqnames, start, end, feature, FDR) 
write.csv(Liver.Tumour_enrich, file="E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/macs2/Results_NLLT/Liver.Tumour_enriched_LTTC.csv")

Similarties.LTTC <-out %>%
  filter(FDR > 0.1) %>%
  select(seqnames, start, end, feature, FDR) 
write.csv(Similarties.LTTC, file="E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/macs2/Results_LTTC/Similarities.LTTC.csv")  
