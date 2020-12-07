library("DiffBind")

dirMacs2 = "E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/macs2/"
Samples <- dir("E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/macs2/NCLT/", pattern="_peaks.narrowPeak")
SampleList <- data.frame(SampleID=substr(Samples,1,8),
                         Tissue=substr(Samples,5,6),
                         Factor="",
                         Condition=substr(Samples,5,6),
                         Treatment="5hmC",
                         Replicate="",
                         bamReads=paste0("E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/", substr(Samples,1,12), "-chr19.bam"),
                         ControlID="",
                         bamControl="E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/Merged_Inp.bam",
                         Peaks=paste0("E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/macs2/NCLT/", substr(Samples,1,12) , "-chr19_macs2_peaks.narrowPeak"),
                         PeakCaller="narrow")
SampleList$Tissue[SampleList$Tissue %in% c("LT")] <- "Liver"
SampleList$Tissue[SampleList$Condition %in% c("NC")] <- "Colon"

SampleList$Condition[SampleList$Condition %in% c("LT")] <- "Tumour"
SampleList$Condition[SampleList$Condition %in% c("NC")] <- "Normal"

SampleList$Replicate <- c(1:4, 1:3)

write.csv(SampleList, file="SampleList_NCLT.csv",quote=F, row.names=F)

SampleList <- read.csv("E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/macs2/SampleList_NCLT.csv", header=T)


myDBA <- dba(sampleSheet="E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/macs2/SampleList_NCLT.csv")
plot(myDBA)
dba.plotHeatmap(myDBA, correlations = TRUE)

myDBA

NCLT <- dba(myDBA)

myDBA_consensus <- dba.peakset(NCLT, consensus=c(DBA_TISSUE,DBA_CONDITION), minOverlap=2)
myDBA_consensus <- dba(myDBA_consensus, mask=myDBA_consensus$masks$Consensus, minOverlap=1)
consensus_peaks <- dba.peakset(myDBA_consensus, bRetrieve=TRUE)
dba.plotHeatmap(NCLT, correlations = FALSE, colScheme = "Reds")

NCLT <- dba.count(NCLT, peaks=consensus_peaks, summits=250, bParallel = 30)
dba.contrast(NCLT)


NCLT.OL <- dba.overlap(myDBA_consensus, myDBA_consensus$masks$Consensus)

dba.plotVenn(myDBA_consensus, myDBA_consensus$masks$Consensus)
dba.plotHeatmap(NCLT, correlations =TRUE, colScheme = "Reds")

NCLT <- dba.contrast(NCLT, categories=DBA_CONDITION)
NCLT <- dba.analyze(NCLT)
NCLT.peaks <- dba.report(NCLT, bCalled=TRUE, th=1)

NCLT.peaks

library(org.Hs.eg.db)
library(ChIPpeakAnno)
data(TSS.human.GRCh38)
NCLT.peaks <- annotatePeakInBatch(NCLT.peaks, AnnotationData=TSS.human.GRCh38)
NCLT.peaks <- addGeneIDs(NCLT.peaks, "org.Hs.eg.db", IDs2Add = c("entrez_id", 'symbol'))

NCLT.peaks

write.csv(NCLT.peaks, "E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/macs2/NCLT.peaks9.csv")

UniqueLT <- NCLT.peaks[which(NCLT.peaks$Called1 >=2 & NCLT.peaks$Called2 < 2 & NCLT.peaks$FDR <= 0.1)]
UniqueNC <- NCLT.peaks[which(NCLT.peaks$Called1 < 2 & NCLT.peaks$Called2 >= 2 & NCLT.peaks$FDR <= 0.1)]
DiffPeak <- NCLT.peaks[which(NCLT.peaks$Called1 >=2 & NCLT.peaks$Called2 >= 2 & NCLT.peaks$FDR <= 0.1)]

write.csv(UniqueLT, "E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/macs2/UniqueLT9.csv")
write.csv(UniqueNC, "E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/macs2/UniqueNC9.csv")
write.csv(DiffPeak, "E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/macs2/DiffPeak9.csv")


