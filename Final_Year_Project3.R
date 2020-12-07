
library(DiffBind)

SampleList <- read.csv("E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/macs2/SampleList1.csv", header=T)
myDBA <- dba(sampleSheet="E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/macs2/SampleList1.csv")
plot(myDBA)
dba.plotHeatmap(myDBA, correlations = TRUE)

myDBA_count <- dba.count(myDBA, summits=250, minOverlap=2)
dba.plotHeatmap(myDBA_count, correlations = TRUE, colScheme = "Greys")
myDBA_norm <- dba.normalize(myDBA_count)
myDBA_norm


LTTC <- dba(myDBA, mask=myDBA$masks$Tumour)

LTTC

LTTC <- dba.count(LTTC, summits=250, bParallel = 30, minOverlap=2)


LTTC <- dba.normalize(LTTC)
LTTC


LTTC.OL <- dba.overlap(myDBA_consensus, myDBA_consensus$masks$Consensus)
LTTC <- dba.contrast(LTTC, categories=DBA_TISSUE)
LTTC
dba.contrast(LTTC)

LTTC <- dba.analyze(LTTC, bBlacklist=FALSE, bGreylist=FALSE)
dba.report(LTTC)

LTTC.peaks <- dba.report(LTTC, bCalled=TRUE, th=1)


library(org.Hs.eg.db)
library(ChIPpeakAnno)
data(TSS.human.GRCh38)
LTTC.peaks <- annotatePeakInBatch(LTTC.peaks, AnnotationData=TSS.human.GRCh38)
LTTC.peaks <- addGeneIDs(LTTC.peaks, "org.Hs.eg.db", IDs2Add = c("entrez_id", 'symbol'))
LTTC.peaks

write.csv(LTTC.peaks, "E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/macs2/LTTC.peaks2.csv")

UniqueLT <- LTTC.peaks[which(LTTC.peaks$Called1 >=2 & LTTC.peaks$Called2 < 2 & LTTC.peaks$FDR <= 0.1)]
UniqueTC <- LTTC.peaks[which(LTTC.peaks$Called1 < 2 & LTTC.peaks$Called2 >= 2 & LTTC.peaks$FDR <= 0.1)]
DiffPeak <- LTTC.peaks[which(LTTC.peaks$Called1 >=2 & LTTC.peaks$Called2 >= 2 & LTTC.peaks$FDR <= 0.1)]

write.csv(UniqueLT, "E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/macs2/UniqueLT_LTTC2.csv")
write.csv(UniqueTC, "E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/macs2/UniqueTC_LTTC2.csv")
write.csv(DiffPeak, "E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/macs2/DiffPeak_LTTC2.csv")
dba.report(LTTC)

