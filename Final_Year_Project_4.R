library(DiffBind)
SampleList <- read.csv("E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/macs2/SampleList1.csv", header=T)
myDBA <- dba(sampleSheet="E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/macs2/SampleList1.csv")
myDBA_count <- dba.count(myDBA, summits=250, minOverlap=2)
NLLT <- dba(myDBA, mask=myDBA$masks$Liver)
NLLT
NLLT <- dba.count(NLLT, minOverlap=2, summits=250)
NLLT
NLLT <- dba.normalize(NLLT)
NLLT
NLLT <- dba.contrast(NLLT, categories=DBA_CONDITION)

dba.contrast(NLLT)
NLLT <- dba.analyze(NLLT, bBlacklist=FALSE, bGreylist=FALSE)

NLLT.peaks <- dba.report(NLLT, bCalled=TRUE, th=1)
library(org.Hs.eg.db)
library(ChIPpeakAnno)
data(TSS.human.GRCh38)
NLLT.peaks <- annotatePeakInBatch(NLLT.peaks, AnnotationData=TSS.human.GRCh38)
NLLT.peaks <- addGeneIDs(NLLT.peaks, "org.Hs.eg.db", IDs2Add = c("entrez_id", 'symbol'))
NLLT.peaks

write.csv(LTTC.peaks, "E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/macs2/LTTC.peaks2.csv")

UniqueTC <- LTTC.peaks[which(LTTC.peaks$Called1 >=2 & LTTC.peaks$Called2 < 2 & LTTC.peaks$FDR <= 0.1)]
UniqueLT <- LTTC.peaks[which(LTTC.peaks$Called1 < 2 & LTTC.peaks$Called2 >= 2 & LTTC.peaks$FDR <= 0.1)]
DiffPeak <- LTTC.peaks[which(LTTC.peaks$Called1 >=2 & LTTC.peaks$Called2 >= 2 & LTTC.peaks$FDR <= 0.1)]

write.csv(UniqueLT, "E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/macs2/UniqueLT_LTTC2.csv")
write.csv(UniqueTC, "E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/macs2/UniqueTC_LTTC2.csv")
write.csv(DiffPeak, "E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/macs2/DiffPeak_LTTC2.csv")
