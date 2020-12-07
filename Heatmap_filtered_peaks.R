
library(DiffBind)
#Create SampleList#
Samples <- dir("E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/macs2/", pattern="_peak2.narrowPeak")
SampleList <- data.frame(SampleID=substr(Samples,1,8),
                         Tissue=substr(Samples,5,6),
                         Factor="",
                         Condition=substr(Samples,5,6),
                         Treatment="5hmC",
                         Replicate="",
                         bamReads=paste0("E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/", substr(Samples,1,12), "-chr19.bam"),
                         ControlID="",
                         bamControl="E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/Merged_Inp.bam",
                         Peaks=paste0("E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/macs2/", substr(Samples,1,12) , "_filtered_peak2.narrowPeak"),
                         PeakCaller="narrow")
SampleList$Tissue[SampleList$Tissue %in% c("LT", "NL")] <- "Liver"
SampleList$Tissue[SampleList$Condition %in% c("NC", "TC")] <- "Colon"
SampleList$Condition[SampleList$Condition %in% c("LT", "TC")] <- "Tumour"
SampleList$Condition[SampleList$Condition %in% c("NL", "NC")] <- "Normal"
SampleList$Replicate <- c(1:4, 1:3, 1:3, 1:3)
write.csv(SampleList, file="E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/Filtered_peaks/SampleList_filteredpeaks4.csv",quote=F, row.names=F)

#Read SampleList into R#
SampleList <- read.csv("E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/Filtered_peaks/SampleList_filteredpeaks4.csv", header=T)

#Create dba (an object in which data can be stored)#
myDBA <- dba(sampleSheet="E:/Final_Year_Project/BathStudent-master/Chr19Fastq/raw-fastq/Fastqc_copy/Filtered_peaks/SampleList_filteredpeaks4.csv")
dba.plotHeatmap(myDBA)
