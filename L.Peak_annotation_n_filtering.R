---
title: "Peak annotation and filtering"
---


```{r}

library(ChIPpeakAnno)
library(edgeR)
library(ggfortify)
library(tidyverse)
library(scales)
library(ggplot2)
library(dendextend)
library(dichromat)
library(reshape2) 
library(dplyr)
library(stringr)
library(GenomicRanges)
library(GenomicFeatures)
library(gplots)
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Dr.eg.db)
library(tidyr)
library(gtools)
library(ggrepel)
library(ggbeeswarm)
library(ReactomePA)
library(VennDiagram)
library(LSD)
library(gridExtra)
library(UpSetR)
library(ggpubr)
library(GenomicScores)
library(phastCons7way.UCSC.hg38)
library(phastCons100way.UCSC.hg38)
library(corrplot)
library(ggforce)
library(idr)
library(GGally)
library(venneuler)
library(eulerr)
library(ggseqlogo)
library(robustbase)




acqua_greens <- c("#04493F","#02655B", "#08756A", "#0C8074", "#179286", "#25A296", "#3EB3A7", "#5FCBC1", "#84E2D9", "#A4F2EA", "#C4FCF6", "#E2FFFD")
acqua_blues  <- c("#02465B","#015666", "#0C6476", "#157589", "#1E8498", "#4397A8", "#4FB0C3", "#70C2D2", "#8DDBEB", "#C1EDF6", "#D7F8FF")
greys <- c("#202020", "#404040", "#606060", "#808080", "#A0A0A0", "#C0C0C0", "#E0E0E0", "#FFFFFF")
reds <- c("#C10303", "#D83535", "#E95F5F", "#F08686", "#FAAEAE")
pinks <- c("#660000", "#990000", "#CC0000", "#FF0000", "#FF6666", "#FF9999","#FFCECE")
yellows <- c( "#9C8803", "#A59213", "#B09D1F", "#BAA82F", "#C6B441", "#D0BF53", "#DACB69", "#E4D782", "#EEE39D", "#FAF2BA", "#FFFBD1")
oranges <- c( "#A15000", "#CE7012", "#E47B12", "#F07D09", "#FF8B15", "#FFA141","#FFCC99")







txdb <- makeTxDbFromGFF(file = "/Users/ionutatanasoai/Documents/RAPseq_Manuscript/Annotations/gencode.v37.annotation.gtf", format = "gtf")
Gencode_v33_IDs <- read.table(file = "/Users/ionutatanasoai/Documents/RAPseq_Manuscript/Annotations/gencode.v37.IDs.txt")
colnames(Gencode_v33_IDs) <- c("gene_ID","transcript_ID","gene_strand","gene_name","gene_type")

####### How the imported "gencode.v37.IDs.txt" file looks like (it is a 5 column tabular file created from parsing the gencode gtf annotation file): #################
####### ENSG00000223972.5	ENST00000456328.2	+	DDX11L1	transcribed_unprocessed_pseudogene
####### ENSG00000223972.5	ENST00000450305.2	+	DDX11L1	transcribed_unprocessed_pseudogene
####### ENSG00000227232.5	ENST00000488147.1	-	WASH7P	unprocessed_pseudogene
####### ENSG00000278267.1	ENST00000619216.1	-	MIR6859-1	miRNA
####### ENSG00000243485.5	ENST00000473358.1	+	MIR1302-2HG	lncRNA
####### ENSG00000243485.5	ENST00000469289.1	+	MIR1302-2HG	lncRNA
####### ENSG00000284332.1	ENST00000607096.1	+	MIR1302-2	miRNA
####### ENSG00000237613.2	ENST00000417324.1	-	FAM138A	lncRNA
####### ENSG00000237613.2	ENST00000461467.1	-	FAM138A	lncRNA
####### ENSG00000268020.3	ENST00000606857.1	+	OR4G4P	unprocessed_pseudogene



Intron_GR <- intronsByTranscript(txdb, use.names = TRUE)
Exon_GR <- exonsBy(txdb, by = "tx", use.names = TRUE)
ThreeUTR_GR <- threeUTRsByTranscript(txdb, use.names = TRUE)
FiveUTR_GR <- fiveUTRsByTranscript(txdb, use.names = TRUE)
CDS_GR <- cdsBy(txdb, by = "tx", use.names = TRUE)
pass_1 <- subsetByOverlaps(Exon_GR,CDS_GR,invert = T)
pass_2 <- subsetByOverlaps(pass_1,ThreeUTR_GR,invert = T)
Exon_GR <- subsetByOverlaps(pass_2,FiveUTR_GR,invert = T)
rm(pass_1)
rm(pass_2)


Introns <- as.data.frame(Intron_GR)[,c(3,4,5,2,7)]
Introns$feature <- rep("intron",nrow(Introns))
Exons <- as.data.frame(Exon_GR)[,c(3,4,5,2,7)]
Exons$feature <- rep("exon",nrow(Exons))
CDSs <- as.data.frame(CDS_GR)[,c(3,4,5,2,7)]
CDSs$feature <- rep("CDS",nrow(CDSs))
FiveUTRs <- as.data.frame(FiveUTR_GR)[,c(3,4,5,2,7)]
FiveUTRs$feature <- rep("5UTR",nrow(FiveUTRs))
ThreeUTRs <- as.data.frame(ThreeUTR_GR)[,c(3,4,5,2,7)]
ThreeUTRs$feature <- rep("3UTR",nrow(ThreeUTRs))
Features <- rbind(Introns, Exons, CDSs, FiveUTRs, ThreeUTRs)
colnames(Features) <- c("chr","start","end","transcript_ID","feature_strand","feature")


Features <- merge(Features,Gencode_v33_IDs,by = "transcript_ID")
Features$gene_ID <- as.character(Features$gene_ID)
Features <- Features[,c(2,3,4,7,6,5,9,10)]
colnames(Features) <- c("chr","start","end","gene_ID","feature","strand","gene_name","gene_type")
Features$IDs <- paste( Features[,1],Features[,2],Features[,3],Features[,4],Features[,5],Features[,6], Features[,7], sep = "_")
Features <- Features[duplicated(Features$IDs) == "FALSE",]
Features <- Features[,1:8]
Features$chr <- as.character(Features$chr)
Features$strand <- as.character(Features$strand)

Features_GR <-  makeGRangesFromDataFrame(Features[,1:6])
colnames(Features) <- c("chr","start","end","gene_ID","feature","gene_strand","gene_name","gene_type")
Features$gene_type <- as.character(Features$gene_type)


```






```{r}




options(scipen = 999,dplyr.summarise.inform = FALSE)

PATH <- "/Users/ionutatanasoai/Documents/RAPseq_Manuscript/NEW_PEAKS/P102/"
PEAKS <- list.files(path = PATH)
Names <- unlist(str_split(PEAKS,"\\.scor"))[grep("bed",unlist(str_split(PEAKS,"\\.scor")),invert = T)]
Tables <- list()
  



for (i in 1:length(Names)){
  Tables[[i]] <- as.data.frame(read.table(paste(PATH,PEAKS[i],sep = ""), stringsAsFactors = F))
  colnames(Tables[[i]]) <- c("chr","Summit_start","Summit_end","start","end","strand","Count_rep1","Count_rep2","minuslog10pval_rep1","minuslog10pval_rep2","Rep1","Rep2","Halo","Input","minusInput_rep1","minusInput_rep2","minusInput_Halo","positive_fa","negative_fa")
}

names(Tables) <- Names

####### the few variables created bellow are required by the idr package to compute the irreproducibility discovery rate #######
## mu <- 2.6
## sigma <- 1.3
## rho <- 0.8
## p <- 0.7


for(i in Names){
  
  ############. If the user wants, the few lines bellow can be used to compute IDRs (both local and global) before any filtering is done #######
  ## x<-Tables[[i]][,c("minuslog10pval_rep1","minuslog10pval_rep2")]
  ## idr.out <- est.IDR(x, mu, sigma, rho, p, eps=0.001, max.ite=30)
  ## Tables[[i]]$local_idr <- idr.out$idr
  ## Tables[[i]]$global_IDR <- idr.out$IDR
  
  Tables[[i]]$RBP <- rep(i,nrow(Tables[[i]]))
  Tables[[i]]$Peak_ID <- paste(  Tables[[i]]$chr,Tables[[i]]$start,Tables[[i]]$end,Tables[[i]]$strand, sep = "_" )
  Tables[[i]]$Summit_start <- Tables[[i]]$Summit_start + 4
  Tables[[i]]$Summit_end <- Tables[[i]]$Summit_start + 1 
  Tables[[i]]$IDs <- paste(  Tables[[i]]$chr,Tables[[i]]$Summit_start,Tables[[i]]$Summit_end,Tables[[i]]$strand, sep = "_" )
  Tables[[i]]$Halo[Tables[[i]]$Halo == 0] <- min(Tables[[i]][Tables[[i]]$Halo!=0,"Halo"])
  Tables[[i]]$Input[Tables[[i]]$Input == 0] <- min(Tables[[i]][Tables[[i]]$Input!=0,"Input"]) 
  Tables[[i]]$minuslog10FDR_rep1 <- round(-log10(p.adjust(10^-Tables[[i]]$minuslog10pval_rep1, method = "BH")),5)
  Tables[[i]]$minuslog10FDR_rep2 <- round(-log10(p.adjust(10^-Tables[[i]]$minuslog10pval_rep2, method = "BH")),5)
  Tables[[i]]$FCH_rep1 <- Tables[[i]]$Rep1/Tables[[i]]$Halo
  Tables[[i]]$FCH_rep2 <- Tables[[i]]$Rep2/Tables[[i]]$Halo
  Tables[[i]]$FCI_rep1 <- Tables[[i]]$Rep1/Tables[[i]]$Input
  Tables[[i]]$FCI_rep2 <- Tables[[i]]$Rep2/Tables[[i]]$Input
  Tables[[i]]$FCmean_rep1 <- (Tables[[i]]$FCH_rep1 + Tables[[i]]$FCI_rep1)/2
  Tables[[i]]$FCmean_rep2 <- (Tables[[i]]$FCH_rep2 + Tables[[i]]$FCI_rep2)/2
  Tables[[i]]$BS_rep1 <- log2(Tables[[i]]$Rep1) * Tables[[i]]$FCmean_rep1
  Tables[[i]]$BS_rep2 <- log2(Tables[[i]]$Rep2) * Tables[[i]]$FCmean_rep2
  Tables[[i]]$BS <- (Tables[[i]]$BS_rep1 + Tables[[i]]$BS_rep2)/2
  Tables[[i]]$Mean_FCH <- (Tables[[i]]$FCH_rep1 + Tables[[i]]$FCH_rep2)/2
  Tables[[i]]$Mean_FCI <- (Tables[[i]]$FCI_rep1 + Tables[[i]]$FCI_rep2)/2
  
  ##### Filtering: FDR <= 0.05; Fold Change: above Halo > 1 and above Input > 1; Fold Change Halo over Input < 2  #############
  Tables[[i]] <- Tables[[i]][Tables[[i]]$minuslog10FDR_rep1 >= 1.30103,]
  Tables[[i]] <- Tables[[i]][Tables[[i]]$minuslog10FDR_rep2 >= 1.30103,]    
  Tables[[i]] <- Tables[[i]][Tables[[i]]$Halo/Tables[[i]]$Input < 2,]
  Tables[[i]] <- Tables[[i]][Tables[[i]]$FCH_rep1 > 1 & Tables[[i]]$FCH_rep2 > 1,]
  Tables[[i]] <- Tables[[i]][Tables[[i]]$FCI_rep1 > 1 & Tables[[i]]$FCI_rep2 > 1,]
  

  ##### Filtering for sequencing complexity to account for spurious and artifactual alignments, GA dinucleotide used for complexity determination ##### 
  Tables[[i]]$positive_fa_check <- str_sub(Tables[[i]]$positive_fa,85,115)
  Gs <- str_count(Tables[[i]]$positive_fa_check,"G") / (115-85) 
  As <- str_count(Tables[[i]]$positive_fa_check,"A") / (115-85) 
  GAs <- Gs + As
  Tables[[i]]$GAs <- GAs
  Tables[[i]] <- Tables[[i]][Tables[[i]]$GAs <= 0.7,]

  
  
  ####### Peak Annotation ########
  GR <- makeGRangesFromDataFrame(Tables[[i]][,c("chr","Summit_start","Summit_end","Peak_ID","RBP","strand")])
    
  Tables[[i]] <- Tables[[i]][as.data.frame(findOverlaps(GR,Features_GR, type = "within"))[,1],]
  Annots <- Features[as.data.frame(findOverlaps(GR,Features_GR, type = "within"))[,2],][,4:8]
  Tables[[i]] <- cbind(Tables[[i]],Annots)
  Tables[[i]]$strand <- as.character(Tables[[i]]$strand)
  Tables[[i]]$gene_strand <- as.character(Tables[[i]]$gene_strand)
  Tables[[i]] <- Tables[[i]][Tables[[i]]$strand == Tables[[i]]$gene_strand, ]
  Tables[[i]]$IDs <- paste( Tables[[i]]$IDs, Tables[[i]]$feature, Tables[[i]]$gene_ID, sep = "_" )
  Tables[[i]] <- Tables[[i]][duplicated(Tables[[i]]$IDs) == "FALSE",]
  Tables[[i]]$IDs <- paste( Tables[[i]]$chr,Tables[[i]]$Summit_start,Tables[[i]]$Summit_end,Tables[[i]]$strand, Tables[[i]]$feature, sep = "_" )
  Tables[[i]] <- Tables[[i]][duplicated(Tables[[i]]$IDs) == "FALSE",]
  
  Tables[[i]]$IDs <- paste( Tables[[i]]$chr,Tables[[i]]$Summit_start,Tables[[i]]$Summit_end,Tables[[i]]$strand, sep = "_" )
  Tables[[i]]$Unique_Anno <-  ave( seq_along(Tables[[i]]$IDs), Tables[[i]]$IDs, FUN = length ) == 1
  Tables[[i]] <- Tables[[i]][(Tables[[i]]$Unique_Anno == "FALSE" & Tables[[i]]$feature == "intron") == "FALSE" ,]
  #print(paste(i,nrow(Tables[[i]])))
  
  Tables[[i]]$Unique_Anno <-  ave( seq_along(Tables[[i]]$IDs), Tables[[i]]$IDs, FUN = length ) == 1
  Tables[[i]] <- Tables[[i]][(Tables[[i]]$Unique_Anno == "FALSE" & Tables[[i]]$feature == "exon") == "FALSE" ,]
  #print(paste(i,nrow(Tables[[i]])))
  
  Tables[[i]]$Unique_Anno <-  ave( seq_along(Tables[[i]]$IDs), Tables[[i]]$IDs, FUN = length ) == 1
  Tables[[i]] <- Tables[[i]][(Tables[[i]]$Unique_Anno == "FALSE" & Tables[[i]]$feature == "5UTR") == "FALSE" ,]
  #print(paste(i,nrow(Tables[[i]])))
  
  Tables[[i]]$Unique_Anno <-  ave( seq_along(Tables[[i]]$IDs), Tables[[i]]$IDs, FUN = length ) == 1
  Tables[[i]] <- Tables[[i]][(Tables[[i]]$Unique_Anno == "FALSE" & Tables[[i]]$feature == "3UTR") == "FALSE" ,]
  #print(paste(i,nrow(Tables[[i]])))
  
  Tables[[i]]$peak_uniqueness <- ave( seq_along(Tables[[i]]$Peak_ID), Tables[[i]]$Peak_ID, FUN = length ) == 1
  #print(paste(i,nrow(Tables[[i]])))
   
  AAA <- Tables[[i]][,c("gene_ID","BS")]
  BBB <- AAA %>% group_by(gene_ID) %>% summarise(Gene_BI = sum(BS))
  BBB <- as.data.frame(BBB)

  Tables[[i]] <- Tables[[i]][grep("pseudogene",Tables[[i]]$gene_type,invert = T),]
  Tables[[i]] <- Tables[[i]][grep("tRNA|rRNA",Tables[[i]]$gene_type,invert = T),]
  #print(paste(i,nrow(Tables[[i]])))
  
  Tables[[i]] <- merge(Tables[[i]],BBB,by = "gene_ID")

}






for(i in Names){
  a <- paste(i,nrow(Tables[[i]]))
  b <- round(cor(Tables[[i]]$FCH_rep1, Tables[[i]]$FCH_rep2,method = "spearman"),2)
  print(paste(a,b,sep=" "))
}


Columns <- c("chr", "start", "end", "Peak_ID", "RBP", "strand", "Summit_start", "Summit_end","Rep1","Rep2","Halo","Input","minuslog10pval_rep1", "minuslog10pval_rep2", "minuslog10FDR_rep1", "minuslog10FDR_rep2", "FCH_rep1", "FCH_rep2", "FCI_rep1", "FCI_rep2", "Mean_FCH", "Mean_FCI", "BS", "gene_ID", "gene_name", "gene_type", "feature", "Gene_BI", "local_idr", "global_IDR", "positive_fa", "negative_fa" )

Columns_rename <- c("chr", "start", "end", "Peak_ID", "RBP", "strand", "Summit_start", "Summit_end","Rep1","Rep2","Halo","Input","minuslog10pval_rep1", "minuslog10pval_rep2", "minuslog10FDR_rep1", "minuslog10FDR_rep2", "FCH_rep1", "FCH_rep2", "FCI_rep1", "FCI_rep2", "Mean_FCH", "Mean_FCI", "BS", "gene_ID", "gene_name", "gene_type", "feature", "Gene_BS", "local_idr", "global_IDR", "positive_fa", "negative_fa" )

for(i in Names){
  Tables[[i]] <- Tables[[i]][,Columns]
  colnames(Tables[[i]]) <- Columns_rename
}



PATH <- "/Users/ionutatanasoai/Documents/RAPseq_Manuscript/RAP_Annotated_PEAKS/"
sapply(names(Tables), function (x) write.table(Tables[[x]], file = paste(paste(PATH,x,sep = ""),".peaks.txt",sep = ""), row.names = F, col.names = T, sep =  "\t") )



```


