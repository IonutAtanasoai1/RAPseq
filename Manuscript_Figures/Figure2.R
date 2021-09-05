---
title: "Figure 2"
output:
  html_document:
    df_print: paged
---

# Authors: Ionut Atanasoai

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
#library(vioplot)
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
library(zoo)



### color codes for colors used in figures


acqua_greens <- c("#04493F","#02655B", "#08756A", "#0C8074", "#179286", "#25A296", "#3EB3A7", "#5FCBC1", "#84E2D9", "#A4F2EA", "#C4FCF6", "#E2FFFD")
acqua_blues  <- c("#02465B","#015666", "#0C6476", "#157589", "#1E8498", "#4397A8", "#4FB0C3", "#70C2D2", "#8DDBEB", "#C1EDF6", "#D7F8FF")
greys <- c("#202020", "#404040", "#606060", "#808080", "#A0A0A0", "#C0C0C0", "#E0E0E0", "#FFFFFF")
reds <- c("#C10303", "#D83535", "#E95F5F", "#F08686", "#FAAEAE")
pinks <- c("#660000", "#990000", "#CC0000", "#FF0000", "#FF6666", "#FF9999","#FFCECE")
yellows <- c( "#9C8803", "#A59213", "#B09D1F", "#BAA82F", "#C6B441", "#D0BF53", "#DACB69", "#E4D782", "#EEE39D", "#FAF2BA", "#FFFBD1")
oranges <- c( "#A15000", "#CE7012", "#E47B12", "#F07D09", "#FF8B15", "#FFA141","#FFCC99")




```




Read in Peak files

```{r}

## read in eCLIP peak files

PATH <- "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/eCLIP_Annotated_PEAKS/"
PEAKS <- list.files(path = PATH)
Names <- unlist(str_split(PEAKS,"\\.p"))[grep("txt",unlist(str_split(PEAKS,"\\.p")), invert = T)]
CLIPs <- list()


for (i in 1:length(Names)){
  
  CLIPs[[i]] <- as.data.frame(read.table(paste(PATH,PEAKS[i],sep = ""), header = T, stringsAsFactors = F))

}

names(CLIPs) <- Names






## read in RAP peak files

PATH <- "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/RAP_Annotated_PEAKS/P102/"
PEAKS <- list.files(path = PATH)
Names <- unlist(str_split(PEAKS,"\\.p"))[grep("txt",unlist(str_split(PEAKS,"\\.p")), invert = T)]
Tables <- list()


for (i in 1:length(Names)){
  
  Tables[[i]] <- as.data.frame(read.table(paste(PATH,PEAKS[i],sep = ""), header = T, stringsAsFactors = F))

}

names(Tables) <- Names


```


Figure 2A
 
```{r}

## the following lines reconstructs the sequence motif logos from PWM published by the relative studies
## the RAP and CLIP PWM are computed by DREME with the parameters specified in the manuscript methods section. For YBX3 in order to obtain the ACAHC motif, one background 4 mers (GAAC, masked to NNNN) needed to be masked. In bash: "cat fasta.file.fa | sed 's/GAAC/NNNN/g' > masked.fasta.file.fa"
## The RBFOX2 GMAUG motif is recomputed by weighing the abundance of GCAUG and GAAUG with the binding scores these motifs have in the peaks they are present
## de novo motif discovery was performed only for RAPseq and eCLIPseq 
## in figure 2A the logos for RBNS were drawn based on the logos published on the ENCODE Project website; PWM were not available



##########################    RNA Compete   ##################

PATH <- "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/Fastas_for_DREME/RNA_Compete_PWMs/"
PEAKS <- list.files(path = PATH)
Names <- unlist(str_split(PEAKS,"\\."))[grep("txt",unlist(str_split(PEAKS,"\\.")),invert = T)]
RNACompete <- list()
  



for (i in 1:length(Names)){
  RNACompete[[i]] <- as.data.frame(read.table(paste(PATH,PEAKS[i],sep = ""), stringsAsFactors = F, header = T))
  RNACompete[[i]] <- t(RNACompete[[i]][,-1])
}

names(RNACompete) <- Names


logos_plots <- list()
for (i in 1:length(Names)){
  logos_plots[[i]] <- ggplot() + geom_logo(RNACompete[[Names[i]]], font="helvetica_bold") + theme_classic() + ylim(0,2)
}
names(logos_plots) <- Names



RNA_Compete <- ggarrange(logos_plots$HNRNPA1,logos_plots$HNRNPC,logos_plots$PTBP1,logos_plots$RBFOX2,logos_plots$YBX1,nrow=length(Names))




##########################    HTR  SELEX   ##################


PATH <- "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/Fastas_for_DREME/HTR-SELEX_PWM/"
PEAKS <- list.files(path = PATH)
Names <- unlist(str_split(PEAKS,"\\."))[grep("txt",unlist(str_split(PEAKS,"\\.")),invert = T)]
HTR_SELEX <- list()
  



for (i in 1:length(Names)){
  HTR_SELEX[[i]] <- as.data.frame(read.table(paste(PATH,PEAKS[i],sep = ""), stringsAsFactors = F, header = F))
  rownames(HTR_SELEX[[i]]) <- HTR_SELEX[[i]][,1]
  colnames(HTR_SELEX[[i]]) <- NULL
  HTR_SELEX[[i]] <- as.matrix(HTR_SELEX[[i]][,-1])
}

names(HTR_SELEX) <- Names
HTR_SELEX$HNRNPA0 <- HTR_SELEX$HNRNPA0[,5:11]
HTR_SELEX$HNRNPC <- HTR_SELEX$HNRNPC[,1:7]
HTR_SELEX$RBFOX1 <- HTR_SELEX$RBFOX1[,8:14]
HTR_SELEX$YBX1 <- HTR_SELEX$YBX1[,3:9]

logos_SELEX <- list()
for (i in 1:length(Names)){
  logos_SELEX[[i]] <- ggplot() + geom_logo(HTR_SELEX[[Names[i]]], font="helvetica_bold") + theme_classic() + ylim(0,2)
}
names(logos_SELEX) <- Names



SELEX <- ggarrange(logos_SELEX$HNRNPA0,logos_SELEX$HNRNPC,NULL,logos_SELEX$RBFOX1,logos_SELEX$YBX1,nrow=length(Names)+1)




############################    RAP SEQ    #################################

###   A   C   G   U   #####


GCATG <- Tables[["RBFOX2"]][grep("GAATG",str_sub(Tables[["RBFOX2"]]$positive_fa,80,120),invert=T),]
GCATG <- Tables[["RBFOX2"]][grep("GCATG",str_sub(Tables[["RBFOX2"]]$positive_fa,80,120)),"BS"]
GCATG <- rep( "GCATG",round(sqrt(length(GCATG))*median(GCATG)) )

GAATG <- Tables[["RBFOX2"]][grep("GCATG",str_sub(Tables[["RBFOX2"]]$positive_fa,80,120),invert=T),]
GAATG <- Tables[["RBFOX2"]][grep("GAATG",str_sub(Tables[["RBFOX2"]]$positive_fa,80,120)),"BS"]
GAATG <- rep( "GAATG",round(sqrt(length(GAATG))*median(GAATG)) )

GMATG <- c(GCATG,GAATG)

GMATG <- DNAStringSet(GMATG)
GMATG <- consensusMatrix(DNAStringSet(GMATG))[1:4,]
rownames(GMATG) <- c("A","C","G","U")
GMAUG <- GMATG
GMAUG <- cbind(GMAUG,c(0,0,0,0),c(0,0,0,0))

GMAUG <- ggplot() + geom_logo(GMAUG, font="helvetica_bold") + theme_classic() + ylim(0,2)



PATH <- "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/Fastas_for_DREME/FINAL_RAP_Motifs/PFMs/"
PEAKS <- list.files(path = PATH)
Names <- unlist(str_split(PEAKS,"\\."))[grep("txt",unlist(str_split(PEAKS,"\\.")),invert = T)]
RAPseq <- list()
  



for (i in 1:length(Names)){
  RAPseq[[i]] <- t(as.data.frame(read.table(paste(PATH,PEAKS[i],sep = ""), stringsAsFactors = F, header = F)))
  rownames(RAPseq[[i]]) <- c("A","C","G","U")
  colnames(RAPseq[[i]]) <- NULL
}

names(RAPseq) <- Names
RAPseq$HNRNPA1 <- cbind(RAPseq$HNRNPA1,c(0,0,0,0),c(0,0,0,0))
RAPseq$HNRNPC <- cbind(RAPseq$HNRNPC,c(0,0,0,0),c(0,0,0,0))
RAPseq$PTBP1 <- cbind(RAPseq$PTBP1,c(0,0,0,0))
RAPseq$YBX3 <- cbind(RAPseq$YBX3,c(0,0,0,0),c(0,0,0,0))

logos_RAPseq <- list()
for (i in 1:length(Names)){
  logos_RAPseq[[i]] <- ggplot() + geom_logo(RAPseq[[Names[i]]], font="helvetica_bold") + theme_classic() + ylim(0,2)
}
names(logos_RAPseq) <- Names




RAP <- ggarrange(logos_RAPseq$HNRNPA1,logos_RAPseq$HNRNPC,logos_RAPseq$PTBP1,GMAUG,logos_RAPseq$YBX3,nrow=length(Names)+1)







############################   eCLIP SEQ    #################################


PATH <- "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/Fastas_for_DREME/FINAL_CLIP_Motifs/PFMs/"
PEAKS <- list.files(path = PATH)
Names <- unlist(str_split(PEAKS,"\\."))[grep("txt",unlist(str_split(PEAKS,"\\.")),invert = T)]
eCLIPseq <- list()
  



for (i in 1:length(Names)){
  eCLIPseq[[i]] <- t(as.data.frame(read.table(paste(PATH,PEAKS[i],sep = ""), stringsAsFactors = F, header = F)))
  rownames(eCLIPseq[[i]]) <- c("A","C","G","U")
  colnames(eCLIPseq[[i]]) <- NULL
}

names(eCLIPseq) <- Names
eCLIPseq$HNRNPA1 <- cbind(eCLIPseq$HNRNPA1,c(0,0,0,0),c(0,0,0,0),c(0,0,0,0))
eCLIPseq$HNRNPC <- cbind(eCLIPseq$HNRNPC,c(0,0,0,0),c(0,0,0,0))
eCLIPseq$PTBP1 <- cbind(eCLIPseq$PTBP1,c(0,0,0,0))
eCLIPseq$RBFOX2 <- cbind(eCLIPseq$RBFOX2,c(0,0,0,0),c(0,0,0,0))
eCLIPseq$YBX3 <- cbind(eCLIPseq$YBX3,c(0,0,0,0),c(0,0,0,0))

logos_eCLIPseq <- list()
for (i in 1:length(Names)){
  logos_eCLIPseq[[i]] <- ggplot() + geom_logo(eCLIPseq[[Names[i]]], font="helvetica_bold") + theme_classic() + ylim(0,2)
}
names(logos_eCLIPseq) <- Names




eCLIP <- ggarrange(logos_eCLIPseq$HNRNPA1,logos_eCLIPseq$HNRNPC,logos_eCLIPseq$PTBP1,logos_eCLIPseq$RBFOX2,logos_eCLIPseq$YBX3,nrow=length(Names))


ggarrange(RAP,eCLIP,RNA_Compete,NULL,SELEX,ncol=5)



```











Figure 2B

```{r, fig.height=9, fig.width=5}


Tables[["HNRNPC"]]$Us <- str_count(str_sub(Tables[["HNRNPC"]]$positive_fa,80,120),"T") / (120-80)
CLIPs[["HNRNPC"]]$Us <- str_count(str_sub(CLIPs[["HNRNPC"]]$positive_fa,43,83),"T") / (83-43)

RAP <- Tables[["HNRNPC"]][,c("Us","BS","Mean_FCI")]
colnames(RAP) <- c("Us","BS","MeanFC")
RAP$Assay <- rep("RAP",nrow(RAP))
CLIP <- CLIPs[["HNRNPC"]][,c("Us","BS","Mean_FC")]
colnames(CLIP) <- c("Us","BS","MeanFC")
CLIP$Assay <- rep("CLIP",nrow(CLIP))
HNRNPC <- rbind(RAP,CLIP)
HNRNPC$Us <- HNRNPC$Us*100
HNRNPC$Us[HNRNPC$Us < 10] <- 10
HNRNPC$Us[HNRNPC$Us > 55] <- 55

a <- ggplot2::ggplot(data=HNRNPC[HNRNPC$Assay == "RAP",], aes(x=Us, y=BS)) + 
  stat_smooth(method="loess", color=acqua_greens[3], fill=acqua_greens[9], level = 0.99) + 
  theme_classic2(base_size = 25) + 
  ylab("Binding Score") + 
  xlab("Uracil content (%)") + 
  scale_x_continuous(breaks = c(10,25,40,55)) 

b <- ggplot2::ggplot(data=HNRNPC[HNRNPC$Assay == "CLIP",], aes(x=Us, y=BS)) + 
  stat_smooth(method="loess", color=greys[1], fill=greys[6], level = 0.99) + 
  theme_classic2(base_size = 25) + 
  ylab("Binding Score") + 
  xlab("Uracil content (%)") + 
  scale_x_continuous(breaks = c(10,25,40,55)) 



ggarrange(a,b,ncol=1)



par(bty="n",mar=c(5,5,1,1), mfrow=c(2,1))
boxplot(data=HNRNPC[HNRNPC$Assay == "RAP",], MeanFC~Us, outline=F, range=1, notch=F, col=acqua_greens[c(10,10,10,9,9,9,8,8,8,7,7,7,6,6,6,5,5,4,4)], ylab="Fold Change", xlab="Uracil content (%)", cex.lab=2, las=1, cex.axis=1.75)
boxplot(data=HNRNPC[HNRNPC$Assay == "CLIP",], MeanFC~Us, outline=F, range=1, notch=F, col=greys[c(7,7,7,6,6,6,5,5,5,4,4,4,3,3,3,2,2,2,1,1)], ylab="Fold Change", xlab="Uracil content (%)", cex.lab=2, las=1, cex.axis=1.75)




```


Figure 2C

```{r, fig.height=10, fig.width=15}



## DeepTools count matrices (bin size = 1nt) are read in
## DeepTools was used with the following command to generate the per nt count matrices
## computeMatrix reference-point -S rep1.bw -R rep1.bed -a 100 -b 100 -bs 1

PATH <- "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/For_DeepTools/Coverages/RAP/"
PEAKS <- list.files(path = PATH)
Names <- unlist(str_split(PEAKS,"\\.matr"))[grep("txt",unlist(str_split(PEAKS,"\\.matr")),invert = T)]
RAPseq <- list()
  



for (i in 1:length(Names)){
  
  RAPseq[[i]] <- as.data.frame(read.table(paste(PATH,PEAKS[i],sep = ""), stringsAsFactors = F))
  
}
names(RAPseq) <- Names



par(mfrow=c(4,5),bty="n")
Names <- c("HNRNPA1","HNRNPC","PTBP1","RBFOX2","YBX3")
for (i in 1:length(Names)){
  
  a <- paste (Names[i],".rep1",sep="")
  b <- paste (Names[i],".rep2",sep="")

  
  a <- rollapply(colMeans(as.matrix(RAPseq[[a]])), width=2, FUN = mean)
  b <- rollapply(colMeans(as.matrix(RAPseq[[b]])), width=2, FUN = mean)

  x <- a/max(c(a,b))
  y <- b/max(c(a,b))

  
  plot(x,type="l",col=acqua_greens[2], lwd=4, ylab=NA, xlab=NA, ylim=c(0,1), las=1, axes=FALSE)
  axis(side=1, labels=T)
  axis(side=2, labels=T)
  points(y,type="l",col=acqua_greens[7], lwd=4)

  
}



Tables[["PTBP1"]]$positive_fa <- gsub("T","Y",Tables[["PTBP1"]]$positive_fa)
Tables[["PTBP1"]]$positive_fa <- gsub("C","Y",Tables[["PTBP1"]]$positive_fa)
Tables[["PTBP1"]]$negative_fa <- gsub("T","Y",Tables[["PTBP1"]]$negative_fa)
Tables[["PTBP1"]]$negative_fa <- gsub("C","Y",Tables[["PTBP1"]]$negative_fa)

Names <- c("HNRNPA1","HNRNPC","PTBP1","RBFOX2","YBX3")
Motifs <- c("TAG","TTTT","YYYY","GCATG|GAATG","AACATC|AACAAC|AACACC")
for (i in 1:length(Names)){
  
  signal <- density(unlist(str_locate_all(Tables[[Names[i]]]$positive_fa,Motifs[i])),bw=2)
  noise <- density(unlist(str_locate_all(Tables[[Names[i]]]$negative_fa,Motifs[i])),bw=2)
  
  plot(signal, col=NA, main=NA, xlab=NA, ylab=NA, ylim=c(0,0.03), axes=FALSE)
  axis(side=1, labels=T)
  axis(side=2, labels=T)
  polygon(signal, col=alpha(acqua_greens[5], 0.9),border = NA)
  polygon(noise, col=alpha(yellows[5], 0.5),border = NA)
  
  
}





## DeepTools count matrices (bin size = 1nt) are read in
## DeepTools was used with the following command to generate the per nt count matrices
## computeMatrix reference-point -S rep1.bw -R rep1.bed -a 100 -b 100 -bs 1


PATH <- "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/For_DeepTools/Coverages/CLIP/"
PEAKS <- list.files(path = PATH)
Names <- unlist(str_split(PEAKS,"\\.matr"))[grep("txt",unlist(str_split(PEAKS,"\\.matr")),invert = T)]
CLIPseq <- list()
  



for (i in 1:length(Names)){
  
  CLIPseq[[i]] <- as.data.frame(read.table(paste(PATH,PEAKS[i],sep = ""), stringsAsFactors = F))
  
}
names(CLIPseq) <- Names




Names <- c("HNRNPA1_eCLIP","HNRNPC_eCLIP","PTBP1_eCLIP","RBFOX2_eCLIP","YBX3_eCLIP")
for (i in 1:length(Names)){
  
  a <- paste (Names[i],".rep1",sep="")
  b <- paste (Names[i],".rep2",sep="")


  
  a <- rollapply(colMeans(as.matrix(CLIPseq[[a]])), width=2, FUN = mean)
  b <- rollapply(colMeans(as.matrix(CLIPseq[[b]])), width=2, FUN = mean)


  x <- a/max(c(a,b))
  y <- b/max(c(a,b))


  
  plot(x,type="l",col=greys[2], lwd=4, ylab=NA, xlab=NA, ylim=c(0,1), las=1, axes=FALSE)
  axis(side=1, labels=T)
  axis(side=2, labels=T)
  points(y,type="l",col=greys[4], lwd=4)


  
}



CLIPs[["PTBP1"]]$positive_fa <- gsub("T","Y",CLIPs[["PTBP1"]]$positive_fa)
CLIPs[["PTBP1"]]$positive_fa <- gsub("C","Y",CLIPs[["PTBP1"]]$positive_fa)
CLIPs[["PTBP1"]]$negative_fa <- gsub("T","Y",CLIPs[["PTBP1"]]$negative_fa)
CLIPs[["PTBP1"]]$negative_fa <- gsub("C","Y",CLIPs[["PTBP1"]]$negative_fa)

Motifs <- c("TAG","TTTT","YYYY","GCATG","CCATC|TCATC")
Names <- c("HNRNPA1","HNRNPC","PTBP1","RBFOX2","YBX3")
for (i in 1:length(Names)){
  
  signal <- density(unlist(str_locate_all(CLIPs[[Names[i]]]$positive_fa,Motifs[i])),bw=2)
  noise <- density(unlist(str_locate_all(CLIPs[[Names[i]]]$negative_fa,Motifs[i])),bw=2)
  
  plot(signal, col=NA, main=NA, xlab=NA, ylab=NA, ylim=c(0,0.03), axes=FALSE)
  axis(side=1, labels=T)
  axis(side=2, labels=T)
  polygon(signal, col=alpha(greys[3], 0.9),border = NA)
  polygon(noise, col=alpha(yellows[5], 0.5),border = NA)
  
  
}








```



Figure 2E

```{r, fig.height=7, fig.width=4}


files <- c("HNRNPA1","HNRNPC","PTBP1","RBFOX2","YBX3")





par(mfrow=c(5,3), mar=c(1,2,1,1))
set.seed(123)
for(i in files){
  a <- log2(Tables[[i]]$FCH_rep1)
  b <- log2(Tables[[i]]$FCH_rep2)
  plot(a,b,xlim=c(0,10),ylim=c(0,10), xlab=NA, ylab=NA, bty="n", xaxt="n", yaxt="n", col = alpha(acqua_greens[4],0.25), cex = 0.6, pch=16)
  axis(1, at = c(0, 5, 10), lwd=1.5, labels = T)
  axis(2, at = c(0, 5, 10), lwd=1.5, labels = T)
  c <- "rho ="
  d <- round(cor(a,b, method = "spearman"), 2)
  cd <- paste(c,d,sep=" ")
  text(3,9.5,labels = d, cex=1.5)
  
  A <- Tables[[i]][sample(1:nrow(Tables[[i]]),(nrow(CLIPs[[i]]))),]
  a <- log2(A$FCH_rep1)
  b <- log2(A$FCH_rep2)
  plot(a,b,xlim=c(0,10),ylim=c(0,10), xlab=NA, ylab=NA, bty="n", xaxt="n", yaxt="n", col = alpha(acqua_greens[7],0.25), cex = 0.6, pch=16)
  axis(1, at = c(0, 5, 10), lwd=1.5, labels = T)
  axis(2, at = c(0, 5, 10), lwd=1.5, labels = T)
  c <- "rho ="
  d <- round(cor(a,b, method = "spearman"), 2)
  cd <- paste(c,d,sep=" ")
  text(3,9.5,labels = d, cex=1.5) 
  
  a <- log2(CLIPs[[i]]$FC_rep1)
  b <- log2(CLIPs[[i]]$FC_rep2)
  plot(a,b,xlim=c(0,10),ylim=c(0,10), xlab=NA, ylab=NA, bty="n", xaxt="n", yaxt="n", col = alpha(greys[2],0.25), cex = 0.6, pch=16)
  axis(1, at = c(0, 5, 10), lwd=1.5, labels = T)
  axis(2, at = c(0, 5, 10), lwd=1.5, labels = T)
  c <- "rho ="
  d <- round(cor(a,b, method = "spearman"), 2)
  cd <- paste(c,d,sep=" ")
  text(3,9.5,labels = d, cex=1.5)
}



RAP <- c()
CLIP <- c()
for(i in files){
  a <- Tables[[i]]$FCH_rep1
  b <- Tables[[i]]$FCH_rep2
  RAP <- c(RAP, cor(a,b,method="spearman"))
  
  a <- CLIPs[[i]]$FC_rep1
  b <- CLIPs[[i]]$FC_rep2
  CLIP <- c(CLIP, cor(a,b,method="spearman"))
}

cors <- c(RAP,CLIP)
assay <- c(rep("aRAP",5),rep("eCLIP",5))

data <- data.frame(cors,assay)
data$cors <- as.numeric(data$cors)
data$assay <- as.factor(data$assay)





set.seed(1)
ggplot(data=data,aes(x=assay,y=cors,fill=assay)) + 
  geom_dotplot(binaxis = "y", 
               stackdir='center', 
               binwidth = .005,
               dotsize = 6,
               position = position_jitter(height=NULL, width = 0.15),
               alpha = 0.9
               ) + 
  theme(panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color="black"),
        axis.ticks = element_blank(),
        axis.text = element_text(color="black"),
        legend.position = "none",
        axis.title = element_blank()
        ) +
  scale_fill_manual(values = c(acqua_greens[4],greys[4]) ) +
  ylim(0.5,1)






```








Figure 2F left panels

```{r, fig.height=3, fig.width=15}



venns <- list()
plot_area <- c()

for (i in 1:length(files)){

AB <- length(intersect(unique(Tables[[files[i]]]$gene_ID),unique(CLIPs[[files[i]]]$gene_ID))) 
A <- length(unique(Tables[[files[i]]]$gene_ID)) - AB
B <- length(unique(CLIPs[[files[i]]]$gene_ID)) - AB
plot_area <- c(plot_area,length(unique(Tables[[files[i]]]$gene_ID))+AB)

v <- euler(c(RAP=A, CLIP=B, "RAP&CLIP"=AB))
venns[[i]] <- plot( v, fills=c(alpha(acqua_greens[i*2],0.7),alpha(greys[i],0.7), yellows[i*2]), main=files[i], quantities=TRUE, edges=F )
}

grid.arrange(grobs = venns, widths=plot_area)





```

Figure 2F right panels

```{r, fig.height=12, fig.width=2.5}





a <- unique(Tables[["HNRNPA1"]][,c("gene_ID","Gene_BS")])
b <- unique(CLIPs[["HNRNPA1"]][,c("gene_ID","Gene_BS_CLIP")])
c <- b[grep(paste(intersect(b$gene_ID,a$gene_ID),collapse="|"),b$gene_ID),]
c$intersection <- rep("common",nrow(c))
d <- b[grep(paste(setdiff(b$gene_ID,a$gene_ID),collapse="|"),b$gene_ID),]
d$intersection <- rep("onlyCLIP",nrow(d))
e <- rbind(c,d)
aa <- ggplot( data = e, aes( x=intersection,y=log2(Gene_BS_CLIP) ) ) +
  geom_violin(aes(fill = intersection), bw=0.85, width=0.8) + 
  scale_fill_manual(values = c(yellows[2],greys[5])) +
  geom_boxplot(color="black", fill=NA, width = 0.2, outlier.shape = NA) + 
  theme_pubr() +
  coord_cartesian(ylim=c(0,13))

HNRNPA1_pval <- wilcox.test(c$Gene_BS_CLIP,d$Gene_BS_CLIP, alternative = "greater")$p.value 

a <- unique(Tables[["HNRNPC"]][,c("gene_ID","Gene_BS")])
b <- unique(CLIPs[["HNRNPC"]][,c("gene_ID","Gene_BS_CLIP")])
c <- b[grep(paste(intersect(b$gene_ID,a$gene_ID),collapse="|"),b$gene_ID),]
c$intersection <- rep("common",nrow(c))
d <- b[grep(paste(setdiff(b$gene_ID,a$gene_ID),collapse="|"),b$gene_ID),]
d$intersection <- rep("onlyCLIP",nrow(d))
e <- rbind(c,d)
bb <- ggplot( data = e, aes( x=intersection,y=log2(Gene_BS_CLIP) ) ) +
  geom_violin(aes(fill = intersection), bw=0.85, width=0.8) + 
  scale_fill_manual(values = c(yellows[2],greys[5])) +
  geom_boxplot(color="black", fill=NA, width = 0.2, outlier.shape = NA) + 
  theme_pubr() +
  coord_cartesian(ylim=c(0,13))

HNRNPC_pval <- wilcox.test(c$Gene_BS_CLIP,d$Gene_BS_CLIP, alternative = "greater")$p.value 




a <- unique(Tables[["PTBP1"]][,c("gene_ID","Gene_BS")])
b <- unique(CLIPs[["PTBP1"]][,c("gene_ID","Gene_BS_CLIP")])
c <- b[grep(paste(intersect(b$gene_ID,a$gene_ID),collapse="|"),b$gene_ID),]
c$intersection <- rep("common",nrow(c))
d <- b[grep(paste(setdiff(b$gene_ID,a$gene_ID),collapse="|"),b$gene_ID),]
d$intersection <- rep("onlyCLIP",nrow(d))
e <- rbind(c,d)
cc <- ggplot( data = e, aes( x=intersection,y=log2(Gene_BS_CLIP) ) ) +
  geom_violin(aes(fill = intersection), bw=0.85, width=0.8) + 
  scale_fill_manual(values = c(yellows[2],greys[5])) +
  geom_boxplot(color="black", fill=NA, width = 0.2, outlier.shape = NA) + 
  theme_pubr() +
  coord_cartesian(ylim=c(0,13))

PTBP1_pval <- wilcox.test(c$Gene_BS_CLIP,d$Gene_BS_CLIP, alternative = "greater")$p.value 


a <- unique(Tables[["RBFOX2"]][,c("gene_ID","Gene_BS")])
b <- unique(CLIPs[["RBFOX2"]][,c("gene_ID","Gene_BS_CLIP")])
c <- b[grep(paste(intersect(b$gene_ID,a$gene_ID),collapse="|"),b$gene_ID),]
c$intersection <- rep("common",nrow(c))
d <- b[grep(paste(setdiff(b$gene_ID,a$gene_ID),collapse="|"),b$gene_ID),]
d$intersection <- rep("onlyCLIP",nrow(d))
e <- rbind(c,d)
dd <- ggplot( data = e, aes( x=intersection,y=log2(Gene_BS_CLIP) ) ) +
  geom_violin(aes(fill = intersection), bw=0.85, width=0.8) + 
  scale_fill_manual(values = c(yellows[2],greys[5])) +
  geom_boxplot(color="black", fill=NA, width = 0.2, outlier.shape = NA) + 
  theme_pubr() +
  coord_cartesian(ylim=c(0,13))

RBFOX2_pval <- wilcox.test(c$Gene_BS_CLIP,d$Gene_BS_CLIP, alternative = "greater")$p.value 


a <- unique(Tables[["YBX3"]][,c("gene_ID","Gene_BS")])
b <- unique(CLIPs[["YBX3"]][,c("gene_ID","Gene_BS_CLIP")])
c <- b[grep(paste(intersect(b$gene_ID,a$gene_ID),collapse="|"),b$gene_ID),]
c$intersection <- rep("common",nrow(c))
d <- b[grep(paste(setdiff(b$gene_ID,a$gene_ID),collapse="|"),b$gene_ID),]
d$intersection <- rep("onlyCLIP",nrow(d))
e <- rbind(c,d)
ee <- ggplot( data = e, aes( x=intersection,y=log2(Gene_BS_CLIP) ) ) +
  geom_violin(aes(fill = intersection), bw=0.85, width=0.8) + 
  scale_fill_manual(values = c(yellows[2],greys[5])) +
  geom_boxplot(color="black", fill=NA, width = 0.2, outlier.shape = NA) + 
  theme_pubr() + 
  coord_cartesian(ylim=c(0,13))

YBX3_pval <- wilcox.test(c$Gene_BS_CLIP,d$Gene_BS_CLIP, alternative = "greater")$p.value 

print("Wilcoxon rank sum test p values")
print(paste("HNRNPA1",HNRNPA1_pval, sep=" "))
print(paste("HNRNPC",HNRNPC_pval,sep=" "))
print(paste("PTBP1",PTBP1_pval,sep=" "))
print(paste("RBFOX2",RBFOX2_pval,sep=" "))
print(paste("YBX3",YBX3_pval,sep=" "))

ggarrange(aa,bb,cc,dd,ee,nrow=5, common.legend = T)




```


Supplementary Figure 2B

```{r, fig.height=15, fig.width=12}

GNWYG <- c("GCATG","GAATG","GGATG","GTTTG","GCTTG","GATTG","GGTTG","GAACG","GCACG","GTATG","GATCG","GCTCG","GTTCG","GGTCG")



GNWYG_scores <- list()
medians <- c()
for (i in 1:length(GNWYG)){

AA <- Tables[["RBFOX2"]][grep(GNWYG[i],str_sub(Tables[["RBFOX2"]]$positive_fa,50,150)),]
aa <- GNWYG[grep(GNWYG[i],GNWYG,invert=T)]

GNWYG_scores[[i]] <- AA[grep(paste(aa,collapse = "|"),str_sub(AA$positive_fa,50,150),invert=T),]$BS



medians <- c(medians,median(GNWYG_scores[[i]]))

}


GNWYG <- gsub("T","U",GNWYG)
names(GNWYG_scores) <- GNWYG




par(bty="n")
boxplot2(GNWYG_scores[order(-medians)], outline=F, horizontal=F, las=3, ylab="Binding Score", range=1, col=c(pinks[2],acqua_blues[c(1:10,10,11,11)]), boxwex=0.95)


GCAUG <- GNWYG_scores$GCAUG
GNWYG_scores <- GNWYG_scores[-1]


for (i in names(GNWYG_scores)){
  aa <- paste("GCAUG vs",i)
  bb <- wilcox.test(GCAUG,GNWYG_scores[[i]], alternative = "greater")$p.value
  cc <- paste("P-value =", bb)
  dd <- paste(aa,cc,sep=" ")
  print(dd)
}




```


Supplementary Figure 2C

```{r}

### have to re-read-in files because for Main Fig2C Cs and Ts were overwritten to Ys for PTBP1

PATH <- "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/RAP_Annotated_PEAKS/P102/"
PEAKS <- list.files(path = PATH)
Names <- unlist(str_split(PEAKS,"\\.p"))[grep("txt",unlist(str_split(PEAKS,"\\.p")), invert = T)]
Tables <- list()
for (i in 1:length(Names)){
  Tables[[i]] <- as.data.frame(read.table(paste(PATH,PEAKS[i],sep = ""), header = T, stringsAsFactors = F))
}
names(Tables) <- Names








GNWYG <- c("GCATG","GAATG","GGATG","GTTTG","GCTTG","GATTG","GGTTG","GAACG","GCACG","GTATG","GATCG","GCTCG","GTTCG","GGTCG")

GNWYG <- GNWYG[order(-medians)]




par(mfrow=c(3,5), mar=c(2,3.5,2,0.5))


for (i in GNWYG){

ymax <- round(nrow(Tables[["RBFOX2"]][grep(i,str_sub(Tables[["RBFOX2"]]$positive_fa,50,150)),]) / 2)
  
hist(unlist(str_locate_all(Tables[["RBFOX2"]]$positive_fa,i))-100, breaks=40, main=gsub("T","U",i), ylab=NA, xlab=NA, col=alpha(acqua_blues[8],0.9), ylim=c(0,ymax), border=NA, xlim=c(-100,100), las=1, xaxt="n", cex.lab=1.5)
axis(side = 1, at = c(-100,-50,0,50,100), labels = c("-100",NA,"0",NA,"100"))
x <- nrow(Tables[["RBFOX2"]][grep(i,str_sub(Tables[["RBFOX2"]]$positive_fa,80,120)),])
n <- nrow(Tables[["RBFOX2"]][grep(i,Tables[["RBFOX2"]]$positive_fa),])

hist(unlist(str_locate_all(Tables[["HNRNPA1"]]$positive_fa,i))-100, breaks=40, type="l",col=alpha(yellows[1],0.25), add=T, border=NA)
n_neg <- nrow(Tables[["HNRNPA1"]][grep(i,Tables[["HNRNPA1"]]$positive_fa),])
x_neg <- nrow(Tables[["HNRNPA1"]][grep(i,str_sub(Tables[["HNRNPA1"]]$positive_fa,80,120)),])
p <- x_neg / n_neg
bt <- binom.test(x,n,p, alternative = "greater")
btp <- bt$p.value
if ( btp <= 0.001 ) {
to_add <- "***" 
} else if ( btp <= 0.01 ) {
to_add <- "**" 
} else if ( btp <= 0.05 ) {
to_add <- "*" 
} else {
to_add <- "n.s." 
}
text(50,ymax/10*9,labels = to_add, cex=2.5, col = yellows[1])


hist(unlist(str_locate_all(Tables[["HNRNPC"]]$positive_fa,i))-100, breaks=40, type="l",col=alpha(yellows[2],0.25), add=T, border=NA)
n_neg <- nrow(Tables[["HNRNPC"]][grep(i,Tables[["HNRNPC"]]$positive_fa),])
x_neg <- nrow(Tables[["HNRNPC"]][grep(i,str_sub(Tables[["HNRNPC"]]$positive_fa,80,120)),])
p <- x_neg / n_neg
bt <- binom.test(x,n,p, alternative = "greater")
btp <- bt$p.value
if ( btp <= 0.001 ) {
to_add <- "***" 
} else if ( btp <= 0.01 ) {
to_add <- "**" 
} else if ( btp <= 0.05 ) {
to_add <- "*" 
} else {
to_add <- "n.s." 
}
text(50,ymax/10*7.75,labels = to_add, cex=2.5, col = yellows[2])



hist(unlist(str_locate_all(Tables[["PTBP1"]]$positive_fa,i))-100, breaks=40, type="l",col=alpha(yellows[3],0.25), add=T, border=NA)
n_neg <- nrow(Tables[["PTBP1"]][grep(i,Tables[["PTBP1"]]$positive_fa),])
x_neg <- nrow(Tables[["PTBP1"]][grep(i,str_sub(Tables[["PTBP1"]]$positive_fa,80,120)),])
p <- x_neg / n_neg
bt <- binom.test(x,n,p, alternative = "greater")
btp <- bt$p.value
if ( btp <= 0.001 ) {
to_add <- "***" 
} else if ( btp <= 0.01 ) {
to_add <- "**" 
} else if ( btp <= 0.05 ) {
to_add <- "*" 
} else {
to_add <- "n.s." 
}
text(50,ymax/10*6.5,labels = to_add, cex=2.5, col = yellows[4])


}





```




Supplementary Figure 2D


```{r, fig.height=7, fig.width=25}


par(mar=c(2,2,4,2),mfrow=c(2,14),bty="n")
GNWYG <- c("GCATG","GAATG","GGATG","GTTTG","GCTTG","GATTG","GGTTG","GAACG","GCACG","GTATG","GATCG","GCTCG","GTTCG","GGTCG")
GNWYG <- GNWYG[order(-medians)]
for (i in GNWYG){
  AA <- Tables[["RBFOX2"]][grep(i,str_sub(Tables[["RBFOX2"]]$positive_fa,50,150)),]
  aa <- GNWYG[grep(i,GNWYG,invert=T)]
  AA <- AA[grep(paste(aa,collapse = "|"),str_sub(AA$positive_fa,50,150),invert=T),]
  aaa <- paste("T",i,sep="")
  TGNWYG <- AA[grep(aaa,str_sub(AA$positive_fa,50,150)),]$BS
  VGNWYG <- AA[grep(aaa,str_sub(AA$positive_fa,50,150), invert = T),]$BS
  NGNWYG <- list(VGNWYG, TGNWYG)
  names(NGNWYG) <- c("V","U") 

  boxplot2(NGNWYG, outline=F, ylab=NA,main=gsub("T","U",i), cex.lab=2, las=1, boxwex=0.95,ylim=c(0,300), col=acqua_blues[c(8,2)])
}



for (i in GNWYG){
  AA <- Tables[["RBFOX2"]][grep(i,str_sub(Tables[["RBFOX2"]]$positive_fa,50,150)),]
  aa <- GNWYG[grep(i,GNWYG,invert=T)]
  AA <- AA[grep(paste(aa,collapse = "|"),str_sub(AA$positive_fa,50,150),invert=T),]
  AA$aa <- str_count(str_sub(AA$positive_fa,50,150),i)
  AA$aa[AA$aa>2] <- 2
  AA$aa <- factor(as.character(AA$aa),levels=c("1","2"))
  
  boxplot2(data=AA,BS~aa, outline=F, ylab=NA,main=gsub("T","U",i), cex.lab=2, las=1, boxwex=0.95, range=1, ylim=c(0,250), col=acqua_blues[c(8,2)])
}


```

Supplemetary Figure 2E

```{r, fig.height=15, fig.width=11}


## the following lines exclude canonical GCAUG as part of the GNWYG to demonstrate that GNWYG acts indipendently of GCAUG 


par(mfrow=c(1,2),bty="n")
GNWYG <- c("GAATG","GGATG","GTTTG","GCTTG","GATTG","GGTTG","GAACG","GCACG","GTATG","GATCG","GCTCG","GTTCG","GGTCG")

aa <- str_count(str_sub(Tables[["RBFOX2"]]$positive_fa,50,150),paste(GNWYG,collapse="|"))
AA <- log2(Tables[["RBFOX2"]][aa>0,"Mean_FCI"])
aa[aa>4] <- 4
aa <- aa[aa>0]
boxplot2(AA~aa,outline=F,boxwex=0.95, ylab="Fold Change (log2)", ylim=c(0,8), xlab="GNWYG without GCATG", col=acqua_blues[c(10,7,4,1)])

 
  x <- AA[aa == 1]
  y <- AA[aa == 2]
  rst <- t.test(x,y,alternative = "less")$p.value
  if ( rst <= 0.001 ) {
  to_add <- "***" 
  } else if ( rst <= 0.01 ) {
  to_add <- "**" 
  } else if ( rst <= 0.05 ) {
  to_add <- "*" 
  } else {
    a <- "P ="
    b <- round(rst,3) 
  to_add <- paste(a,b,sep=" ") 
  }
  text(1.5,7.25,labels = to_add, cex=1.5)
 
 
  x <- AA[aa == 2]
  y <- AA[aa == 3]
  rst <- t.test(x,y,alternative = "less")$p.value
  if ( rst <= 0.001 ) {
  to_add <- "***" 
  } else if ( rst <= 0.01 ) {
  to_add <- "**" 
  } else if ( rst <= 0.05 ) {
  to_add <- "*" 
  } else {
    a <- "P ="
    b <- round(rst,3) 
  to_add <- paste(a,b,sep=" ") 
  }
  text(2.5,7.5,labels = to_add, cex=1.5)

 
  x <- AA[aa == 3]
  y <- AA[aa == 4]
  rst <- t.test(x,y,alternative = "less")$p.value
  if ( rst <= 0.001 ) {
  to_add <- "***" 
  } else if ( rst <= 0.01 ) {
  to_add <- "**" 
  } else if ( rst <= 0.05 ) {
  to_add <- "*" 
  } else {
    a <- "P ="
    b <- round(rst,3) 
  to_add <- paste(a,b,sep=" ") 
  }
  text(3.5,7.75,labels = to_add, cex=1.5)




aa <- str_count(str_sub(CLIPs[["RBFOX2"]]$positive_fa,20,120),paste(GNWYG,collapse="|"))
AA <- CLIPs[["RBFOX2"]][aa>0,"MeanFC"]
aa[aa>4] <- 4
aa <- aa[aa>0]
boxplot2(AA~aa,outline=F,boxwex=0.95, ylab="Fold Change (log2)", ylim=c(2,8), xlab="GNWYG without GCATG", col=greys[c(7,5,3,1)])


 
  x <- AA[aa == 1]
  y <- AA[aa == 2]
  rst <- t.test(x,y,alternative = "less")$p.value
  if ( rst <= 0.001 ) {
  to_add <- "***" 
  } else if ( rst <= 0.01 ) {
  to_add <- "**" 
  } else if ( rst <= 0.05 ) {
  to_add <- "*" 
  } else {
    a <- "P ="
    b <- round(rst,3) 
  to_add <- paste(a,b,sep=" ")
  }
  text(1.5,7.25,labels = to_add, cex=1.5)
 
  x <- AA[aa == 2]
  y <- AA[aa == 3]
  rst <- t.test(x,y,alternative = "less")$p.value
  if ( rst <= 0.001 ) {
  to_add <- "***" 
  } else if ( rst <= 0.01 ) {
  to_add <- "**" 
  } else if ( rst <= 0.05 ) {
  to_add <- "*" 
  } else {
    a <- "P ="
    b <- round(rst,3) 
  to_add <- paste(a,b,sep=" ")
  }
  text(2.5,7.5,labels = to_add, cex=1.5)

 
  x <- AA[aa == 3]
  y <- AA[aa == 4]
  rst <- t.test(x,y,alternative = "less")$p.value
  if ( rst <= 0.001 ) {
  to_add <- "***" 
  } else if ( rst <= 0.01 ) {
  to_add <- "**" 
  } else if ( rst <= 0.05 ) {
  to_add <- "*" 
  } else {
    a <- "P ="
    b <- round(rst,3) 
  to_add <- paste(a,b,sep=" ")
  }
  text(3.5,7.75,labels = to_add, cex=1.5)


```

Supplementary Figure 2F

```{r}



AACAHC <- c("AACATC","AACACC","AACAAC")
ACAHC <- c("ACATC","ACACC","ACAAC")
CAHC <- c("CATC","CACC","CAAC")

AACAHC_scores <- list()
for (i in 1:3){
  AA <- Tables[["YBX3"]][grep(AACAHC[i],str_sub(Tables[["YBX3"]]$positive_fa,50,150)),]
  aa <- CAHC[-i]
  AACAHC_scores[[i]] <- AA[grep(paste(aa,collapse = "|"),str_sub(AA$positive_fa,50,150),invert=T),]$BS
  AACAHC_scores[[i]] <- AA$BS
}
names(AACAHC_scores) <- AACAHC

ACAHC_scores <- list()
for (i in 1:3){
  AA <- Tables[["YBX3"]][grep(AACAHC[i],str_sub(Tables[["YBX3"]]$positive_fa,50,150), invert=T),]
  AA <- AA[grep(ACAHC[i],str_sub(AA$positive_fa,50,150)),]
  aa <- CAHC[-i]
  ACAHC_scores[[i]] <- AA[grep(paste(aa,collapse = "|"),str_sub(AA$positive_fa,50,150),invert=T),]$BS
  ACAHC_scores[[i]] <- AA$BS
}
names(ACAHC_scores) <- ACAHC


CAHC_scores <- list()
for (i in 1:3){
  AA <- Tables[["YBX3"]][grep(AACAHC[i],str_sub(Tables[["YBX3"]]$positive_fa,50,150), invert=T),]
  AA <- AA[grep(ACAHC[i],str_sub(AA$positive_fa,50,150), invert=T),]
  AA <- AA[grep(CAHC[i],str_sub(AA$positive_fa,50,150)),]
  aa <- CAHC[-i]
  CAHC_scores[[i]] <- AA[grep(paste(aa,collapse = "|"),str_sub(AA$positive_fa,50,150),invert=T),]$BS
  CAHC_scores[[i]] <- AA$BS
}
names(CAHC_scores) <- CAHC

ALL_CAHCs_scores <- c(AACAHC_scores[1],ACAHC_scores[1],CAHC_scores[1], AACAHC_scores[2],ACAHC_scores[2],CAHC_scores[2], AACAHC_scores[3],ACAHC_scores[3],CAHC_scores[3])








par(bty="n")
boxplot2(ALL_CAHCs_scores,outline=F, range=1, col=acqua_blues[c(2,6,10)], las=3, ylim=c(0,80), ylab="Binding Score")


 
  x <- ALL_CAHCs_scores$AACATC
  y <- ALL_CAHCs_scores$ACATC
  rst <- t.test(x,y,alternative = "greater")$p.value
  if ( rst <= 0.001 ) {
  to_add <- "***" 
  } else if ( rst <= 0.01 ) {
  to_add <- "**" 
  } else if ( rst <= 0.05 ) {
  to_add <- "*" 
  } else {
  to_add <- "ns"
  }
  text(1.5,77.5,labels = to_add, cex=1.5)

  x <- ALL_CAHCs_scores$ACATC
  y <- ALL_CAHCs_scores$CATC
  rst <- t.test(x,y,alternative = "greater")$p.value
  if ( rst <= 0.001 ) {
  to_add <- "***" 
  } else if ( rst <= 0.01 ) {
  to_add <- "**" 
  } else if ( rst <= 0.05 ) {
  to_add <- "*" 
  } else {
  to_add <- "ns"
  }
  text(2.5,75,labels = to_add, cex=1.5)


   
  x <- ALL_CAHCs_scores$AACACC
  y <- ALL_CAHCs_scores$ACACC
  rst <- t.test(x,y,alternative = "greater")$p.value
  if ( rst <= 0.001 ) {
  to_add <- "***" 
  } else if ( rst <= 0.01 ) {
  to_add <- "**" 
  } else if ( rst <= 0.05 ) {
  to_add <- "*" 
  } else {
  to_add <- "ns"
  }
  text(4.5,77.5,labels = to_add, cex=1.5)

  x <- ALL_CAHCs_scores$ACACC
  y <- ALL_CAHCs_scores$CACC
  rst <- t.test(x,y,alternative = "greater")$p.value
  if ( rst <= 0.001 ) {
  to_add <- "***" 
  } else if ( rst <= 0.01 ) {
  to_add <- "**" 
  } else if ( rst <= 0.05 ) {
  to_add <- "*" 
  } else {
  to_add <- "ns"
  }
  text(5.5,75,labels = to_add, cex=1.5)

  
  
  x <- ALL_CAHCs_scores$AACAAC
  y <- ALL_CAHCs_scores$ACAAC
  rst <- t.test(x,y,alternative = "greater")$p.value
  if ( rst <= 0.001 ) {
  to_add <- "***" 
  } else if ( rst <= 0.01 ) {
  to_add <- "**" 
  } else if ( rst <= 0.05 ) {
  to_add <- "*" 
  } else {
  to_add <- "ns"
  }
  text(7.5,77.5,labels = to_add, cex=1.5)

  x <- ALL_CAHCs_scores$ACAAC
  y <- ALL_CAHCs_scores$CAAC
  rst <- t.test(x,y,alternative = "greater")$p.value
  if ( rst <= 0.001 ) {
  to_add <- "***" 
  } else if ( rst <= 0.01 ) {
  to_add <- "**" 
  } else if ( rst <= 0.05 ) {
  to_add <- "*" 
  } else {
  to_add <- "ns"
  }
  text(8.5,75,labels = to_add, cex=1.5)
  

  




```

Supplementary Figure 2G

```{r,fig.height=10, fig.width=6}


files <- c("HNRNPA1","HNRNPC","PTBP1","RBFOX2","YBX3","IGF2BP1")

a <- table(CLIPs[["RBFOX2"]][CLIPs[["RBFOX2"]]$feature != "exon","feature"])
for (i in files){
  
  b <- table(CLIPs[[i]][CLIPs[[i]]$feature != "exon","feature"])/nrow(CLIPs[[i]][CLIPs[[i]]$feature != "exon",])
  a <- rbind(a,b)
}
a <- a[-1,]
a <- t(a)
colnames(a) <- files




files <- c("HNRNPA1","HNRNPC","PTBP1","RBFOX2","YBX3","IGF2BP1")

c <- table(Tables[["IRP1"]][Tables[["IRP1"]]$feature != "exon","feature"])
for (i in files){
  
  b <- table(Tables[[i]][Tables[[i]]$feature != "exon","feature"])/nrow(Tables[[i]][Tables[[i]]$feature != "exon",])
  c <- rbind(c,b)
}
c <- c[-1,]
c <- t(c)
colnames(c) <- files





par(mfrow=c(2,1))
barplot(a, col=greys[c(1,3,5,7)], las=3, space=0.05, ylab="Fraction of Peaks")
barplot(c, col=acqua_blues[c(1,4,7,11)], las=3, space=0.05, ylab="Fraction of Peaks")


a
c



```


