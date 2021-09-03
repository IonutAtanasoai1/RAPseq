---
title: "Figure 5"
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
library(rrvgo)



### color codes for colors used in figures


acqua_greens <- c("#04493F","#02655B", "#08756A", "#0C8074", "#179286", "#25A296", "#3EB3A7", "#5FCBC1", "#84E2D9", "#A4F2EA", "#C4FCF6", "#E2FFFD")
acqua_blues  <- c("#02465B","#015666", "#0C6476", "#157589", "#1E8498", "#4397A8", "#4FB0C3", "#70C2D2", "#8DDBEB", "#C1EDF6", "#D7F8FF")
greys <- c("#202020", "#404040", "#606060", "#808080", "#A0A0A0", "#C0C0C0", "#E0E0E0", "#FFFFFF")
reds <- c("#C10303", "#D83535", "#E95F5F", "#F08686", "#FAAEAE")
pinks <- c("#660000", "#990000", "#CC0000", "#FF0000", "#FF6666", "#FF9999","#FFCECE")
yellows <- c( "#9C8803", "#A59213", "#B09D1F", "#BAA82F", "#C6B441", "#D0BF53", "#DACB69", "#E4D782", "#EEE39D", "#FAF2BA", "#FFFBD1")
oranges <- c( "#A15000", "#CE7012", "#E47B12", "#F07D09", "#FF8B15", "#FFA141","#FFCC99")

```

FIGURE 5B
```{r, fig.height=4, fig.width=4}

Inputs <- read.table(file="/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/INPUTs.TPMs.txt", stringsAsFactors = F, header=T)


fits <- summary(lm(Inputs$TPM_T7RAP~Inputs$TPM_RAP))
R2 <- as.character(round(fits$adj.r.squared,2))
R2 <- paste("R2 = ", R2, sep="")
spearman <- as.character(round(cor(Inputs$TPM_T7RAP,Inputs$TPM_RAP, method = "spearman"),2))
spearman <- paste("Spearman = ", spearman, sep="")
pearson <- as.character(round(cor(Inputs$TPM_T7RAP,Inputs$TPM_RAP, method = "pearson"),2))
pearson <- paste("Pearson = ", pearson, sep="")
N <- paste("n = ", nrow(Inputs), sep="")

TPM_RAP <- Inputs$TPM_RAP
TPM_T7RAP <- Inputs$TPM_T7RAP


#pdf(file = "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_5/scatter_Inputs.pdf",4,4)
layout.matrix <- matrix(c(2, 1, 0, 3), nrow = 2, ncol = 2)
layout(mat = layout.matrix, heights = c(0.5, 2),  widths = c(2, 0.5)) 

par(mar = c(5, 5, 0, 0))
heatscatter(log2(TPM_T7RAP),log2(TPM_RAP), colpal=greys[2:8], alpha = 80, cex=0.8, bty="l", las=1, xlim=c(0,18), ylim=c(0,18), main="", pch=16)
text(x=4.4,y=16,labels=spearman)
text(x=3.9,y=15,labels=pearson)
text(x=2.5,y=14,labels=R2)
text(x=2.6,y=13,labels=N)
par(mar = c(0.5, 5, 0.5, 0))
d1 <- density(log2(TPM_T7RAP), bw=0.2)
plot(d1$x,d1$y, main=NA, bty="n", xlab=NA, type="l", ylab="density", ylim=c(0,0.3), bty="n", xlim=c(0,18), las=1)
abline(v=median(log2(TPM_T7RAP)))
text(x=5,y=0.29,labels=round(median(log2(TPM_T7RAP)),1))
par(mar = c(5, 0.5, 0, 0.5))
d2 <- density(log2(TPM_RAP), bw=0.2)
plot(d2$y,d2$x, main=NA, bty="n", ylab=NA, type="l", xlab="density", xlim=c(0,0.3), bty="n", ylim=c(0,18), las=2)
abline(h=median(log2(TPM_RAP)))
text(x=0.275,y=5,labels=round(median(log2(TPM_RAP)),1))
#dev.off()

```

FIGURE 5C
```{r, fig.width=3, fig.height=3}

YTHDF1<- read.table(file="/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/RAP_Annotated_PEAKS/YTHDF1.peaks.txt", header=T, stringsAsFactors = F)
YTHDF1$T7Input[YTHDF1$T7Input == 0] <- min(YTHDF1[YTHDF1$T7Input != 0,"T7Input"])
YTHDF1$T7Halo[YTHDF1$T7Halo == 0] <- min(YTHDF1[YTHDF1$T7Halo != 0,"T7Halo"])

RAP_enrichments <- YTHDF1$Peak_ID
T7RAP_enrichments <- YTHDF1[YTHDF1$T7YTHDF1_rep1/YTHDF1$T7Input > 1 & YTHDF1$T7YTHDF1_rep2/YTHDF1$T7Input > 1 & YTHDF1$T7YTHDF1_rep1/YTHDF1$T7Halo > 1 & YTHDF1$T7YTHDF1_rep2/YTHDF1$T7Halo > 1, "Peak_ID"]

AB <- list(RAP_enrichments, T7RAP_enrichments)
names(AB) <- c("with modif", "without modif")
v <- euler(AB, shape="circle")

#pdf(file = "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_5/venn_sites_RAP.vs.T7RAP.pdf",3,3)
plot(v, fills=c(alpha(acqua_greens[3],0.75),alpha(acqua_greens[9],0.9)), quantities=TRUE, edges=F)
#dev.off()

```

FIGURE 5D

```{r, fig.height=5,fig.width=5}

YTHDF1$FCH_T7rep1 <- YTHDF1$T7YTHDF1_rep1 / YTHDF1$T7Halo
YTHDF1$FCH_T7rep2 <- YTHDF1$T7YTHDF1_rep2 / YTHDF1$T7Halo

YTHDF1$FCI_T7rep1 <- YTHDF1$T7YTHDF1_rep1 / YTHDF1$T7Input
YTHDF1$FCI_T7rep2 <- YTHDF1$T7YTHDF1_rep2 / YTHDF1$T7Input


YTHDF1$Mean_FCs_rep1 <- (YTHDF1$FCH_rep1 + YTHDF1$FCI_rep1)/2
YTHDF1$Mean_FCs_rep2 <- (YTHDF1$FCH_rep2 + YTHDF1$FCI_rep2)/2

YTHDF1$Mean_FCs_T7rep1 <- (YTHDF1$FCH_T7rep1 + YTHDF1$FCI_T7rep1)/2
YTHDF1$Mean_FCs_T7rep2 <- (YTHDF1$FCH_T7rep2 + YTHDF1$FCI_T7rep2)/2

YTHDF1$BS_rep1 <- YTHDF1$Mean_FCs_rep1 * log2( (YTHDF1$Rep1)/2 +1)
YTHDF1$BS_rep2 <- YTHDF1$Mean_FCs_rep2 * log2( (YTHDF1$Rep2)/2 +1)
YTHDF1$BS_T7rep1 <- YTHDF1$Mean_FCs_T7rep1 * log2( (YTHDF1$T7YTHDF1_rep1)/2 +1)
YTHDF1$BS_T7rep2 <- YTHDF1$Mean_FCs_T7rep2 * log2( (YTHDF1$T7YTHDF1_rep2)/2 +1)

YTHDF1$BS_T7 <- (YTHDF1$BS_T7rep1 + YTHDF1$BS_T7rep2)/2
YTHDF1$Mean_FCH_T7 <- (YTHDF1$FCH_T7rep1 + YTHDF1$FCH_T7rep2)/2
YTHDF1$Mean_FCI_T7 <- (YTHDF1$FCI_T7rep1 + YTHDF1$FCI_T7rep2)/2
YTHDF1$MeanFC <- (YTHDF1$Mean_FCH + YTHDF1$Mean_FCI)/2
YTHDF1$MeanFC_T7 <- (YTHDF1$Mean_FCH_T7 + YTHDF1$Mean_FCI_T7)/2

FCs <- YTHDF1[,c("BS_rep1","BS_rep2","BS_T7rep1","BS_T7rep2")]
FCs <- melt(FCs, value.name = "Binding_Score", variable.name = "Assay")
A <- rep("RAP", nrow(FCs)/2)
B <- rep("T7RAP", nrow(FCs)/2)
FCs$Assay_Type <- c(A,B)


aa <- ggplot( data=FCs, aes(x=Assay, y=log2(Binding_Score+1)) ) +
  geom_jitter(aes(color=Assay_Type), pch=16, alpha=0.25) +
  scale_color_manual(values = c(acqua_greens[6],greys[6])) +
  theme_classic(base_size = 12.5)

bb <- aa + geom_violin(trim = T, bw=0.75, scale="width", lwd=1, fill=NA) +
  geom_boxplot(outlier.shape = NA, width=0.15, lwd=1, fill=NA) 
  



wt <- wilcox.test(YTHDF1$BS_rep1, YTHDF1$BS_T7rep1, alternative = "greater")
wt <- wt$p.value
if ( wt < 0.001 ) {
ONE <- "***" 
} else if ( wt < 0.01 ) {
ONE <- "**" 
} else if ( wt < 0.05 ) {
ONE <- "*" 
} else {
ONE <- "n.s." 
}

wt <- wilcox.test(YTHDF1$BS_rep1, YTHDF1$BS_T7rep2, alternative = "greater")
wt <- wt$p.value
if ( wt < 0.001 ) {
TWO <- "***" 
} else if ( wt < 0.01 ) {
TWO <- "**" 
} else if ( wt < 0.05 ) {
TWO <- "*" 
} else {
TWO <- "n.s." 
}

wt <- wilcox.test(YTHDF1$BS_rep2, YTHDF1$BS_T7rep1, alternative = "greater")
wt <- wt$p.value
if ( wt < 0.001 ) {
THREE <- "***" 
} else if ( wt < 0.01 ) {
THREE <- "**" 
} else if ( wt < 0.05 ) {
THREE <- "*" 
} else {
THREE <- "n.s." 
}

wt <- wilcox.test(YTHDF1$BS_rep2, YTHDF1$BS_T7rep2, alternative = "greater")
wt <- wt$p.value
if ( wt < 0.001 ) {
FOUR <- "***" 
} else if ( wt < 0.01 ) {
FOUR <- "**" 
} else if ( wt < 0.05 ) {
FOUR <- "*" 
} else {
FOUR <- "n.s." 
}



#pdf(file = "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_5/violins_RAP.vs.T7RAP.pdf",5,5)
bb + geom_text(aes(x=c(2.9),y=c(8.5)), label=ONE, size = 15) +
  geom_text(aes(x=c(4.1),y=c(8.5)), label=TWO, size = 15) +
  geom_text(aes(x=c(2.9),y=c(9.5)), label=THREE, size = 15) +
  geom_text(aes(x=c(4.1),y=c(9.5)), label=FOUR, size = 15) 
#dev.off()

bb

```

FIGURE 5E

```{r, fig.width=3, fig.height=4.5}


DRACH <- c("GGACT", "GAACT", "GGACA", "GAACA", "GGACC", "GAACC", "AAACT", "AGACT", "AAACA", "AGACA", "AAACC", "AGACC", "TGACT", "TAACT")

COUNTS <- c()
for (i in DRACH){
  COUNTS <- c(COUNTS,length(grep(i,str_sub(YTHDF1$positive_fa,85,115))))
}

names(COUNTS) <- DRACH
COUNTS <- COUNTS[order(COUNTS)]

NEG_COUNTS <- c()
for (i in names(COUNTS)){
  NEG_COUNTS <- c(NEG_COUNTS,length(grep(i,str_sub(YTHDF1$negative_fa,85,115))))
}
names(NEG_COUNTS) <- names(COUNTS)

FCs <- COUNTS/NEG_COUNTS
FCs <- FCs[order(FCs)]

#pdf(file = "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_5/barplot_m6A_DRACH.pdf",3,4.5)
barplot(FCs, horiz = T, las=1, xlab="Counts: Bound sites / Control sites", col=acqua_greens[c(1,1,2,2,3,3,4,5,6,7,8,9,10,11)], xlim=c(0,10), space=0.1)
abline(v=c(0,1,2,3,4,5,6,7,8,9,10), lty=2, lwd=1)
#dev.off()


```

FIGURE 5F

```{r, fig.width=4, fig.height=8}



A <- YTHDF1[grep("GGACT|GGAC|GACT",str_sub(YTHDF1$positive_fa,70,130)),]

AB <- list(log2( A$MeanFC+1  ),log2( A$MeanFC_T7+1  ))
AB <- list(A$MeanFC,A$MeanFC_T7)

names(AB) <- c("YES m6A","NO m6A")


#pdf(file = "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_5/boxes_m6A_dep_for_GGACT.pdf",3,5)
par(bty="n")
boxplot2(AB, outline=F, range=1, boxwex=0.35, notch=T, las=2, ylab="Mean Fold Change", ylim=c(0,25), lty=1, main="m6A dependency in GGACT peaks", col=c(acqua_greens[9],greys[6]))
abline(h=c(1,2))
#dev.off()

wt <- wilcox.test(AB$`YES m6A`, AB$`NO m6A`, alternative = "greater")
wt <- wt$p.value
if ( wt < 0.001 ) {
add <- "***" 
} else if ( wt < 0.01 ) {
add <- "**" 
} else if ( wt < 0.05 ) {
add <- "*" 
} else {
add <- "n.s." 
}

text(x=1.5, y=23, labels=add, cex=3)

```

Supplementary FIGURE S5 D

```{r, fig.width=16, fig.height=10}

gene_IDs <- unique(unlist(str_split(unique(YTHDF1$gene_ID),"\\."))[seq(1,length(unique(YTHDF1$gene_ID))*2,2)])
expressed_genes <- unique(unlist(str_split(unique(Inputs$gene_ID),"\\."))[seq(1,length(unique(Inputs$gene_ID))*2,2)])
expressed_genes <- unique(c(expressed_genes, gene_IDs))

ALL_GOs <- enrichGO(
                gene = gene_IDs,
                universe = expressed_genes, 
                keyType       = "ENSEMBL",
                OrgDb         = org.Hs.eg.db,
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                minGSSize     = 15,
                readable      = TRUE
                )
BP <- as.data.frame(ALL_GOs)


MF <- BP[BP$ONTOLOGY == "MF",]
### Collapse child GOs into parent GOs
simMatrix <- calculateSimMatrix(MF$ID,
                                orgdb="org.Hs.eg.db",
                                ont="MF",
                                method="Rel")
scores <- setNames(-log10(MF$qvalue), MF$ID)
reducedTerms_MF <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.7,
                                orgdb="org.Hs.eg.db")
MFs <- reducedTerms_MF[,c("go","parentTerm")]
MFs$Ont <- rep("GO:MF",nrow(MFs))
colnames(MFs) <- c("ID","parentTerm","Ont")

#pdf(file = "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_5/treeMAPplot_MFs.pdf",10,10)
treemapPlot(reducedTerms_MF)
#dev.off()


CC <- BP[BP$ONTOLOGY == "CC",]
### Collapse child GOs into parent GOs
simMatrix <- calculateSimMatrix(CC$ID,
                                orgdb="org.Hs.eg.db",
                                ont="CC",
                                method="Rel")
scores <- setNames(-log10(CC$qvalue), CC$ID)
reducedTerms_CC <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.7,
                                orgdb="org.Hs.eg.db")
CCs <- reducedTerms_CC[,c("go","parentTerm")]
CCs$Ont <- rep("GO:CC",nrow(CCs))
colnames(CCs) <- c("ID","parentTerm","Ont")

#pdf(file = "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_5/treeMAPplot_CCs.pdf",10,10)
treemapPlot(reducedTerms_CC)
#dev.off()


BP <- BP[BP$ONTOLOGY == "BP",]
### Collapse child GOs into parent GOs
simMatrix <- calculateSimMatrix(BP$ID,
                                orgdb="org.Hs.eg.db",
                                ont="BP",
                                method="Rel")
scores <- setNames(-log10(BP$qvalue), BP$ID)
reducedTerms_BP <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.7,
                                orgdb="org.Hs.eg.db")
BPs <- reducedTerms_BP[,c("go","parentTerm")]
BPs$Ont <- rep("GO:BP",nrow(BPs))
colnames(BPs) <- c("ID","parentTerm","Ont")

#pdf(file = "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_5/treeMAPplot_BPs.pdf",10,10)
treemapPlot(reducedTerms_BP)
#dev.off()

GOs <- rbind(BPs,CCs,MFs)
Child_GOs <- as.data.frame(ALL_GOs)
Child_GOs <- Child_GOs[,c("ID","p.adjust","geneID")]

AA <- merge(GOs,Child_GOs,by="ID")

oo <- c("O","O")
pT <- c("pT","pT")
pp <- c("pp","pp")
gg <- c("gg","gg")
CC <- data.frame(oo,pT,pp,gg)
colnames(CC) <- c("Ontology","parentTerm","FDR","gene_name")

for (i in unique(AA$parentTerm)){
  
  gg <- paste(AA[grep(i,AA$parentTerm),"geneID"],collapse = "/")
  gg <- unique(unlist(str_split(gg,"\\/")))
  
  pp <- rep(median(AA[grep(i,AA$parentTerm),"p.adjust"]),length(gg))
  
  pT <- rep(i,length(gg)) 
  
  oo <- unique(AA[grep(i,AA$parentTerm),"Ont"])
  oo <- rep(oo,length(gg))
  
  BB <- data.frame(oo,pT,pp,gg)
  colnames(BB) <- c("Ontology","parentTerm","FDR","gene_name")
  CC <- rbind(CC,BB)
  
}
CC <- CC[-c(1,2),]
DD <- YTHDF1[,c("gene_name","feature","BS","BS_T7")]
CC <- merge(CC,DD,by="gene_name")
CC$BS_T7[CC$BS_T7==0] <- min(CC$BS_T7[CC$BS_T7!=0])




CCC <- CC[,c("parentTerm","BS")]
CCC <- CCC %>% group_by(parentTerm) %>% summarize(Pathway_BS = sum(BS))
CCC <- as.data.frame(CCC)
DDD <- CC[,c("parentTerm","BS_T7")]
DDD <- DDD %>% group_by(parentTerm) %>% summarize(Pathway_BST7 = sum(BS_T7))
DDD <- as.data.frame(DDD)
EEE <- merge(CCC,DDD,by="parentTerm")
FFF <- unique(CC[,c("parentTerm","Ontology")])
EEE <- merge(EEE,FFF,by="parentTerm")
EEE$m6A_dependency <- EEE$Pathway_BS/EEE$Pathway_BST7
GGG <- as.data.frame(table(as.character(CC$parentTerm)))
colnames(GGG) <- c("parentTerm","Binding_Sites")
EEE <- merge(EEE,GGG,by="parentTerm")
GGG <- unique(CC[,c("parentTerm","FDR")])
EEE <- merge(EEE,GGG,by="parentTerm")
EEE$parentTerm <- factor(EEE$parentTerm, levels = EEE[order(EEE$m6A_dependency),"parentTerm"])
EEE$FDR <- as.numeric(EEE$FDR)

ee <- ggplot(data=EEE, aes(x=parentTerm, y=m6A_dependency, fill=FDR)) + 
  geom_bar(stat="identity") +
  coord_flip() +
  theme_pubr() +
  scale_fill_gradient(high=acqua_greens[9], low=acqua_greens[2])

ff <- ggplot(data=EEE, aes(x=parentTerm, y=Binding_Sites, fill=Ontology)) + 
  geom_bar(stat="identity") +
  coord_flip() +
  theme_pubr() +
  scale_fill_manual(values=c(acqua_blues[3],greys[4],yellows[3]))

aa <- ggarrange(ee,ff)

#pdf(file = "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_5/GOs_m6A_dependencies.pdf",16,10)
aa
#dev.off()



```
Figure 5G

```{r, fig.width=12, fig.height=4}

AAA <- unique(CC[grep("stem cell population maintenance|chromatin binding",CC$parentTerm),c("gene_name","parentTerm")])
AAA <- merge(AAA,YTHDF1,by="gene_name")
AAA <- unique(AAA[,c("gene_name","parentTerm","Peak_ID","BS","BS_T7","MeanFC","MeanFC_T7","feature","Rep1","Rep2", "minuslog10pval_rep1", "minuslog10pval_rep2")])
AAA$Signal <- (AAA$Rep1 + AAA$Rep2)/2
AAA$BS_T7[AAA$BS_T7<1] <- 1
AAA$m6A_dependency_sites <- AAA$BS/AAA$BS_T7
AAA <- AAA[AAA$m6A_dependency_sites >= 10,]
AAA$Rank <- rank(-AAA$Signal)
AAA <- AAA[order(AAA$Rank),]
AAA$COLOR <- gsub("chromatin binding",yellows[2],AAA$parentTerm)
AAA$COLOR <- gsub("stem cell population maintenance",acqua_greens[3],AAA$COLOR)
BBB <- head(AAA, 4)



zz <- ggplot( data=AAA,aes(x=Rank,y=log2(Signal), color=parentTerm) ) +
  geom_point(size=7, alpha=0.8) +
  theme_classic() +
  ylab("log2(RAPseq signal)") 
  
yy <- zz + geom_text(data=BBB, aes(x=Rank+10,y=log2(Signal), color=parentTerm, label=paste(gene_name,feature,sep="-")), size=7) +
  scale_color_manual(values=c(yellows[3],acqua_greens[3]))

#pdf(file = "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_5/topGOs_topBSites_rank_dotplot.pdf",12,4)
yy
#dev.off()

BBB

```


Supplementary FIGURE 5B

```{r, fig.width=6, fig.height=4}



YTHDF1$positive_fa_check <- str_sub(YTHDF1$positive_fa,85,115)
YTHDF1$negative_fa_check <- str_sub(YTHDF1$negative_fa)
STRING_posT7 <- as.character(paste(YTHDF1$positive_fa_check,collapse = 'NN'))
STRING_negT7 <- as.character(paste(YTHDF1$negative_fa_check,collapse = 'NN'))
#STRING_pos_YTHDF1 <- as.character(paste(RAP_Tables[["YTHDF1"]]$positive_fa,collapse = 'NN')) 

k=5
bases <- c("A","C","T","G")
kmers <- unite(as.data.frame(permutations(n=4,r=k,v=bases,repeats.allowed = T)),col = kmers,sep = "")

kmers<-cbind(kmers,str_count(STRING_posT7,kmers[,1]))
kmers<-cbind(kmers,str_count(STRING_negT7,kmers[,1]))
#kmers<-cbind(kmers,str_count(STRING_pos_YTHDF1,kmers[,1]))
colnames(kmers) <- c("K","pos","neg")

kmers$pos_fraction <- kmers$pos/sum(kmers$pos)
kmers$neg_fraction <- kmers$neg/sum(kmers$neg)
kmers <- kmers[kmers$pos > 10,]
kmers$Enr <- kmers$pos_fraction/kmers$neg_fraction
kmers$RBP <- rep("YTHDF1",nrow(kmers))
kmers$K_length <- rep("FIVE",nrow(kmers))
kmers <- kmers[(kmers$pos_fraction - kmers$neg_fraction) >0,]
kmers$Norm_enr <- kmers$pos_fraction * kmers$Enr
kmers5 <- kmers

noGGACT <- kmers5[grep("GGACT",kmers5$K,invert = T),]
GGACT <- kmers5[grep("GGACT",kmers5$K),]
GGACT$motif <- "GGACT"
GGAC <- noGGACT[grep("GGAC",noGGACT$K),]
GGAC$motif <- rep("GGAC",nrow(GGAC))
GACT <- noGGACT[grep("GACT",noGGACT$K),]
GACT$motif <- rep("GACT",nrow(GACT))
noMotif <- noGGACT[grep("GACT|GGAC|GGACT",noGGACT$K,invert=T),]
noMotif$motif <- rep("None",nrow(noMotif))
kmers5 <- rbind(GGACT,GGAC,GACT,noMotif)




k=6
bases <- c("A","C","T","G")
kmers <- unite(as.data.frame(permutations(n=4,r=k,v=bases,repeats.allowed = T)),col = kmers,sep = "")

kmers<-cbind(kmers,str_count(STRING_posT7,kmers[,1]))
kmers<-cbind(kmers,str_count(STRING_negT7,kmers[,1]))
#kmers<-cbind(kmers,str_count(STRING_pos_YTHDF1,kmers[,1]))
colnames(kmers) <- c("K","pos","neg")

kmers$pos_fraction <- kmers$pos/sum(kmers$pos)
kmers$neg_fraction <- kmers$neg/sum(kmers$neg)
kmers <- kmers[kmers$pos > 10,]
kmers$Enr <- kmers$pos_fraction/kmers$neg_fraction
kmers$RBP <- rep("YTHDF1",nrow(kmers))
kmers$K_length <- rep("SIX",nrow(kmers))
kmers <- kmers[(kmers$pos_fraction - kmers$neg_fraction) >0,]
kmers$Norm_enr <- kmers$pos_fraction * kmers$Enr
kmers6 <- kmers

noGGACT <- kmers6[grep("GGACT",kmers6$K,invert = T),]
GGACT <- kmers6[grep("GGACT",kmers6$K),]
GGACT$motif <- rep("GGACT",nrow(GGACT))
GGAC <- noGGACT[grep("GGAC",noGGACT$K),]
GGAC$motif <- rep("GGAC",nrow(GGAC))
GACT <- noGGACT[grep("GACT",noGGACT$K),]
GACT$motif <- rep("GACT",nrow(GACT))
noMotif <- noGGACT[grep("GACT|GGAC|GGACT",noGGACT$K,invert=T),]
noMotif$motif <- rep("None",nrow(noMotif))
kmers6 <- rbind(GGACT,GGAC,GACT,noMotif)



k=7
bases <- c("A","C","T","G")
kmers <- unite(as.data.frame(permutations(n=4,r=k,v=bases,repeats.allowed = T)),col = kmers,sep = "")

kmers<-cbind(kmers,str_count(STRING_posT7,kmers[,1]))
kmers<-cbind(kmers,str_count(STRING_negT7,kmers[,1]))
#kmers<-cbind(kmers,str_count(STRING_pos_YTHDF1,kmers[,1]))
colnames(kmers) <- c("K","pos","neg")

kmers$pos_fraction <- kmers$pos/sum(kmers$pos)
kmers$neg_fraction <- kmers$neg/sum(kmers$neg)
kmers <- kmers[kmers$pos > 10,]
kmers$Enr <- kmers$pos_fraction/kmers$neg_fraction
kmers$RBP <- rep("YTHDF1",nrow(kmers))
kmers$K_length <- rep("SEVEN",nrow(kmers))
kmers <- kmers[(kmers$pos_fraction - kmers$neg_fraction) >0,]
kmers$Norm_enr <- kmers$pos_fraction * kmers$Enr
kmers7 <- kmers

noGGACT <- kmers7[grep("GGACT",kmers7$K,invert = T),]
GGACT <- kmers7[grep("GGACT",kmers7$K),]
GGACT$motif <- rep("GGACT",nrow(GGACT))
GGAC <- noGGACT[grep("GGAC",noGGACT$K),]
GGAC$motif <- rep("GGAC",nrow(GGAC))
GACT <- noGGACT[grep("GACT",noGGACT$K),]
GACT$motif <- rep("GACT",nrow(GACT))
noMotif <- noGGACT[grep("GACT|GGAC|GGACT",noGGACT$K,invert=T),]
noMotif$motif <- rep("None",nrow(noMotif))
kmers7 <- rbind(GGACT,GGAC,GACT,noMotif)




kmers <- rbind(kmers5,kmers6,kmers7)
kmers$K_length <- factor(kmers$K_length, levels = c("FIVE","SIX","SEVEN"))



all <- ggplot() + 
  geom_jitter(data = kmers[kmers$motif != "None",], aes(x=K_length, y=Norm_enr, size = log2(Norm_enr+2), color = motif), width = 0.15, alpha = 0.75) +
  scale_color_manual(values = c(yellows[c(2,8)],acqua_blues[1])) +
  ylim(0,0.06)


all1 <- all +
  geom_jitter(data = kmers[kmers$motif == "None",],aes(x=K_length, y=Norm_enr, size = log2(Norm_enr+2)), width = 0.25, alpha = 0.25, color = "grey75") +
    theme(panel.background = element_blank(), 
        panel.grid.major = element_line(size = 0.25, linetype = 2, color=alpha("grey40",0.5)),
        panel.grid.minor = element_line(size = 0.10, linetype = 2, color=alpha("grey40",0.5)),
        axis.text = element_text(size=20,color = "black"),
        panel.border = element_rect(fill=NA,color="black", size = 1),
        axis.ticks.length = unit(3,"mm"),
        axis.ticks = element_line(size = 0.75, color = "black"))

#pdf(file = "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_5/dotplot_kmers.pdf",6,4)
all1
#dev.off()

```

Supplementary FIGURE 5C

```{r, fig.width=4, fig.height=4}

YTHDF1$BS_T7[YTHDF1$BS_T7 == 0] <- min(YTHDF1$BS_T7[YTHDF1$BS_T7 != 0])
YTHDF1$m6A_dependency <- YTHDF1$BS/YTHDF1$BS_T7
YTHDF1$positive_fa_check <- str_sub(YTHDF1$positive_fa,70,130)

GGACT <- YTHDF1[grep("GGACT",YTHDF1$positive_fa_check),c("Peak_ID","positive_fa_check","m6A_dependency","MeanFC","BS", "MeanFC_T7","BS_T7")]
GGACT$motif <- rep("GGACT",nrow(GGACT))
nrow(GGACT)

noGGACT <- YTHDF1[grep("GGACT",YTHDF1$positive_fa_check, invert = T),c("Peak_ID","positive_fa_check","m6A_dependency","MeanFC","BS", "MeanFC_T7","BS_T7")]
noGGACT$motif <- rep("noGGACT",nrow(noGGACT))
nrow(noGGACT)
nrow(noGGACT[grep("GACT|GGAC|GAAC",noGGACT$positive_fa_check),])


GACT <- noGGACT[grep("GACT",noGGACT$positive_fa_check),]
GACT <- GACT[grep("GGAC",GACT$positive_fa,invert = T),]
GACT$motif <- rep("GACT",nrow(GACT))
nrow(GACT)

GGAC <- noGGACT[grep("GGAC",noGGACT$positive_fa_check),]
GGAC <- GGAC[grep("GACT",GGAC$positive_fa,invert = T),]
GGAC$motif <- rep("GGAC",nrow(GGAC))
nrow(GGAC)



all_4mers <- noGGACT[grep("GACT|GGAC",noGGACT$positive_fa_check),]
unique_4mers <- rbind(GACT,GGAC)

two_or_more_4mers <- all_4mers[grep(paste(setdiff(all_4mers$Peak_ID,unique_4mers$Peak_ID),collapse = "|"),all_4mers$Peak_ID),]
two_or_more_4mers$motif <- rep("2orMORE", nrow(two_or_more_4mers))

YTH_motifs <- rbind(GGACT,unique_4mers,two_or_more_4mers)
YTH_motifs$motif <- factor(YTH_motifs$motif, levels = c("GGACT","2orMORE","GACT","GGAC"))



dd <- ggplot(data=YTH_motifs) + 
  stat_ecdf(aes(log2(MeanFC), color=motif), geom = "line", lwd=2) +
  theme(panel.background = element_blank(), 
        panel.grid.major = element_line(size = 0.25, linetype = 2, color=alpha("grey50",0.5)),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size=20,color = "black"),
        panel.border = element_rect(fill=NA,color="black", size = 1),
        axis.ticks.length = unit(3,"mm"),
        axis.ticks = element_line(size = 0.75, color = "black")) +
  scale_color_manual(values = c(acqua_blues[c(2,5)],yellows[c(2,4)])) 





YTH_motifs$FoldChange <- log2(YTH_motifs$MeanFC)
# Kernel density estimate of data
dens1 = density(
  YTH_motifs[YTH_motifs$motif == "GGACT","FoldChange"], 
  bw=0.4
)
dens1 = data.frame(x=dens1$x, y=dens1$y)
# Kernel density estimate of data
dens2 = density(
  YTH_motifs[YTH_motifs$motif == "2orMORE","FoldChange"], 
  bw=0.4
)
dens2 = data.frame(x=dens2$x, y=dens2$y)
# Kernel density estimate of data
dens3 = density(
  YTH_motifs[YTH_motifs$motif == "GACT","FoldChange"], 
  bw=0.4
)
dens3 = data.frame(x=dens3$x, y=dens3$y)
# Kernel density estimate of data
dens4 = density(
  YTH_motifs[YTH_motifs$motif == "GGAC","FoldChange"], 
  bw=0.4
)
dens4 = data.frame(x=dens4$x, y=dens4$y)



#pdf(file = "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_5/ECDFs.pdf",3,4)
par(bty="o")
plot(data=dens1, cumsum(y)/sum(y)~x, lwd=2, col=acqua_blues[2], las=1, xlab="log2FC", type="l")
points(data=dens2, cumsum(y)/sum(y)~x, lwd=2, col=acqua_blues[5], las=1, xlab="log2FC", type="l")
points(data=dens3, cumsum(y)/sum(y)~x, lwd=2, col=yellows[2], las=1, xlab="log2FC", type="l")
points(data=dens4, cumsum(y)/sum(y)~x, lwd=2, col=yellows[4], las=1, xlab="log2FC", type="l")
#dev.off()



wilcox.test(YTH_motifs[YTH_motifs$motif == "GGACT","FoldChange"], YTH_motifs[YTH_motifs$motif == "GACT","FoldChange"], alternative = "greater")
wilcox.test(YTH_motifs[YTH_motifs$motif == "GGACT","FoldChange"], YTH_motifs[YTH_motifs$motif == "GGAC","FoldChange"], alternative = "greater")

wilcox.test(YTH_motifs[YTH_motifs$motif == "2orMORE","FoldChange"], YTH_motifs[YTH_motifs$motif == "GACT","FoldChange"], alternative = "greater")
wilcox.test(YTH_motifs[YTH_motifs$motif == "2orMORE","FoldChange"], YTH_motifs[YTH_motifs$motif == "GGAC","FoldChange"], alternative = "greater")

wilcox.test(YTH_motifs[YTH_motifs$motif == "2orMORE","FoldChange"], YTH_motifs[YTH_motifs$motif == "GGACT","FoldChange"])

```

```{r}
sessionInfo()
```

