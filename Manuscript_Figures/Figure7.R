---
title: "R Notebook"
output:
  html_document:
    df_print: paged
---
## Authors: Sofia Papavasileiou & Ionut Atanasoai
```{r}

library(g3viz)
library(png)
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
library(ggvenn)
library(ggseqlogo)
library(robustbase)
library(rrvgo)
library(scatterpie)


### color codes for colors used in figures


acqua_greens <- c("#04493F","#02655B", "#08756A", "#0C8074", "#179286", "#25A296", "#3EB3A7", "#5FCBC1", "#84E2D9", "#A4F2EA", "#C4FCF6", "#E2FFFD")
acqua_blues  <- c("#02465B","#015666", "#0C6476", "#157589", "#1E8498", "#4397A8", "#4FB0C3", "#70C2D2", "#8DDBEB", "#C1EDF6", "#D7F8FF")
greys <- c("#202020", "#404040", "#606060", "#808080", "#A0A0A0", "#C0C0C0", "#E0E0E0", "#FFFFFF")
reds <- c("#C10303", "#D83535", "#E95F5F", "#F08686", "#FAAEAE")
pinks <- c("#660000", "#990000", "#CC0000", "#FF0000", "#FF6666", "#FF9999","#FFCECE")
yellows <- c( "#9C8803", "#A59213", "#B09D1F", "#BAA82F", "#C6B441", "#D0BF53", "#DACB69", "#E4D782", "#EEE39D", "#FAF2BA", "#FFFBD1")
oranges <- c( "#A15000", "#CE7012", "#E47B12", "#F07D09", "#FF8B15", "#FFA141","#FFCC99")

```


```{r, fig.width=12, fig.height=4.5}


###### data downloaded from https://hive.biochemistry.gwu.edu/biomuta on 2020 - 04 - 12

IGF2BP1_BioMuta <- read.csv(file = "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_7/BioMuta/1.biomuta-proteinview-2020-04-12-15-35-25.csv", header = T, sep = ",")
IGF2BP1_BioMuta <- IGF2BP1_BioMuta[,c(1,6,7,11)]
IGF2BP2_BioMuta <- read.csv(file = "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_7/BioMuta/2.biomuta-proteinview-2020-04-12-15-50-03.csv", header = T, sep = ",")
IGF2BP2_BioMuta <- IGF2BP2_BioMuta[,c(1,6,7,11)]
IGF2BP3_BioMuta <- read.csv(file = "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_7/BioMuta/3.biomuta-proteinview-2020-04-12-15-51-10.csv", header = T, sep = ",")
IGF2BP3_BioMuta <- IGF2BP3_BioMuta[,c(1,6,7,11)]



IGF2BPs_MUTs <- rbind(IGF2BP1_BioMuta,IGF2BP2_BioMuta,IGF2BP3_BioMuta)
IGF2BPs_MUTs$UniProtKB.AC <- gsub("Q9NZI8-1","IGF2BP1",IGF2BPs_MUTs$UniProtKB.AC)
IGF2BPs_MUTs$UniProtKB.AC <- gsub("O00425-1","IGF2BP3",IGF2BPs_MUTs$UniProtKB.AC)
IGF2BPs_MUTs$UniProtKB.AC <- gsub("Q9Y6M1-2","IGF2BP2",IGF2BPs_MUTs$UniProtKB.AC)

IGF2BPs_MUTs <- IGF2BPs_MUTs[,c("UniProtKB.AC","Protein.Position","Frequency")]
colnames(IGF2BPs_MUTs) <- c("gene_name","Protein.Position","Frequency")
IGF2BPs_MUTs <- tibble(IGF2BPs_MUTs)
IGF2BPs_MUTs <- IGF2BPs_MUTs %>% group_by(gene_name,Protein.Position) %>% summarise(Freq = sum(Frequency))


aa <- ggplot(data = IGF2BPs_MUTs) +
  geom_bar(aes(x=Protein.Position, y=Freq), stat="identity", size=0.5, color="black") +  
  geom_point(aes(x=Protein.Position, y=Freq, fill=gene_name), size=6, pch=21, color="black") +
  facet_wrap(~gene_name) +
  scale_fill_manual(values = acqua_greens[c(3,6,9)]) +
  theme_pubr() +
  ylim(0,40)

#pdf(file = "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_7/Biomuta_lollipops.pdf",15,4)
aa
#dev.off()

```





```{r, fig.width=12, fig.height=3}

IGF2BP1s <- read.table(file = "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_7/PEAKS/IGF2BP1s.metafile.txt", header = T, stringsAsFactors = F, sep = "\t")

IGF2BP2s <- read.table(file = "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_7/PEAKS/IGF2BP2s.metafile.txt", header = T, stringsAsFactors = F, sep = "\t")

IGF2BP3s <- read.table(file = "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_7/PEAKS/IGF2BP3s.metafile.txt", header = T, stringsAsFactors = F, sep = "\t")



#pdf(file = "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_7/densities_mutans_paraloggs_Signals.pdf",15,3)
par(mfrow=c(1,3), bty="n")
plot(density(log10(IGF2BP1s$wt_rep1+1),bw=0.1) , xlab="log10(Norm Counts)", ylab="Density (bw=0.1)", col=greys[2], xlim=c(-0.2, 4.5), ylim=c(0,1.5), lwd=3, main="IGF2BP1", las=1, cex.axis=2, cex.lab=2, cex.main=1.5)
points(density(log10(IGF2BP1s$wt_rep1+1),bw=0.1), type="l", lwd=3, col=greys[4])
points(density(log10(IGF2BP1s$R167C_rep1+1),bw=0.1), type="l", lwd=3, col=acqua_greens[3])
points(density(log10(IGF2BP1s$R167C_rep2+1),bw=0.1), type="l", lwd=3, col=acqua_greens[5])
points(density(log10(IGF2BP1s$R167H_rep1+1),bw=0.1), type="l", lwd=3, col=yellows[2])
points(density(log10(IGF2BP1s$R167H_rep2+1),bw=0.1), type="l", lwd=3, col=yellows[4])
text(4,1.25,"WT", col=greys[2], cex=2)
text(4,1,"R167C", col=acqua_greens[4], cex=2)
text(4,0.75,"R167H", col=yellows[2], cex=2)
text(4,0.5,paste("n=",nrow(IGF2BP1s),sep=""),col="black",cex=2)
abline(v=median( log10(c(IGF2BP1s$wt_rep1+1, IGF2BP1s$wt_rep2+1)) ), lwd=3, col=greys[2], lty=2)
abline(v=median( log10(c(IGF2BP1s$R167C_rep1+1, IGF2BP1s$R167C_rep2+1)) ), lwd=3, col=acqua_greens[4], lty=2)
abline(v=median( log10(c(IGF2BP1s$R167H_rep1+1, IGF2BP1s$R167H_rep2+1)) ), lwd=3, col=yellows[2], lty=2)

plot(density(log10(IGF2BP2s$A_rep1+1),bw=0.1) , xlab="log10(Norm Counts)", ylab="Density (bw=0.1)", col=greys[2], xlim=c(-0.2, 4.5), ylim=c(0,1.5), lwd=3, main="IGF2BP2", las=1, cex.axis=2, cex.lab=2, cex.main=1.5)
points(density(log10(IGF2BP2s$A_rep2+1),bw=0.1), type="l", lwd=3, col=greys[4])
points(density(log10(IGF2BP2s$B_rep1+1),bw=0.1), type="l", lwd=3, col=oranges[2])
points(density(log10(IGF2BP2s$B_rep2+1),bw=0.1), type="l", lwd=3, col=oranges[4])
text(4,1.25,"Isoform A", col=greys[2], cex=2)
text(4,1,"Isoform B", col=oranges[3], cex=2)
text(4,0.75,paste("n=",nrow(IGF2BP2s),sep=""),col="black",cex=2)
abline(v=median( log10(c(IGF2BP2s$A_rep1+1, IGF2BP2s$A_rep2+1)) ), lwd=3, col=greys[2], lty=2)
abline(v=median( log10(c(IGF2BP2s$B_rep1+1, IGF2BP2s$B_rep2+1)) ), lwd=3, col=oranges[3], lty=2)

plot(density(log10(IGF2BP3s$wt_rep1+1),bw=0.1) , xlab="log10(Norm Counts)", ylab="Density (bw=0.1)", col=greys[2], xlim=c(-0.2, 4.5), ylim=c(0,1.5), lwd=3, main="IGF2BP3", las=1, cex.axis=2, cex.lab=2, cex.main=1.5)
points(density(log10(IGF2BP3s$wt_rep2+1),bw=0.1), type="l", lwd=3, col=greys[4])
points(density(log10(IGF2BP3s$R525C_rep1+1),bw=0.1), type="l", lwd=3, col=acqua_blues[5])
points(density(log10(IGF2BP3s$R525C_rep2+1),bw=0.1), type="l", lwd=3, col=acqua_blues[7])
points(density(log10(IGF2BP3s$I474M_rep1+1),bw=0.1), type="l", lwd=3, col=pinks[2])
points(density(log10(IGF2BP3s$I474M_rep2+1),bw=0.1), type="l", lwd=3, col=pinks[3])
text(4,1.25,"WT", col=greys[2], cex=2)
text(4,1,"R525C", col=acqua_blues[5], cex=2)
text(4,0.75,"I474M", col=pinks[3], cex=2)
text(4,0.5,paste("n=",nrow(IGF2BP3s),sep=""),col="black",cex=2)
abline(v=median( log10(c(IGF2BP3s$wt_rep1+1, IGF2BP3s$wt_rep2+1)) ), lwd=3, col=greys[2], lty=2)
abline(v=median( log10(c(IGF2BP3s$R525C_rep1+1, IGF2BP3s$R525C_rep2+1)) ), lwd=3, col=acqua_blues[5], lty=2)
abline(v=median( log10(c(IGF2BP3s$I474M_rep1+1, IGF2BP3s$I474M_rep2+1)) ), lwd=3, col=pinks[3], lty=2)
#dev.off()



```



```{r, fig.width=16, fig.height=3.75}

SUB <- IGF2BP1s[,7:12] 
SUB <- SUB[SUB$wt_rep1 >= 5 & SUB$wt_rep2 >= 5 & SUB$R167C_rep1 >= 5 & SUB$R167C_rep2 >= 5 & SUB$R167H_rep1 >= 5 & SUB$R167H_rep2 >= 5,]
nrow(SUB)

WTvsRtoC <- SUB[,c("wt_rep1", "wt_rep2", "R167C_rep1", "R167C_rep2")] 
x <- WTvsRtoC
group <- c(1,1,2,2)
y <- DGEList(counts=x, group = group)
design <- model.matrix(~group)
y <- estimateDisp(y,design)
y$samples$lib.size <- c(1000000,1000000,1000000,1000000)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit)
X <-qlf$table
X$p.adjust <- -log10(p.adjust(X$PValue, method = "BH"))
X$PValue <- -log10(X$PValue)
WTvsRtoC$FDR <- X$p.adjust
WTvsRtoC$Significant <- WTvsRtoC$FDR >= 1.30103
WTvsRtoC$log2FC <- X$logFC
WTvsRtoC$log2CPM <- X$logCPM
WTvsRtoC$wt_won <- WTvsRtoC$log2FC <= -1 & WTvsRtoC$FDR >= 1.30103
WTvsRtoC$RtoC_won <- WTvsRtoC$log2FC >= 1 & WTvsRtoC$FDR >= 1.30103
WTvsRtoC$DiffBind <- paste(WTvsRtoC$wt_won,WTvsRtoC$RtoC_won,sep = "_")
WTvsRtoC$DiffBind <- gsub("FALSE_FALSE","N.S.",WTvsRtoC$DiffBind)
WTvsRtoC$DiffBind <- gsub("FALSE_TRUE","R167C",WTvsRtoC$DiffBind)
WTvsRtoC$DiffBind <- gsub("TRUE_FALSE","wt",WTvsRtoC$DiffBind)

WTvsRtoH <- SUB[,c("wt_rep1", "wt_rep2", "R167H_rep1", "R167H_rep2")] 
x <- WTvsRtoH
group <- c(1,1,2,2)
y <- DGEList(counts=x, group = group)
design <- model.matrix(~group)
y <- estimateDisp(y,design)
y$samples$lib.size <- c(1000000,1000000,1000000,1000000)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit)
X <-qlf$table
X$p.adjust <- -log10(p.adjust(X$PValue, method = "BH"))
X$PValue <- -log10(X$PValue)
WTvsRtoH$FDR <- X$p.adjust
WTvsRtoH$Significant <- WTvsRtoH$FDR >= 1.30103
WTvsRtoH$log2FC <- X$logFC
WTvsRtoH$log2CPM <- X$logCPM
WTvsRtoH$wt_won <- WTvsRtoH$log2FC <= -1 & WTvsRtoH$FDR >= 1.30103
WTvsRtoH$RtoH_won <- WTvsRtoH$log2FC >= 1 & WTvsRtoH$FDR >= 1.30103
WTvsRtoH$DiffBind <- paste(WTvsRtoH$wt_won,WTvsRtoH$RtoH_won,sep = "_")
WTvsRtoH$DiffBind <- gsub("FALSE_FALSE","N.S.",WTvsRtoH$DiffBind)
WTvsRtoH$DiffBind <- gsub("FALSE_TRUE","R167H",WTvsRtoH$DiffBind)
WTvsRtoH$DiffBind <- gsub("TRUE_FALSE","wt",WTvsRtoH$DiffBind)
AA <- WTvsRtoC[,c(5,7,8,11)]
AA$comparison <- rep("wtBP1 vs R167C", nrow(AA))
BB <- WTvsRtoH[,c(5,7,8,11)]
BB$comparison <- rep("wtBP1 vs R167H", nrow(BB))
COMPARISONS_1 <- rbind(AA,BB)
colnames(COMPARISONS_1) <- c("minuslog10FDR", "log2FC", "log2CPM", "DiffBind", "comparison")




SUB <- IGF2BP2s[,7:10] 
SUB <- SUB[SUB$A_rep1 >= 5 & SUB$A_rep2 >= 5 & SUB$B_rep1 >= 5 & SUB$B_rep2 >= 5 ,]
nrow(SUB)

AvsB <- SUB
x <- AvsB
group <- c(1,1,2,2)
y <- DGEList(counts=x, group = group)
design <- model.matrix(~group)
y <- estimateDisp(y,design)
y$samples$lib.size <- c(1000000,1000000,1000000,1000000)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit)
X <-qlf$table
X$p.adjust <- -log10(p.adjust(X$PValue, method = "BH"))
X$PValue <- -log10(X$PValue)
AvsB$FDR <- X$p.adjust
AvsB$Significant <- AvsB$FDR >= 1.30103
AvsB$log2FC <- X$logFC
AvsB$log2CPM <- X$logCPM
AvsB$wt_won <- AvsB$log2FC <= -1 & AvsB$FDR >= 1.30103
AvsB$RtoC_won <- AvsB$log2FC >= 1 & AvsB$FDR >= 1.30103
AvsB$DiffBind <- paste(AvsB$wt_won,AvsB$RtoC_won,sep = "_")
AvsB$DiffBind <- gsub("FALSE_FALSE","N.S.",AvsB$DiffBind)
AvsB$DiffBind <- gsub("FALSE_TRUE","Isoform_B",AvsB$DiffBind)
AvsB$DiffBind <- gsub("TRUE_FALSE","Isoform_A",AvsB$DiffBind)
AvsB$comparison <- rep("Isof_A vs Isof_B", nrow(AvsB))
COMPARISONS_2 <- AvsB[,c(5,7,8,11,12)]
colnames(COMPARISONS_2) <- c("minuslog10FDR", "log2FC", "log2CPM", "DiffBind", "comparison")
  







SUB <- IGF2BP3s[,7:12] 
SUB <- SUB[SUB$wt_rep1 >= 5 & SUB$wt_rep2 >= 5 & SUB$R525C_rep1 >= 5 & SUB$R525C_rep2 >= 5 & SUB$I474M_rep1 >= 5 & SUB$I474M_rep2 >= 5,]
nrow(SUB)

WTvsR525C <- SUB[,c("wt_rep1", "wt_rep2", "R525C_rep1", "R525C_rep2")] 
x <- WTvsR525C
group <- c(1,1,2,2)
y <- DGEList(counts=x, group = group)
design <- model.matrix(~group)
y <- estimateDisp(y,design)
y$samples$lib.size <- c(1000000,1000000,1000000,1000000)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit)
X <-qlf$table
X$p.adjust <- -log10(p.adjust(X$PValue, method = "BH"))
X$PValue <- -log10(X$PValue)
WTvsR525C$FDR <- X$p.adjust
WTvsR525C$Significant <- WTvsR525C$FDR >= 1.30103
WTvsR525C$log2FC <- X$logFC
WTvsR525C$log2CPM <- X$logCPM
WTvsR525C$wt_won <- WTvsR525C$log2FC <= -1 & WTvsR525C$FDR >= 1.30103
WTvsR525C$R525C_won <- WTvsR525C$log2FC >= 1 & WTvsR525C$FDR >= 1.30103
WTvsR525C$DiffBind <- paste(WTvsR525C$wt_won,WTvsR525C$R525C_won,sep = "_")
WTvsR525C$DiffBind <- gsub("FALSE_FALSE","N.S.",WTvsR525C$DiffBind)
WTvsR525C$DiffBind <- gsub("FALSE_TRUE","R525C",WTvsR525C$DiffBind)
WTvsR525C$DiffBind <- gsub("TRUE_FALSE","wt",WTvsR525C$DiffBind)

WTvsI474M <- SUB[,c("wt_rep1", "wt_rep2", "I474M_rep1", "I474M_rep2")] 
x <- WTvsI474M
group <- c(1,1,2,2)
y <- DGEList(counts=x, group = group)
design <- model.matrix(~group)
y <- estimateDisp(y,design)
y$samples$lib.size <- c(1000000,1000000,1000000,1000000)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit)
X <-qlf$table
X$p.adjust <- -log10(p.adjust(X$PValue, method = "BH"))
X$PValue <- -log10(X$PValue)
WTvsI474M$FDR <- X$p.adjust
WTvsI474M$Significant <- WTvsI474M$FDR >= 1.30103
WTvsI474M$log2FC <- X$logFC
WTvsI474M$log2CPM <- X$logCPM
WTvsI474M$wt_won <- WTvsI474M$log2FC <= -1 & WTvsI474M$FDR >= 1.30103
WTvsI474M$I474M_won <- WTvsI474M$log2FC >= 1 & WTvsI474M$FDR >= 1.30103
WTvsI474M$DiffBind <- paste(WTvsI474M$wt_won,WTvsI474M$I474M_won,sep = "_")
WTvsI474M$DiffBind <- gsub("FALSE_FALSE","N.S.",WTvsI474M$DiffBind)
WTvsI474M$DiffBind <- gsub("FALSE_TRUE","I474M",WTvsI474M$DiffBind)
WTvsI474M$DiffBind <- gsub("TRUE_FALSE","wt",WTvsI474M$DiffBind)
AA <- WTvsR525C[,c(5,7,8,11)]
AA$comparison <- rep("wtBP3 vs R525C", nrow(AA))
BB <- WTvsI474M[,c(5,7,8,11)]
BB$comparison <- rep("wtBP3 vs I474M", nrow(BB))
COMPARISONS_3 <- rbind(AA,BB)
colnames(COMPARISONS_3) <- c("minuslog10FDR", "log2FC", "log2CPM", "DiffBind", "comparison")



ALL_COMPARISONS <- rbind(COMPARISONS_1, COMPARISONS_2, COMPARISONS_3)
ALL_COMPARISONS$comparison <- factor(ALL_COMPARISONS$comparison, levels=c("wtBP1 vs R167C","wtBP1 vs R167H","Isof_A vs Isof_B","wtBP3 vs I474M","wtBP3 vs R525C"))

table(ALL_COMPARISONS[,c("DiffBind","comparison")])

dd <- ggplot(data=ALL_COMPARISONS) +
  geom_point(aes(x=log2CPM, y=log2FC, color=DiffBind), size=3) +
  scale_color_manual(values = c( pinks[2],greys[2],oranges[3],alpha(greys[6],0.3),yellows[2],acqua_blues[5],greys[2]  ) ) +
  theme_pubr() +
  geom_hline(yintercept = 0, size=1, color="black", lty=2) +
  ylim(-4,4) +
  xlim(2,15) + 
  facet_wrap(~comparison, ncol=5)
dd

#pdf(file = "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_7/diff_bind_MDplots.pdf",15,4.28)
#dd
#dev.off()









```



```{r, fig.width=10, fig.height=2}




SUB <- IGF2BP2s
SUB$feature <- gsub("exon","ncRNA",SUB$feature)
SUB <- SUB[SUB$A_rep1 >= 5 & SUB$A_rep2 >= 5 & SUB$B_rep1 >= 5 & SUB$B_rep2 >= 5 ,]
SUB$DiffBind <- AvsB$DiffBind


options(scipen=99)
aa <- (t(table(SUB[,c("feature","DiffBind")])) / rowSums(t(table(SUB[,c("feature","DiffBind")])))) * 100
aa <- aa[1:2,c(1,5)]
aa
chisq.test(aa)




#pdf(file = "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_7/barplots_features_for_diffboundsites.pdf",15,3)


par(mfrow=c(1,5))
SUB <- IGF2BP1s
SUB$feature <- gsub("exon","ncRNA",SUB$feature)
SUB <- SUB[SUB$wt_rep1 >= 5 & SUB$wt_rep2 >= 5 & SUB$R167C_rep1 >= 5 & SUB$R167C_rep2 >= 5 & SUB$R167H_rep1 >= 5 & SUB$R167H_rep2 >= 5,]
SUB$DiffBind <- WTvsRtoC$DiffBind
SUB$DiffBind <- factor(SUB$DiffBind, levels = c("N.S."))
aa <- (t(table(SUB[,c("feature","DiffBind")])) / rowSums(t(table(SUB[,c("feature","DiffBind")]))))*100
barplot(aa, las=2, col=greys[6], ylim=c(0,180))



SUB <- IGF2BP1s
SUB$feature <- gsub("exon","ncRNA",SUB$feature)
SUB <- SUB[SUB$wt_rep1 >= 5 & SUB$wt_rep2 >= 5 & SUB$R167C_rep1 >= 5 & SUB$R167C_rep2 >= 5 & SUB$R167H_rep1 >= 5 & SUB$R167H_rep2 >= 5,]
SUB$DiffBind <- WTvsRtoH$DiffBind
SUB$DiffBind <- factor(SUB$DiffBind, levels = c("wt","R167H","N.S."))
barplot((t(table(SUB[,c("feature","DiffBind")])) / rowSums(t(table(SUB[,c("feature","DiffBind")]))))*100, las=2, col=c(greys[2],yellows[2], greys[6]), ylim=c(0,180))



SUB <- IGF2BP2s
SUB$feature <- gsub("exon","ncRNA",SUB$feature)
SUB <- SUB[SUB$A_rep1 >= 5 & SUB$A_rep2 >= 5 & SUB$B_rep1 >= 5 & SUB$B_rep2 >= 5 ,]
SUB$DiffBind <- AvsB$DiffBind
SUB$DiffBind <- factor(SUB$DiffBind, levels = c("Isoform_A","Isoform_B","N.S."))
barplot((t(table(SUB[,c("feature","DiffBind")])) / rowSums(t(table(SUB[,c("feature","DiffBind")]))))*100, las=2, col=c(greys[2],oranges[3], greys[6]), ylim=c(0,180))



SUB <- IGF2BP3s
SUB$feature <- gsub("exon","ncRNA",SUB$feature)
SUB <- SUB[SUB$wt_rep1 >= 5 & SUB$wt_rep2 >= 5 & SUB$R525C_rep1 >= 5 & SUB$R525C_rep2 >= 5 & SUB$I474M_rep1 >= 5 & SUB$I474M_rep2 >= 5,]
SUB$DiffBind <- WTvsI474M$DiffBind
SUB$DiffBind <- factor(SUB$DiffBind, levels = c("wt","I474M","N.S."))
barplot((t(table(SUB[,c("feature","DiffBind")])) / rowSums(t(table(SUB[,c("feature","DiffBind")]))))*100, las=2, col=c(greys[2],pinks[2], greys[6]), ylim=c(0,180))


SUB <- IGF2BP3s
SUB$feature <- gsub("exon","ncRNA",SUB$feature)
SUB <- SUB[SUB$wt_rep1 >= 5 & SUB$wt_rep2 >= 5 & SUB$R525C_rep1 >= 5 & SUB$R525C_rep2 >= 5 & SUB$I474M_rep1 >= 5 & SUB$I474M_rep2 >= 5,]
SUB$DiffBind <- WTvsR525C$DiffBind
SUB$DiffBind <- factor(SUB$DiffBind, levels = c("wt","R525C","N.S."))
barplot((t(table(SUB[,c("feature","DiffBind")])) / rowSums(t(table(SUB[,c("feature","DiffBind")]))))*100, las=2, col=c(greys[2],acqua_blues[5], greys[6]), ylim=c(0,180))



#dev.off()






```






Pathway Diff Bind for BP1s

```{r, fig.width=15, fig.height=20}

Inputs <- read.table(file="/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/INPUTs.TPMs.txt", stringsAsFactors = F, header=T)

SUB <- IGF2BP1s
SUB <- SUB[SUB$wt_rep1 >= 5 & SUB$wt_rep2 >= 5 & SUB$R167C_rep1 >= 5 & SUB$R167C_rep2 >= 5 & SUB$R167H_rep1 >= 5 & SUB$R167H_rep2 >= 5,]

hits <- unique(SUB$gene_ID)
gene_IDs_hits <- unique(unlist(str_split(unique(hits),"\\."))[seq(1,length(unique(hits))*2,2)])
Input_genes <- unique(unlist(str_split(unique(Inputs$gene_ID),"\\."))[seq(1,length(unique(Inputs$gene_ID))*2,2)])
Input_genes <- c(Input_genes, gene_IDs_hits)
IDs <- mapIds(x=org.Hs.eg.db, keys=gene_IDs_hits, column='ENTREZID', keytype = 'ENSEMBL')
IDs <- as.character(IDs)[is.na(as.character(IDs)) == FALSE]
Reactomes <- enrichPathway(
                gene = IDs,
                organism      = "human",
                pAdjustMethod = "BH",
                readable      = TRUE
                )
Reactomes_BP1s <- as.data.frame(Reactomes)
AA <- Reactomes_BP1s[,c("Description","geneID")]
BB <- data.frame(c("NO","NO"), c("NO","NO"))
colnames(BB) <- c("Pathway","gene_name")
for (i in 1:nrow(AA)){
  
  gene_name <- unique(unlist(str_split(AA[i,"geneID"],"\\/")))
  Pathway <- rep(AA[i,"Description"],length(gene_name))
  CC <- data.frame(Pathway,gene_name)
  BB <- rbind(BB,CC)
  
}
BB <- BB[-c(1:2),]
BB <- merge(BB,SUB[,c("gene_name","wt_rep1","wt_rep2","R167C_rep1","R167C_rep2","R167H_rep1","R167H_rep2")], by="gene_name")
BB <- BB[,-1]
BB <- melt(BB)
colnames(BB) <- c("Pathway","RBP","NormCounts")
BB <- BB %>% group_by(Pathway,RBP) %>% summarise(PathCounts = sum(NormCounts))
BB <- dcast(BB, formula = Pathway~RBP)


par(mfrow=c(1,2))
vsRtoH <- BB[,c(2,3,4,5)] 
x <- vsRtoH
group <- c(1,1,2,2)
y <- DGEList(counts=x, group = group)
design <- model.matrix(~group)
y <- estimateDisp(y,design)
#y$samples$lib.size <- c(1000000,1000000,1000000,1000000)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit)
X <-qlf$table
X$p.adjust <- -log10(p.adjust(X$PValue, method = "BH"))
X$PValue <- -log10(X$PValue)
BB$FDR <- X$p.adjust
BB$Significant <- BB$FDR >= 1.30103
BB$log2FC <- X$logFC
BB$log2CPM <- X$logCPM
BB$wt_won <- BB$log2FC <= -0 & BB$FDR >= 1.30103
BB$RtoC_won <- BB$log2FC >= 0 & BB$FDR >= 1.30103
BB$DiffBind <- paste(BB$wt_won,BB$RtoC_won,sep = "_")
BB$DiffBind <- gsub("FALSE_FALSE","N.S.",BB$DiffBind)
BB$DiffBind <- gsub("FALSE_TRUE","R167C",BB$DiffBind)
BB$DiffBind <- gsub("TRUE_FALSE","wt",BB$DiffBind)
BB$COLOR <- gsub("N.S.",greys[6],BB$DiffBind)
BB$COLOR <- gsub("wt","LoF",BB$COLOR)
BB$COLOR <- gsub("R167C","GoF",BB$COLOR)
YY <- BB
YY$Comparison <- rep("wt_R167C", nrow(YY)) 

vsRtoH <- BB[,c(2,3,6,7)] 
x <- vsRtoH
group <- c(1,1,2,2)
y <- DGEList(counts=x, group = group)
design <- model.matrix(~group)
y <- estimateDisp(y,design)
#y$samples$lib.size <- c(1000000,1000000,1000000,1000000)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit)
X <-qlf$table
X$p.adjust <- -log10(p.adjust(X$PValue, method = "BH"))
X$PValue <- -log10(X$PValue)
BB$FDR <- X$p.adjust
BB$Significant <- BB$FDR >= 1.30103
BB$log2FC <- X$logFC
BB$log2CPM <- X$logCPM
BB$wt_won <- BB$log2FC <= -0 & BB$FDR >= 1.30103
BB$RtoC_won <- BB$log2FC >= 0 & BB$FDR >= 1.30103
BB$DiffBind <- paste(BB$wt_won,BB$RtoC_won,sep = "_")
BB$DiffBind <- gsub("FALSE_FALSE","N.S.",BB$DiffBind)
BB$DiffBind <- gsub("FALSE_TRUE","R167H",BB$DiffBind)
BB$DiffBind <- gsub("TRUE_FALSE","wt",BB$DiffBind)
BB$COLOR <- gsub("N.S.",greys[6],BB$DiffBind)
BB$COLOR <- gsub("wt","LoF",BB$COLOR)
BB$COLOR <- gsub("R167H","GoF", BB$COLOR)
YYY <- BB
YYY$Comparison <- rep("wt_R167H", nrow(YYY)) 


#### Pathway Diff Bind for BP2s



Inputs <- read.table(file="/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/INPUTs.TPMs.txt", stringsAsFactors = F, header=T)


SUB <- IGF2BP2s
SUB <- SUB[SUB$A_rep1 >= 5 & SUB$A_rep2 >= 5 & SUB$B_rep1 >= 5 & SUB$B_rep2 >= 5 ,]

hits <- unique(SUB$gene_ID)
gene_IDs_hits <- unique(unlist(str_split(unique(hits),"\\."))[seq(1,length(unique(hits))*2,2)])
Input_genes <- unique(unlist(str_split(unique(Inputs$gene_ID),"\\."))[seq(1,length(unique(Inputs$gene_ID))*2,2)])
Input_genes <- c(Input_genes, gene_IDs_hits)
IDs <- mapIds(x=org.Hs.eg.db, keys=gene_IDs_hits, column='ENTREZID', keytype = 'ENSEMBL')
IDs <- as.character(IDs)[is.na(as.character(IDs)) == FALSE]
Reactomes <- enrichPathway(
                gene = IDs,
                organism      = "human",
                pAdjustMethod = "BH",
                readable      = TRUE
                )
Reactomes_BP1s <- as.data.frame(Reactomes)
AA <- Reactomes_BP1s[,c("Description","geneID")]
BB <- data.frame(c("NO","NO"), c("NO","NO"))
colnames(BB) <- c("Pathway","gene_name")
for (i in 1:nrow(AA)){
  
  gene_name <- unique(unlist(str_split(AA[i,"geneID"],"\\/")))
  Pathway <- rep(AA[i,"Description"],length(gene_name))
  CC <- data.frame(Pathway,gene_name)
  BB <- rbind(BB,CC)
  
}
BB <- BB[-c(1:2),]
BB <- merge(BB,SUB[,c("gene_name","A_rep1","A_rep2","B_rep1","B_rep2")], by="gene_name")
BB <- BB[,-1]
BB <- melt(BB)
colnames(BB) <- c("Pathway","RBP","NormCounts")
BB <- BB %>% group_by(Pathway,RBP) %>% summarise(PathCounts = sum(NormCounts))
BB <- dcast(BB, formula = Pathway~RBP)




AvsB <- BB[,c(2,3,4,5)]
x <- AvsB
group <- c(1,1,2,2)
y <- DGEList(counts=x, group = group)
design <- model.matrix(~group)
y <- estimateDisp(y,design)
#y$samples$lib.size <- c(1000000,1000000,1000000,1000000)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit)
X <-qlf$table
X$p.adjust <- -log10(p.adjust(X$PValue, method = "BH"))
X$PValue <- -log10(X$PValue)
BB$FDR <- X$p.adjust
BB$Significant <- BB$FDR >= 1.30103
BB$log2FC <- X$logFC
BB$log2CPM <- X$logCPM
BB$wt_won <- BB$log2FC <= -0 & BB$FDR >= 1.30103
BB$RtoC_won <- BB$log2FC >= 0 & BB$FDR >= 1.30103
BB$DiffBind <- paste(BB$wt_won,BB$RtoC_won,sep = "_")
BB$DiffBind <- gsub("FALSE_FALSE","N.S.",BB$DiffBind)
BB$DiffBind <- gsub("FALSE_TRUE","Isoform_B",BB$DiffBind)
BB$DiffBind <- gsub("TRUE_FALSE","Isoform_A",BB$DiffBind)
BB$COLOR <- gsub("N.S.",greys[6],BB$DiffBind)
BB$COLOR <- gsub("Isoform_A","LoF",BB$COLOR)
BB$COLOR <- gsub("Isoform_B","GoF",BB$COLOR)
BBB <- BB
BBB$Comparison <- rep("A_vs_B", nrow(BBB)) 




####### Pathway Diff Bind for BP2s



Inputs <- read.table(file="/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/INPUTs.TPMs.txt", stringsAsFactors = F, header=T)


SUB <- IGF2BP3s
SUB <- SUB[SUB$wt_rep1 >= 5 & SUB$wt_rep2 >= 5 & SUB$R525C_rep1 >= 5 & SUB$R525C_rep2 >= 5 & SUB$I474M_rep1 >= 5 & SUB$I474M_rep2 >= 5,]


hits <- unique(SUB$gene_ID)
gene_IDs_hits <- unique(unlist(str_split(unique(hits),"\\."))[seq(1,length(unique(hits))*2,2)])
Input_genes <- unique(unlist(str_split(unique(Inputs$gene_ID),"\\."))[seq(1,length(unique(Inputs$gene_ID))*2,2)])
Input_genes <- c(Input_genes, gene_IDs_hits)
IDs <- mapIds(x=org.Hs.eg.db, keys=gene_IDs_hits, column='ENTREZID', keytype = 'ENSEMBL')
IDs <- as.character(IDs)[is.na(as.character(IDs)) == FALSE]
Reactomes <- enrichPathway(
                gene = IDs,
                organism      = "human",
                pAdjustMethod = "BH",
                readable      = TRUE
                )
Reactomes_BP1s <- as.data.frame(Reactomes)
AA <- Reactomes_BP1s[,c("Description","geneID")]
BB <- data.frame(c("NO","NO"), c("NO","NO"))
colnames(BB) <- c("Pathway","gene_name")
for (i in 1:nrow(AA)){
  
  gene_name <- unique(unlist(str_split(AA[i,"geneID"],"\\/")))
  Pathway <- rep(AA[i,"Description"],length(gene_name))
  CC <- data.frame(Pathway,gene_name)
  BB <- rbind(BB,CC)
  
}
BB <- BB[-c(1:2),]
BB <- merge(BB,SUB[,c("gene_name","wt_rep1","wt_rep2","I474M_rep1","I474M_rep2","R525C_rep1","R525C_rep2")], by="gene_name")
BB <- BB[,-1]
BB <- melt(BB)
colnames(BB) <- c("Pathway","RBP","NormCounts")
BB <- BB %>% group_by(Pathway,RBP) %>% summarise(PathCounts = sum(NormCounts))
BB <- dcast(BB, formula = Pathway~RBP)




BP3<- BB[,c(2,3,4,5)]
x <- BP3
group <- c(1,1,2,2)
y <- DGEList(counts=x, group = group)
design <- model.matrix(~group)
y <- estimateDisp(y,design)
#y$samples$lib.size <- c(1000000,1000000,1000000,1000000)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit)
X <-qlf$table
X$p.adjust <- -log10(p.adjust(X$PValue, method = "BH"))
X$PValue <- -log10(X$PValue)
BB$FDR <- X$p.adjust
BB$Significant <- BB$FDR >= 1.30103
BB$log2FC <- X$logFC
BB$log2CPM <- X$logCPM
BB$wt_won <- BB$log2FC <= -0 & BB$FDR >= 1.30103
BB$RtoC_won <- BB$log2FC >= 0 & BB$FDR >= 1.30103
BB$DiffBind <- paste(BB$wt_won,BB$RtoC_won,sep = "_")
BB$DiffBind <- gsub("FALSE_FALSE","N.S.",BB$DiffBind)
BB$DiffBind <- gsub("FALSE_TRUE","I474M",BB$DiffBind)
BB$DiffBind <- gsub("TRUE_FALSE","wt",BB$DiffBind)
BB$COLOR <- gsub("N.S.",greys[6],BB$DiffBind)
BB$COLOR <- gsub("wt","LoF",BB$COLOR)
BB$COLOR <- gsub("I474M","GoF",BB$COLOR)
DD <- BB
DD$Comparison <- rep("wt_I474M", nrow(DD)) 



BP3<- BB[,c(2,3,6,7)]
x <- BP3
group <- c(1,1,2,2)
y <- DGEList(counts=x, group = group)
design <- model.matrix(~group)
y <- estimateDisp(y,design)
#y$samples$lib.size <- c(1000000,1000000,1000000,1000000)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit)
X <-qlf$table
X$p.adjust <- -log10(p.adjust(X$PValue, method = "BH"))
X$PValue <- -log10(X$PValue)
BB$FDR <- X$p.adjust
BB$Significant <- BB$FDR >= 1.30103
BB$log2FC <- X$logFC
BB$log2CPM <- X$logCPM
BB$wt_won <- BB$log2FC <= -0 & BB$FDR >= 1.30103
BB$RtoC_won <- BB$log2FC >= 0 & BB$FDR >= 1.30103
BB$DiffBind <- paste(BB$wt_won,BB$RtoC_won,sep = "_")
BB$DiffBind <- gsub("FALSE_FALSE","N.S.",BB$DiffBind)
BB$DiffBind <- gsub("FALSE_TRUE","R525C",BB$DiffBind)
BB$DiffBind <- gsub("TRUE_FALSE","wt",BB$DiffBind)
BB$COLOR <- gsub("N.S.",greys[6],BB$DiffBind)
BB$COLOR <- gsub("wt","LoF",BB$COLOR)
BB$COLOR <- gsub("R525C","GoF",BB$COLOR)
DDD <- BB
DDD$Comparison <- rep("wt_R525C", nrow(DDD)) 



YY <- YY[,c("Pathway","FDR","log2FC","COLOR","Comparison")]
YYY <- YYY[,c("Pathway","FDR","log2FC","COLOR","Comparison")]
BBB <- BBB[,c("Pathway","FDR","log2FC","COLOR","Comparison")]
DD <- DD[,c("Pathway","FDR","log2FC","COLOR","Comparison")]
DDD <- DDD[,c("Pathway","FDR","log2FC","COLOR","Comparison")]


#YY <- YY[YY$FDR >= 1.30103,]
#YY1 <- YY[order(YY$log2FC),][1:5,]
#YY2 <- YY[order(-YY$log2FC),][1:5,]
#YY <- rbind(YY1,YY2)
#YY$Pathway <- factor(YY$Pathway, levels=YY[order(YY$log2FC),"Pathway"])

YYY <- YYY[YYY$FDR >= 1.30103,]
YYY1 <- YYY[order(YYY$log2FC),][1:5,]
YYY2 <- YYY[order(-YYY$log2FC),][1:5,]
YYY <- rbind(YYY1,YYY2)
YYY$Pathway <- factor(YYY$Pathway, levels=YYY[order(YYY$log2FC),"Pathway"])

BBB <- BBB[BBB$FDR >= 1.30103,]
BBB1 <- BBB[order(BBB$log2FC),][1:5,]
BBB2 <- BBB[order(-BBB$log2FC),][1:5,]
BBB <- rbind(BBB1,BBB2)
BBB$Pathway <- factor(BBB$Pathway, levels=BBB[order(BBB$log2FC),"Pathway"])

DD <- DD[DD$FDR >= 1.30103,]
DD1 <- DD[order(DD$log2FC),][1:5,]
DD2 <- DD[order(-DD$log2FC),][1:5,]
DD <- rbind(DD1,DD2)
DD$Pathway <- factor(DD$Pathway, levels=DD[order(DD$log2FC),"Pathway"])

DDD <- DDD[DDD$FDR >= 1.30103,]
DDD1 <- DDD[order(DDD$log2FC),][1:5,]
DDD2 <- DDD[order(-DDD$log2FC),][1:5,]
DDD <- rbind(DDD1,DDD2)
DDD$Pathway <- factor(DDD$Pathway, levels=DDD[order(DDD$log2FC),"Pathway"])


ALL <- rbind(YYY,BBB,DD,DDD)
ALL$Pathway <- as.character(ALL$Pathway)
ALL$Pathway <- factor(ALL$Pathway, levels=unique(ALL[order(ALL$Comparison),"Pathway"]))

aa <- ggplot(data=YYY, aes(x=Pathway, y=log2FC, alpha=-(10^-FDR), fill=COLOR)) +
  geom_bar(stat="identity") +
  coord_flip() +
  scale_alpha_binned(range = c(0.5,1)) +
  scale_fill_manual(values = c(yellows[3],greys[2])) +
  ylim(-0.45,0.45) +
  facet_wrap(~Comparison) + 
  theme_pubr(legend = "left")

bb <- ggplot(data=BBB, aes(x=Pathway, y=log2FC, alpha=-(10^-FDR), fill=COLOR)) +
  geom_bar(stat="identity") +
  coord_flip() +
  scale_alpha_binned(range = c(0.5,1)) +
  scale_fill_manual(values = c(oranges[3],greys[2])) +
  ylim(-0.45,0.45) +
  facet_wrap(~Comparison) + 
  theme_pubr(legend = "left")

cc <- ggplot(data=DD, aes(x=Pathway, y=log2FC, alpha=-(10^-FDR), fill=COLOR)) +
  geom_bar(stat="identity") +
  coord_flip() +
  scale_alpha_binned(range = c(0.5,1)) +
  scale_fill_manual(values = c(pinks[2],greys[2])) +
  ylim(-0.45,0.45) +
  facet_wrap(~Comparison) + 
  theme_pubr(legend = "left")

dd <- ggplot(data=DDD, aes(x=Pathway, y=log2FC, alpha=-(10^-FDR), fill=COLOR)) +
  geom_bar(stat="identity") +
  coord_flip() +
  scale_alpha_binned(range = c(0.5,1)) +
  scale_fill_manual(values = c(acqua_blues[5],greys[2])) +
  ylim(-0.45,0.45) +
  facet_wrap(~Comparison) + 
  theme_pubr(legend = "left")

#pdf(file = "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_7/Pathways_barplots_GoFs_n_LoFs.pdf",15,20)
ggarrange(aa,bb,cc,dd, nrow=4) 
#dev.off()

```


kmers k=5

```{r, fig.width=15, fig.height=15}

PATH <- "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/RAP_Annotated_PEAKS/IGF2BPs/"
PEAKS <- list.files(path = PATH)
Names <- unlist(str_split(PEAKS,"\\.peak"))[grep("txt",unlist(str_split(PEAKS,"\\.peak")),invert = T)]
Tables <- list()
for (i in 1:length(Names)){
  Tables[[i]] <- as.data.frame(read.table(paste(PATH,PEAKS[i],sep = ""), stringsAsFactors = F, header = T))
  Tables[[i]]$positive_fa_sub <- str_sub(Tables[[i]]$positive_fa,50,150)
}
names(Tables) <- Names
Names


kmer_table <- as.data.frame(Names)
STRINGS <- c()
for(i in Names){
  STRINGS <- c(STRINGS,as.character(paste(Tables[[i]]$positive_fa_sub,collapse = 'NN')))
}
k=5
kmer_table <- cbind(kmer_table,STRINGS)
kmer_table[] <- lapply(kmer_table, as.character)
bases <- c("A","C","T","G")
kmers <- unite(as.data.frame(permutations(n=4,r=k,v=bases,repeats.allowed = T)),col = kmers,sep = "")
for (i in 1:nrow(kmer_table)){
kmers<-cbind(kmers,str_count(kmer_table[i,2],kmers[,1])/length(kmer_table[i,2]))
}
rownames(kmers) <- kmers$kmers
pos_kmers_5 <- as.data.frame(kmers[,-1], row.names = rownames(kmers))
colnames(pos_kmers_5) <- kmer_table$Names









PATH <- "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/RAP_Annotated_PEAKS/P102/"
PEAKS <- list.files(path = PATH)
Names <- unlist(str_split(PEAKS,"\\.peak"))[grep("txt",unlist(str_split(PEAKS,"\\.peak")),invert = T)]
Tables <- list()
for (i in 1:length(Names)){
  Tables[[i]] <- as.data.frame(read.table(paste(PATH,PEAKS[i],sep = ""), stringsAsFactors = F, header = T))
  Tables[[i]]$negative_fa_sub <- str_sub(Tables[[i]]$positive_fa,50,150)
}
names(Tables) <- Names
Names

kmer_table <- as.data.frame(Names)
STRINGS <- c()
for(i in Names){
  STRINGS <- c(STRINGS,as.character(paste(Tables[[i]]$negative_fa_sub,collapse = 'NN')))
}
k=5
kmer_table <- cbind(kmer_table,STRINGS)
kmer_table[] <- lapply(kmer_table, as.character)
bases <- c("A","C","T","G")
kmers <- unite(as.data.frame(permutations(n=4,r=k,v=bases,repeats.allowed = T)),col = kmers,sep = "")
for (i in 1:nrow(kmer_table)){
kmers<-cbind(kmers,str_count(kmer_table[i,2],kmers[,1])/length(kmer_table[i,2]))
}
rownames(kmers) <- kmers$kmers
neg_kmers_5 <- as.data.frame(kmers[,-1], row.names = rownames(kmers))
colnames(neg_kmers_5) <- kmer_table$Names




Z_k5 <- as.data.frame( scale(t((t(pos_kmers_5)/colSums(pos_kmers_5)))/rowMeans( t(t(neg_kmers_5)/colSums(neg_kmers_5)))) )


Mots_YTH <- Z_k5[grep("GGACT|TGGAC|GACTG|GACTC|CGGAC", rownames(Z_k5)),]
Mots_YTH$COLOR <- rep(pinks[2],nrow(Mots_YTH))
Mots_YTH$SIZE <- rep(1.5,nrow(Mots_YTH))
Mots_IGFs <- Z_k5[grep("ACAAC|CAAAC|AACAC|CACAA|CAACA", rownames(Z_k5)),]
Mots_IGFs$COLOR <- rep(acqua_blues[6],nrow(Mots_IGFs))
Mots_IGFs$SIZE <- rep(1.5,nrow(Mots_IGFs))
Mots_RBF <- Z_k5[grep("GCATG|TGCAT|GAATG|GCACG", rownames(Z_k5)),]
Mots_RBF$COLOR <- rep(yellows[2],nrow(Mots_RBF))
Mots_RBF$SIZE <- rep(1.5,nrow(Mots_RBF))
AA <- Z_k5[grep("CGGAC|GGACT|TGGAC|GACTG|GACTC|ACAAC|CAAAC|AACAC|CACAA|CAACA|GCATG|TGCAT|GAATG|GCACG", rownames(Z_k5), invert=T),]
AA$COLOR <- rep(alpha(greys[6],0.25),nrow(AA))
AA$SIZE <- rep(1,nrow(AA))
Z_k5_bis <- rbind(AA, Mots_IGFs, Mots_YTH, Mots_RBF)
# Correlation panel
panel.cor <- function(x, y){
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    rho <- round(cor(x, y, method = "spearman"), digits=2)
    txt <- as.character(rho)
    cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = 3)
}
# Customize upper panel
upper.panel<-function(x, y){
  points(x,y, pch = 19, col = Z_k5_bis$COLOR, cex = Z_k5_bis$SIZE)
}

#pdf(file = "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_7/Kmers_Scatterplots.pdf",12.5,12.5)
pairs(Z_k5_bis[,c(9,1:8,10)], lower.panel = panel.cor, upper.panel = upper.panel)
#dev.off()

```

```{r}

IGF2BP2s <- read.table(file = "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_7/PEAKS/IGF2BP2s.metafile.txt", header = T, stringsAsFactors = F, sep = "\t")



SUB <- IGF2BP2s[,7:10] 
SUB <- SUB[SUB$A_rep1 >= 5 & SUB$A_rep2 >= 5 & SUB$B_rep1 >= 5 & SUB$B_rep2 >= 5 ,]
nrow(SUB)

AvsB <- SUB
x <- AvsB
group <- c(1,1,2,2)
y <- DGEList(counts=x, group = group)
design <- model.matrix(~group)
y <- estimateDisp(y,design)
y$samples$lib.size <- c(1000000,1000000,1000000,1000000)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit)
X <-qlf$table
X$p.adjust <- -log10(p.adjust(X$PValue, method = "BH"))
X$PValue <- -log10(X$PValue)
AvsB$FDR <- X$p.adjust
AvsB$Significant <- AvsB$FDR >= 1.30103
AvsB$log2FC <- X$logFC
AvsB$log2CPM <- X$logCPM
AvsB$wt_won <- AvsB$log2FC <= -1 & AvsB$FDR >= 1.30103
AvsB$RtoC_won <- AvsB$log2FC >= 1 & AvsB$FDR >= 1.30103
AvsB$DiffBind <- paste(AvsB$wt_won,AvsB$RtoC_won,sep = "_")
AvsB$DiffBind <- gsub("FALSE_FALSE","N.S.",AvsB$DiffBind)
AvsB$DiffBind <- gsub("FALSE_TRUE","Isoform_B",AvsB$DiffBind)
AvsB$DiffBind <- gsub("TRUE_FALSE","Isoform_A",AvsB$DiffBind)
AvsB$comparison <- rep("Isof_A vs Isof_B", nrow(AvsB))
COMPARISONS_2 <- AvsB[,c(5,7,8,11,12)]
colnames(COMPARISONS_2) <- c("minuslog10FDR", "log2FC", "log2CPM", "DiffBind", "comparison")
  


SUB <- IGF2BP2s
SUB$feature <- gsub("exon","ncRNA",SUB$feature)
SUB <- SUB[SUB$A_rep1 >= 5 & SUB$A_rep2 >= 5 & SUB$B_rep1 >= 5 & SUB$B_rep2 >= 5 ,]
SUB$DiffBind <- AvsB$DiffBind
SUB$minuslog10FDR <- AvsB$FDR
SUB$log2FC <- AvsB$log2FC
SUB <- SUB[SUB$DiffBind != "N.S.",]
SUB$DiffBind <- factor(SUB$DiffBind, levels = c("Isoform_A","Isoform_B"))
SUB <- SUB[SUB$gene_type != "protein_coding",]
aa <- t(table(SUB[,c("gene_type","DiffBind")]))

#pdf(file = "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_7/Barplot_ncRNAs_diff_bound.pdf",5,5)
barplot(aa[,order(-aa[2,])], las=2, col=c(greys[2],oranges[3]))
#dev.off()

SUB[,c("DiffBind","minuslog10FDR","log2FC","gene_type","gene_name","chr","Summit_start","Summit_end","strand")]

```




```{r}
library(stringr)
library(ggplot2)
library(factoextra)
library(tidyverse)
library(VennDiagram)

library(ggfortify)
library(eulerr)
library(dplyr)
library(tidyr)
library(ggplot2)

library(ggpubr)
```
# Import tables
```{r}
HOM_list <- read.table(file = "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_7/SOFIA/HOM_MouseHuman.txt", header=T)
ensembl <- 
IGFs_homo <- HOM_list[HOM_list$Symbol %in% c("IGF2BP1", "IGF2BP2", "IGF2BP3"),8]
IGFs_mouse  <- HOM_list[HOM_list$Symbol %in% c("IGF2BP1", "IGF2BP2", "IGF2BP3"),9]
```

# Import and format the data
## Human development
```{r}
# File with the RPKM values from the Cardaso et al study downloaded from https://apps.kaessmannlab.org/evodevoapp/
human_rpkm <- read.table(file = "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_7/SOFIA/Human_rpkm.txt", header = T, row.names = 1)
human_rpkm <- human_rpkm[row.names(human_rpkm) %in% IGFs_homo,]
Hs_Liver <- human_rpkm[,str_starts(colnames(human_rpkm), "Liver") == TRUE]
gene_names <- row.names(Hs_Liver)
Hs_Liver <- as.data.frame(lapply(Hs_Liver, function(x) as.numeric(as.character(x))))
row.names(Hs_Liver) <- gene_names

# Mean the RPKM of each technical replicate
Hs_Liver <- as.data.frame(t(Hs_Liver)) 
group <- sapply(str_split(as.character(row.names(Hs_Liver)), "\\."), `[`, 2)
Hs_Liver <- cbind(group, Hs_Liver)
Hs_Liver <- Hs_Liver %>% group_by(group) %>% summarise_all(funs(as.numeric(mean(., na.rm=TRUE))))

## Birth = timepoint between 13 and 14
human_stages <- read.table("/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_7/SOFIA/Human_stages_key.txt")
names(human_stages) <- c("group", "stage")
Hs_Liver <- merge(human_stages, Hs_Liver, by="group")
row.names(Hs_Liver) <- paste(Hs_Liver$stage, Hs_Liver$group, sep="_")
Hs_Liver <- Hs_Liver[,-c(1,2)]
Hs_Liver <- as.data.frame(t(Hs_Liver))
Hs_Liver$gene <- row.names(Hs_Liver)

plotting_Hs <- tidyr::gather(Hs_Liver, "Developmental_Stage", "expression", -c(gene))
plotting_Hs $Stage <- sapply(str_split(as.character(plotting_Hs $Developmental_Stage), "_"), `[`, 1)
plotting_Hs $Developmental_Stage <- sapply(str_split(as.character(plotting_Hs $Developmental_Stage), "_"), `[`, 2)
plotting_Hs$Study <- rep("Cardosso", nrow(plotting_Hs))
plotting_Hs$species <- rep("H.sapiens", nrow(plotting_Hs))
```

## Mouse development
```{r}
# File with the RPKM values from the Cardaso et al study downloaded from https://apps.kaessmannlab.org/evodevoapp/
mouse_rpkm <- read.table("/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_7/SOFIA/Mouse_rpkm.txt", header=T, row.names=1)
mouse_rpkm <- mouse_rpkm[row.names(mouse_rpkm) %in% IGFs_mouse,]
Mm_Liver <- mouse_rpkm[,str_starts(colnames(mouse_rpkm), "Liver") == TRUE]
gene_names <- row.names(Mm_Liver)
Mm_Liver <- as.data.frame(lapply(Mm_Liver, function(x) as.numeric(as.character(x))))
row.names(Mm_Liver) <- gene_names

# Mean the RPKM of each technical replicate
Mm_Liver <- as.data.frame(t(Mm_Liver)) 
group <- sapply(str_split(as.character(row.names(Mm_Liver)), "\\."), `[`, 2)
Mm_Liver <- cbind(group, Mm_Liver)
Mm_Liver <- Mm_Liver %>% group_by(group) %>% summarise_all(funs(as.numeric(mean(., na.rm=TRUE))))

mouse_stages <- read.table("/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_7/SOFIA/Mouse_stages_key.txt")
names(mouse_stages) <- c("group", "stage")
Mm_Liver <- merge(mouse_stages, Mm_Liver, by="group")
row.names(Mm_Liver) <- paste(Mm_Liver$stage, Mm_Liver$group, sep="_")
Mm_Liver <- Mm_Liver[,-c(1,2)]
Mm_Liver <- as.data.frame(t(Mm_Liver))
Mm_Liver$gene <- row.names(Mm_Liver)
plotting_Mm <- tidyr::gather(Mm_Liver, "Developmental_Stage", "expression", -c(gene))
plotting_Mm $Stage <- sapply(str_split(as.character(plotting_Mm$Developmental_Stage), "_"), `[`, 1)
plotting_Mm $Developmental_Stage <- sapply(str_split(as.character(plotting_Mm$Developmental_Stage), "_"), `[`, 2)
plotting_Mm$Study <- rep("Cardosso", nrow(plotting_Mm))
plotting_Mm$species <- rep("M.musculus", nrow(plotting_Mm))
```

## Human cancer cell lines
```{r}
# Downoladed table of normalized counts from the Schmitt et al study PLOS ONE
human_cell_lines <- read.table("/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_7/SOFIA/journal.pgen.1006024.s014.TSV", header = T)
human_cell_lines <- human_cell_lines[,c(1,3:5)]
human_cell_lines <- gather(human_cell_lines, "Sample", "Normalized_count", 2:4)
human_cell_lines <- human_cell_lines[human_cell_lines$Gene %in% IGFs_homo,]
```

## Mouse cancer cell lines
```{r}
mouse_cell_lines <- read.table("/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_7/SOFIA/journal.pgen.1006024.s015.TSV", header = T)
mouse_cell_lines <- mouse_cell_lines[,c(1,4:6)]
mouse_cell_lines <- gather(mouse_cell_lines, "Sample", "Normalized_count", 2:4)
mouse_cell_lines <- mouse_cell_lines[mouse_cell_lines$Gene %in% IGFs_mouse,]
```

##### Plotting expression during Development

```{r}
making_plots_3_suppl <- function(gene_list){
  colors <- c(acqua_blues[c(2,6,9)])
  genes <- gene_list
  mouse_gene <- HOM_list[HOM_list$homo_ensembl %in% genes, "mouse_ensembl"]
  data_hs <- merge(plotting_Hs[plotting_Hs$gene %in% genes,], HOM_list[,c(7,8)], by.x=1, by.y=2)
  hs_data <- arrange(data_hs, Symbol)
  hs_data$Developmental_Stage <- factor(hs_data$Developmental_Stage, levels = unique(arrange(hs_data, as.numeric(Stage))[,2]))
  data_mm <- merge(plotting_Mm[plotting_Mm$gene %in% mouse_gene,], HOM_list[,c(4,9)], by.x=1, by.y=2)
  mm_data <- arrange(data_mm, mouse_gene)
  mm_data$Developmental_Stage <- factor(mm_data$Developmental_Stage, levels = unique(arrange(mm_data, as.numeric(Stage))[,2]))

    
  Hs_dev <- ggplot(data=hs_data) +
    geom_line(aes(x=Developmental_Stage, y=(as.numeric(expression)), group=Symbol, col=Symbol), size=2, alpha=0.75) +
    geom_point(aes(x=Developmental_Stage, y=(as.numeric(expression)), fill=Symbol), size=5, shape=21) +
    geom_vline(xintercept = 13.5, linetype = 2, color="gray", size=0.8) +
    scale_color_manual(values=colors) +
    xlab("Developmental Stage") +
    ylab("RPKM") +
    ggtitle("Homo_sapiens") +
    ylim(0,80) +
    theme_pubr() +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors)
  Hs_dev <- Hs_dev + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position = "bottom")
      
  Mm_dev <- ggplot(data=mm_data) +
    geom_line(aes(x=Developmental_Stage, y=(as.numeric(expression)), group=mouse_gene, col=mouse_gene), size=2, alpha=0.75) +
    geom_point(aes(x=Developmental_Stage, y=(as.numeric(expression)), fill=mouse_gene), size=5, shape=21) +
    geom_vline(xintercept = 9.5, linetype=2, color="gray", size=0.8) +
    scale_color_manual(values=colors) +
    xlab("Developmental Stage") +
    ylab("RPKM") +
    ggtitle("Mus musculus") +
    ylim(0,80) +
    theme_pubr() + 
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors)
  Mm_dev <- Mm_dev + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "bottom")
  
  gene_plots <- ggarrange(Hs_dev, Mm_dev, ncol=2, align = "hv")
  
  
  
  return(gene_plots)
}
```

```{r, fig.height=8,fig.width=8}
#pdf(file = "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_7/IGF2BPs_development.pdf",8,8)
making_plots_3_suppl(IGFs_homo)
#dev.off()
```


```{r}
making_plots_suppl <- function(x){
  colors <- c(acqua_blues[c(2,6,9)])
  human_cell_lines <- merge(human_cell_lines, HOM_list[,c(7,8)], by.x=1, by.y=2)
  mouse_cell_lines <- merge(mouse_cell_lines, HOM_list[,c(4,9)], by.x=1, by.y=2)
  
  cancer_Hs <-ggplot(human_cell_lines)+
    geom_bar(aes(x=Sample, y=Normalized_count, color=Symbol, fill=Symbol), stat="identity", position = position_dodge(), width=0.75) +
    xlab(NULL) +
    ylab("Normalized Counts") +
    scale_color_manual(values=colors) +
    scale_fill_manual(values=colors) +
    ggtitle("Human cell lines") +
    theme_pubr() +
    ylim(0,250)
  cancer_Hs <- cancer_Hs + theme(legend.position = "bottom")
  
  cancer_Mm <- ggplot(mouse_cell_lines)+
    geom_bar(aes(x=Sample, y=Normalized_count, color=mouse_gene, fill=mouse_gene), stat="identity", position = position_dodge(), width=0.75) +
    xlab(NULL) +
    ylab("Normalized Counts") +
    scale_color_manual(values=colors) +
    scale_fill_manual(values=colors) +
    ggtitle("Mouse cell lines") +
    theme_pubr() +
    ylim(0,250)
   cancer_Mm <- cancer_Mm + theme(legend.position = "bottom")
   
  cancer_plots <- ggarrange(cancer_Hs, cancer_Mm, ncol=2, align = "hv")
  return(cancer_plots)
}
```

```{r, fig.height=8,fig.width=8}
#pdf(file = "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_7/IGF2BPs_cell_lines.pdf",8,8)
making_plots_suppl()
#dev.off()
```



