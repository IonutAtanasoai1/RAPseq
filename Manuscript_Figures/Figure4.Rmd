---
title: "Figure 4"
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


read in files
```{r}

HUR <- read.table("/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/RAP_Annotated_PEAKS/P102/hsHuRHUMAN.peaks.txt", stringsAsFactors = F, header = T)
PTBP1 <- read.table("/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/RAP_Annotated_PEAKS/P102/PTBP1.peaks.txt", stringsAsFactors = F, header = T)
HUR_PTBP1 <- read.table("/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/RAP_Annotated_PEAKS/P102/HuRPTBP1.peaks.txt", stringsAsFactors = F, header = T)

print(nrow(HUR))
print(nrow(PTBP1))
print(nrow(HUR_PTBP1))

```

Figure 4A
```{r,fig.width=4,fig.height=3}

HH <- HUR[,c(1,7,8,4:6)]
HH$Summit_start <- HH$Summit_start - 50
HH$Summit_end <- HH$Summit_end + 50
PP <- PTBP1[,c(1,7,8,4:6)]
PP$Summit_start <- PP$Summit_start - 50
PP$Summit_end <- PP$Summit_end + 50

hur <- makeGRangesFromDataFrame(HH)
ptbp1 <- makeGRangesFromDataFrame(PP)



HuR <- unique(HUR[-as.data.frame(findOverlaps(hur,ptbp1,type = "any"))[,1],"Peak_ID"])
common <- as.character(seq(1,nrow( as.data.frame(findOverlaps(hur,ptbp1,type = "any")) )))
HuR <- c(HuR,common)
Ptbp1 <- unique(PTBP1[-as.data.frame(findOverlaps(hur,ptbp1,type = "any"))[,2],"Peak_ID"])
Ptbp1 <- c(Ptbp1,common)

ABC <- list(HuR, Ptbp1)
names(ABC) <- c("HuR", "Ptbp1")
v <- euler(ABC, shape="circle")




h <- HUR[as.data.frame(findOverlaps(hur,ptbp1,type = "any"))[,1],"Summit_start"]
p <- PTBP1[as.data.frame(findOverlaps(hur,ptbp1,type = "any"))[,2],"Summit_start"]

#pdf(file = "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_4/venn_sites_HuRvsPTBP1.pdf",4,3)
plot(v, fills=c("white", "white", acqua_greens[5]), quantities=TRUE, edges=T, col=c(acqua_greens[c(3,8)]), lwd=4)
#dev.off()

#pdf(file = "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_4/competitive_cooperative_density_plot.pdf",4,3)
par(bty="n")
plot(density(log10(abs(p-h)), bw=0.15), xlab="log10(distance between binding sites)", main="bw = 0.15", ylim=c(0,1), las=1, lwd=4, xaxt="n")
axis(1,at=c(0,0.5,1,1.5,2))
abline(v=log10(30), lty=2)
text(x=1.7,y=0.95,label="30nt")
#dev.off()


```


FIGURE 4B
```{r,fig.width=6, fig.height=5}


HU <- HUR[as.data.frame(findOverlaps(hur,ptbp1,type = "any"))[,1],][,c("Rep1", "Rep2")]
PT <- PTBP1[as.data.frame(findOverlaps(hur,ptbp1,type = "any"))[,2],][,c("Rep1", "Rep2")]
HU_comp <- HU[log10(abs(p-h)+1)<1.492,]
PT_comp <- PT[log10(abs(p-h)+1)<1.492,]



HU_rep1 <- HU_comp$Rep1 
HU_rep2 <- HU_comp$Rep2 
PT_rep1 <- PT_comp$Rep1
PT_rep2 <- PT_comp$Rep2



HU_PT <- data.frame(HU_rep1, HU_rep2, PT_rep1, PT_rep2)


x <- HU_PT
group <- c(1,1,2,2)
y <- DGEList(counts=x, group = group)
design <- model.matrix(~group)
y <- estimateDisp(y,design)

fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit)


#fit <- glmFit(y, design)
#qlf <- glmLRT(fit)

X <-qlf$table
X$p.adjust <- -log10(p.adjust(X$PValue, method = "BH"))
X$PValue <- -log10(X$PValue)


HU_PT$HuR <- log2((HU_PT$HU_rep1 + HU_PT$HU_rep2)/2)
HU_PT$PTBP1 <- log2((HU_PT$PT_rep1 + HU_PT$PT_rep2)/2)
HU_PT$delta_BS <- abs(HU_PT$HuR - HU_PT$PTBP1)
HU_PT$FDR <- X$p.adjust
HU_PT$Significant <- HU_PT$FDR >= 1.30103
HU_PT$logFC <- X$logFC

HU_PT$HuR_won <- HU_PT$logFC <= -1 & HU_PT$FDR >= 1.30103
HU_PT$PTBP1_won <- HU_PT$logFC >= 1 & HU_PT$FDR >= 1.30103
HU_PT$DiffBind <- paste(HU_PT$HuR_won,HU_PT$PTBP1_won,sep = "_")
HU_PT$DiffBind <- gsub("FALSE_FALSE","N.S.",HU_PT$DiffBind)
HU_PT$DiffBind <- gsub("FALSE_TRUE","PTBP1",HU_PT$DiffBind)
HU_PT$DiffBind <- gsub("TRUE_FALSE","HUR",HU_PT$DiffBind)

aa <- ggplot(data=HU_PT) +
  geom_point(aes(x=PTBP1,y=HuR,color=DiffBind, alpha = abs(logFC), size=abs(logFC))) +
  scale_color_manual(values = c(acqua_greens[3],"grey75",acqua_greens[8])) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(size=2,fill=NA,color="black"),
        axis.ticks.length = unit(4,"mm"),
        axis.ticks = element_line(color="black",size=1),
        axis.text = element_text(color="black",size=20)) +
  xlim(2,15) +
  ylim(2,15) +
  xlab("PTBP1 Normalized Counts (log2)") +
  ylab("HUR Normalized Counts (log2)")

#pdf(file = "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_4/Diff_bind_competitive_sites.pdf",6,5)
aa + geom_text(aes(x=13,y=3), label = table(HU_PT$DiffBind)["PTBP1"], size=15, color=acqua_greens[8]) +
  geom_text(aes(x=4,y=14), label = table(HU_PT$DiffBind)["HUR"], size=15, color=acqua_greens[3])
#dev.off()

```


```{r, fig.width=4, fig.height=4}


# median sequenced fragment size HUR  =  50
# median sequenced fragment size HURPTBP1_coRAP  =  51
# median sequenced fragment size PTBP1  =  35


a <- PTBP1$end - PTBP1$start
A <- a - 35
b <- HUR$end - HUR$start
B <- b - 50
c <- HUR_PTBP1$end - HUR_PTBP1$start
C <- c - 51


abc <- list(a,b,c)
names(abc) <- c("PTBP1","HUR","HUR & PTBP1")
ABC <- list(A,B,C)
names(ABC) <- c("PTBP1","HUR","HUR & PTBP1")

boxplot(abc, outline=F, ylab="Peak width")
boxplot(ABC, outline=F, ylab="Peak width normalized by lib size")
abline(h=c(median(A),median(B),median(C)))


```

FIGURE 4E
```{r, fig.width=4, fig.height=5}

# median sequenced fragment size HUR  =  50
# median sequenced fragment size HURPTBP1  =  51
# median sequenced fragment size PTBP1  =  35

HH <- HUR[,c(1,7,8,4:6)]
HH$Summit_start <- HH$Summit_start - 50
HH$Summit_end <- HH$Summit_end + 50
PP <- PTBP1[,c(1,7,8,4:6)]
PP$Summit_start <- PP$Summit_start - 50
PP$Summit_end <- PP$Summit_end + 50
HP <- HUR_PTBP1[,c(1,7,8,4:6)]
HP$Summit_start <- HP$Summit_start - 50
HP$Summit_end <- HP$Summit_end + 50

hur <- makeGRangesFromDataFrame(HH)
ptbp1 <- makeGRangesFromDataFrame(PP)
hurptbp1 <- makeGRangesFromDataFrame(HP)


HuR <- unique(HUR[-as.data.frame(findOverlaps(hur,hurptbp1,type = "any"))[,1],"Peak_ID"])
common <- unique(HUR_PTBP1[as.data.frame(findOverlaps(hur,hurptbp1,type = "any"))[,2],"Peak_ID"])
HuR <- c(HuR,common)
coRAP1 <- unique(HUR_PTBP1[-as.data.frame(findOverlaps(hur,hurptbp1,type = "any"))[,2],"Peak_ID"])
coRAP1 <- c(coRAP1,common)
ABC <- list(HuR, coRAP1)
names(ABC) <- c("HuR", "coRAP")
v1 <- euler(ABC, shape="circle")

PT <- unique(PTBP1[-as.data.frame(findOverlaps(ptbp1,hurptbp1,type = "any"))[,1],"Peak_ID"])
common <- unique(HUR_PTBP1[as.data.frame(findOverlaps(ptbp1,hurptbp1,type = "any"))[,2],"Peak_ID"])
PT <- c(PT,common)
coRAP2 <- unique(HUR_PTBP1[-as.data.frame(findOverlaps(ptbp1,hurptbp1,type = "any"))[,2],"Peak_ID"])
coRAP2 <- c(coRAP2,common)
ABC <- list(PT, coRAP2)
names(ABC) <- c("PTBP1", "coRAP")
v2 <- euler(ABC, shape="circle")

#pdf(file = "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_4/venn_sites_HUR.vs.coRAP.pdf",3,4)
plot(v1, fills=c("white", "white", acqua_greens[3]), quantities=TRUE, edges=T, col=c(acqua_greens[c(3,3)]), lwd=4)
#dev.off()
#pdf(file = "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_4/venn_sites_PTBP1.vs.coRAP.pdf",3,4)
plot(v2, fills=c("white", "white", acqua_greens[8]), quantities=TRUE, edges=T, col=c(acqua_greens[c(8,8)]), lwd=4)
#dev.off()


CC <- HUR_PTBP1[as.data.frame(findOverlaps(hur,hurptbp1,type="any"))[,2],]
EE <- HUR_PTBP1[as.data.frame(findOverlaps(ptbp1,hurptbp1,type="any"))[,2],]

PT <- PTBP1[as.data.frame(findOverlaps(ptbp1,hurptbp1,type="any"))[,1],]
HU <- HUR[as.data.frame(findOverlaps(hur,hurptbp1,type="any"))[,1],]

a <- (PT$end - PT$start) - 35
b <- (EE$end-EE$start)-51
c <- (HU$end - HU$start) - 50
d <- (CC$end-CC$start)-51

ABC <- list(c,d,NULL,a,b)
names(ABC) <- c("HUR","coRAP",NA,"PTBP1","coRAP")
#pdf(file = "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_4/Peak_widths.pdf",4,5)
par(bty="n")
boxplot(ABC, outline=F, ylab="Peak width normalized by lib size", las=2, ylim=c(-50,200), col=c(acqua_greens[c(4,2,3,9,7)]), lwd=2)
text(x=1,y=median(c)+9, label=median(c), cex=1.5)
text(x=2,y=median(d)+9, label=median(d), cex=1.5)
text(x=4,y=median(a)+9, label=median(a), cex=1.5)
text(x=5,y=median(b)+9, label=median(b), cex=1.5)

pw <- wilcox.test(a,b, alternative = "less")
pwp <- pw$p.value
if ( pwp <= 0.001 ) {
to_add <- "***" 
} else if ( btp <= 0.01 ) {
to_add <- "**" 
} else if ( btp <= 0.05 ) {
to_add <- "*" 
} else {
to_add <- "n.s." 
}
text(1.5,200,labels = to_add, cex=2.5)

pw <- wilcox.test(c,d, alternative = "less")
pwp <- pw$p.value
if ( pwp <= 0.001 ) {
to_add <- "***" 
} else if ( btp <= 0.01 ) {
to_add <- "**" 
} else if ( btp <= 0.05 ) {
to_add <- "*" 
} else {
to_add <- "n.s." 
}
text(4.5,200,labels = to_add, cex=2.5)
#dev.off()

```
FIGURE 4F
```{r, fig.width=4,fig.height=3}


uniq1 <- unique(HUR_PTBP1[-as.data.frame(findOverlaps(hur,hurptbp1,type = "any"))[,2],])

HH <- uniq1[,c(1,7,8,4:6)]
HH$Summit_start <- HH$Summit_start - 50
HH$Summit_end <- HH$Summit_end + 50
uniq2 <- makeGRangesFromDataFrame(HH)


UNIQs <- unique(uniq1[-as.data.frame(findOverlaps(ptbp1,uniq2,type = "any"))[,2],])

UNIQs$positive_fa <- str_sub(UNIQs$positive_fa,43,157)


bases <- c("A","C","T","G")
kmers <- unite(as.data.frame(permutations(n=4,r=3,v=bases,repeats.allowed = T)),col = kmers,sep = "")$kmers
k3 <- paste("TTT",kmers,sep="")
k3 <- paste(k3,"TTT",sep="")
kmers <- unite(as.data.frame(permutations(n=4,r=2,v=bases,repeats.allowed = T)),col = kmers,sep = "")$kmers
k2 <- paste("TTT",kmers,sep="")
k2 <- paste(k2,"TTT",sep="")
k1 <- c("TTTTTTT","TTTATTT","TTTGTTT","TTTCTTT")
hur_ks <- c(k1,k2,k3)
#ptbp1_ks <- c("TTCTCT","TTTTCT","CTTTCT","TCTTCT","CTCTCT","TCCTCT","TGTTCT","CCTTCT","TGCTCT","CCCTCT","CGTTCT","CGCTCT")
ptbp1_ks <- c("TCTCT","TTTCT","TTTCT","CTTCT","TCTCT","CCTCT","GTTCT","CTTCT","GCTCT","CCTCT","GTTCT","GCTCT")


for (i in hur_ks){
  UNIQs$positive_fa <- gsub(i,"U",UNIQs$positive_fa)
}

for (i in ptbp1_ks){
  UNIQs$positive_fa <- gsub(i,"Y",UNIQs$positive_fa)
}

UNIQs$Ys <- str_count(UNIQs$positive_fa,"Y")
UNIQs$Us <- str_count(UNIQs$positive_fa,"U")
UNIQs$TTTs <- str_count(UNIQs$positive_fa,"TTT")

#table(UNIQs[,c("Ys","Us")])
#as.matrix(table(UNIQs[,c("Ys","Us")]))[-1,-1]

AAA <- UNIQs[UNIQs$Ys == 1 & UNIQs$Us == 1,]
AAA$Y_loc <- unlist(str_locate_all(AAA$positive_fa,"Y"))[seq(1,nrow(AAA)*2,2)]
AAA$U_loc <- unlist(str_locate_all(AAA$positive_fa,"U"))[seq(1,nrow(AAA)*2,2)]
AAA$distance <- abs(AAA$Y_loc - AAA$U_loc)


cc <- ggplot2::ggplot(data=AAA, aes(x=distance, y=log2(Mean_FCH))) + 
  stat_smooth(method="loess", color=acqua_greens[3], fill=acqua_greens[9], span=0.6, level=0.95, lwd=2) + 
  theme_classic2(base_size = 15) + 
  coord_cartesian(xlim=c(0,100), ylim=c(2.4,4.6)) +
  geom_vline(xintercept = c(27), color=acqua_greens[3], lty=2, lwd=1) +
  geom_text(aes(x=16,y=4.55), label = 27, size=7.5, col=acqua_greens[3]) +
  #geom_text(aes(x=37,y=4.55), label = 30, size=7.5, col=acqua_greens[3]) +
  xlab("distance between motifs (nt)") +
  ylab("log2(FC)")

#pdf(file = "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_4/FC~vdistance_between_motfs.pdf",4,3)
cc
#dev.off()

```

Supplementary FIGURE S4.E
```{r, fig.width=4,fig.height=3}


uniq1 <- unique(HUR_PTBP1[-as.data.frame(findOverlaps(hur,hurptbp1,type = "any"))[,2],])

HH <- uniq1[,c(1,7,8,4:6)]
HH$Summit_start <- HH$Summit_start - 50
HH$Summit_end <- HH$Summit_end + 50
uniq2 <- makeGRangesFromDataFrame(HH)


UNIQs <- unique(uniq1[-as.data.frame(findOverlaps(ptbp1,uniq2,type = "any"))[,2],])

UNIQs$negative_fa <- str_sub(UNIQs$negative_fa,43,157)


bases <- c("A","C","T","G")
kmers <- unite(as.data.frame(permutations(n=4,r=3,v=bases,repeats.allowed = T)),col = kmers,sep = "")$kmers
k3 <- paste("TTT",kmers,sep="")
k3 <- paste(k3,"TTT",sep="")
kmers <- unite(as.data.frame(permutations(n=4,r=2,v=bases,repeats.allowed = T)),col = kmers,sep = "")$kmers
k2 <- paste("TTT",kmers,sep="")
k2 <- paste(k2,"TTT",sep="")
k1 <- c("TTTTTTT","TTTATTT","TTTGTTT","TTTCTTT")
hur_ks <- c(k1,k2,k3)
#ptbp1_ks <- c("TTCTCT","TTTTCT","CTTTCT","TCTTCT","CTCTCT","TCCTCT","TGTTCT","CCTTCT","TGCTCT","CCCTCT","CGTTCT","CGCTCT")
ptbp1_ks <- c("TCTCT","TTTCT","TTTCT","CTTCT","TCTCT","CCTCT","GTTCT","CTTCT","GCTCT","CCTCT","GTTCT","GCTCT")


for (i in hur_ks){
  UNIQs$negative_fa <- gsub(i,"U",UNIQs$negative_fa)
}

for (i in ptbp1_ks){
  UNIQs$negative_fa <- gsub(i,"Y",UNIQs$negative_fa)
}

UNIQs$Ys <- str_count(UNIQs$negative_fa,"Y")
UNIQs$Us <- str_count(UNIQs$negative_fa,"U")
UNIQs$TTTs <- str_count(UNIQs$negative_fa,"TTT")

#table(UNIQs[,c("Ys","Us")])
#as.matrix(table(UNIQs[,c("Ys","Us")]))[-1,-1]

AAA <- UNIQs[UNIQs$Ys == 1 & UNIQs$Us == 1,]
AAA$Y_loc <- unlist(str_locate_all(AAA$negative_fa,"Y"))[seq(1,nrow(AAA)*2,2)]
AAA$U_loc <- unlist(str_locate_all(AAA$negative_fa,"U"))[seq(1,nrow(AAA)*2,2)]
AAA$distance <- abs(AAA$Y_loc - AAA$U_loc)


cc <- ggplot2::ggplot(data=AAA, aes(x=distance, y=log2(Mean_FCH))) + 
  stat_smooth(method="loess", color=greys[1], fill=greys[5], span=0.6, level=0.95, lwd=2) + 
  theme_classic2(base_size = 15) + 
  coord_cartesian(xlim=c(0,100), ylim=c(2.4,4.6)) +
  geom_vline(xintercept = c(27), color=greys[1], lty=2, lwd=1) +
  geom_text(aes(x=16,y=4.55), label = 27, size=7.5, col=greys[1]) +
  #geom_text(aes(x=37,y=4.55), label = 30, size=7.5, col=greys[3]) +
  xlab("distance between motifs (nt)") +
  ylab("log2(FC)")

#pdf(file = "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_4/FC~vdistance_between_motfs_CTRL_sites.pdf",4,3)
cc
#dev.off()

```
Supplementary FIGURE S4.F
```{r, fig.width=4,fig.height=3}

uniq1 <- unique(HUR[-as.data.frame(findOverlaps(hur,hurptbp1,type = "any"))[,1],])

HH <- uniq1[,c(1,7,8,4:6)]
HH$Summit_start <- HH$Summit_start - 50
HH$Summit_end <- HH$Summit_end + 50
uniq2 <- makeGRangesFromDataFrame(HH)


UNIQs <- unique(uniq1[-as.data.frame(findOverlaps(ptbp1,uniq2,type = "any"))[,2],])

UNIQs$positive_fa <- str_sub(UNIQs$positive_fa,43,157)


bases <- c("A","C","T","G")
kmers <- unite(as.data.frame(permutations(n=4,r=3,v=bases,repeats.allowed = T)),col = kmers,sep = "")$kmers
k3 <- paste("TTT",kmers,sep="")
k3 <- paste(k3,"TTT",sep="")
kmers <- unite(as.data.frame(permutations(n=4,r=2,v=bases,repeats.allowed = T)),col = kmers,sep = "")$kmers
k2 <- paste("TTT",kmers,sep="")
k2 <- paste(k2,"TTT",sep="")
k1 <- c("TTTTTTT","TTTATTT","TTTGTTT","TTTCTTT")
hur_ks <- c(k1,k2,k3)
#ptbp1_ks <- c("TTCTCT","TTTTCT","CTTTCT","TCTTCT","CTCTCT","TCCTCT","TGTTCT","CCTTCT","TGCTCT","CCCTCT","CGTTCT","CGCTCT")
ptbp1_ks <- c("TCTCT","TTTCT","TTTCT","CTTCT","TCTCT","CCTCT","GTTCT","CTTCT","GCTCT","CCTCT","GTTCT","GCTCT")


for (i in hur_ks){
  UNIQs$positive_fa <- gsub(i,"U",UNIQs$positive_fa)
}

for (i in ptbp1_ks){
  UNIQs$positive_fa <- gsub(i,"Y",UNIQs$positive_fa)
}

UNIQs$Ys <- str_count(UNIQs$positive_fa,"Y")
UNIQs$Us <- str_count(UNIQs$positive_fa,"U")
UNIQs$TTTs <- str_count(UNIQs$positive_fa,"TTT")

#table(UNIQs[,c("Ys","Us")])
#as.matrix(table(UNIQs[,c("Ys","Us")]))[-1,-1]

AAA <- UNIQs[UNIQs$Ys == 1 & UNIQs$Us == 1,]
AAA$Y_loc <- unlist(str_locate_all(AAA$positive_fa,"Y"))[seq(1,nrow(AAA)*2,2)]
AAA$U_loc <- unlist(str_locate_all(AAA$positive_fa,"U"))[seq(1,nrow(AAA)*2,2)]
AAA$distance <- abs(AAA$Y_loc - AAA$U_loc)


cc <- ggplot2::ggplot(data=AAA, aes(x=distance, y=log2(Mean_FCH))) + 
  stat_smooth(method="loess", color=acqua_blues[5], fill=greys[7], span=0.6, level=0.95, lwd=2) + 
  theme_classic2(base_size = 15) + 
  coord_cartesian(xlim=c(0,100), ylim=c(2.4,4.6)) +
  geom_vline(xintercept = c(27), color=acqua_blues[5], lty=2, lwd=1) +
  geom_text(aes(x=16,y=4.55), label = 27, size=7.5, col=acqua_blues[5]) +
  #geom_text(aes(x=37,y=4.55), label = 30, size=7.5, col=acqua_blues[3]) +
  xlab("distance between motifs (nt)") +
  ylab("log2(FC)")

#pdf(file = "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_4/FC~vdistance_between_motfs_HURonly_sites.pdf",4,3)
cc
#dev.off()

```

Supplementary FIGURE S4.G
```{r, fig.width=4,fig.height=3}

uniq1 <- unique(PTBP1[-as.data.frame(findOverlaps(ptbp1,hurptbp1,type = "any"))[,1],])

HH <- uniq1[,c(1,7,8,4:6)]
HH$Summit_start <- HH$Summit_start - 50
HH$Summit_end <- HH$Summit_end + 50
uniq2 <- makeGRangesFromDataFrame(HH)


UNIQs <- unique(uniq1[-as.data.frame(findOverlaps(hur,uniq2,type = "any"))[,2],])

UNIQs$positive_fa <- str_sub(UNIQs$positive_fa,43,157)


bases <- c("A","C","T","G")
kmers <- unite(as.data.frame(permutations(n=4,r=3,v=bases,repeats.allowed = T)),col = kmers,sep = "")$kmers
k3 <- paste("TTT",kmers,sep="")
k3 <- paste(k3,"TTT",sep="")
kmers <- unite(as.data.frame(permutations(n=4,r=2,v=bases,repeats.allowed = T)),col = kmers,sep = "")$kmers
k2 <- paste("TTT",kmers,sep="")
k2 <- paste(k2,"TTT",sep="")
k1 <- c("TTTTTTT","TTTATTT","TTTGTTT","TTTCTTT")
hur_ks <- c(k1,k2,k3)
#ptbp1_ks <- c("TTCTCT","TTTTCT","CTTTCT","TCTTCT","CTCTCT","TCCTCT","TGTTCT","CCTTCT","TGCTCT","CCCTCT","CGTTCT","CGCTCT")
ptbp1_ks <- c("TCTCT","TTTCT","TTTCT","CTTCT","TCTCT","CCTCT","GTTCT","CTTCT","GCTCT","CCTCT","GTTCT","GCTCT")


for (i in hur_ks){
  UNIQs$positive_fa <- gsub(i,"U",UNIQs$positive_fa)
}

for (i in ptbp1_ks){
  UNIQs$positive_fa <- gsub(i,"Y",UNIQs$positive_fa)
}

UNIQs$Ys <- str_count(UNIQs$positive_fa,"Y")
UNIQs$Us <- str_count(UNIQs$positive_fa,"U")
UNIQs$TTTs <- str_count(UNIQs$positive_fa,"TTT")

#table(UNIQs[,c("Ys","Us")])
#as.matrix(table(UNIQs[,c("Ys","Us")]))[-1,-1]

AAA <- UNIQs[UNIQs$Ys == 1 & UNIQs$Us == 1,]
AAA$Y_loc <- unlist(str_locate_all(AAA$positive_fa,"Y"))[seq(1,nrow(AAA)*2,2)]
AAA$U_loc <- unlist(str_locate_all(AAA$positive_fa,"U"))[seq(1,nrow(AAA)*2,2)]
AAA$distance <- abs(AAA$Y_loc - AAA$U_loc)


cc <- ggplot2::ggplot(data=AAA, aes(x=distance, y=log2(Mean_FCH))) + 
  stat_smooth(method="loess", color=acqua_blues[9], fill=greys[7], span=0.6, level=0.95, lwd=2) + 
  theme_classic2(base_size = 15) + 
  coord_cartesian(xlim=c(0,100), ylim=c(1,5)) +
  geom_vline(xintercept = c(27), color=acqua_blues[9], lty=2, lwd=1) +
  geom_text(aes(x=16,y=5), label = 27, size=7.5, col=acqua_blues[9]) +
  #geom_text(aes(x=37,y=4.55), label = 30, size=7.5, col=acqua_blues[3]) +
  xlab("distance between motifs (nt)") +
  ylab("log2(FC)")

#pdf(file = "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_4/FC~vdistance_between_motfs_PTBP1only_sites.pdf",4,3)
cc
#dev.off()


```


Supplementary FIGURE S4.A
```{r, fig.width=3, fig.height=3}


PT <- unique(PTBP1[,"gene_name"])
HU <- unique(HUR[,"gene_name"])
coRAP <- unique(HUR_PTBP1[,"gene_name"])
ABC <- list(HU,PT,coRAP)

names(ABC) <- c("HUR","PTBP1","coRAP")
v1 <- euler(ABC, shape="circle")

#pdf(file = "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_4/Gene_OVRLPS.pdf",3,3)
plot(v1, fills=c(acqua_blues[c(6,10)], acqua_blues[1]), quantities=TRUE, edges=T, col=c("white","white","white"), lwd=1, main="gene overlaps")
#dev.off()

PT <- unique(PTBP1[,"gene_name"])
HU <- unique(HUR[,"gene_name"])
coRAP <- unique(HUR_PTBP1[,"gene_name"])

coRAP_only <- setdiff(coRAP,HU)
coRAP_only <- setdiff(coRAP_only,PT)

```

Supplementary FIGURE S4.B left panel
```{r, fig.width=6, fig.height=4.5}

aa <- as.data.frame(intersect(HU, coRAP))
colnames(aa) <- "gene_name"

HU_C <- merge(HUR,aa,by="gene_name")
HU_C <- HU_C[,c("Rep1","Rep2","gene_name")]
BBB <- melt(HU_C) %>% group_by(gene_name,variable) %>% summarise(Gene_Counts = sum(value))
BBB <- as.data.frame(BBB)
HU_C <- reshape2::dcast(BBB, gene_name~variable)
colnames(HU_C) <- c("gene_name","HUR_1","HUR_2")

HUPT_C <- merge(HUR_PTBP1,aa,by="gene_name")
HUPT_C <- HUPT_C[,c("Rep1","Rep2","gene_name")]
BBB <- melt(HUPT_C) %>% group_by(gene_name,variable) %>% summarise(Gene_Counts = sum(value))
BBB <- as.data.frame(BBB)
HUPT_C <- reshape2::dcast(BBB, gene_name~variable)
colnames(HUPT_C) <- c("gene_name","coRAP_1","coRAP_2")

HUR_coRAP_commongenes <- merge(HU_C,HUPT_C,by="gene_name")



x <- HUR_coRAP_commongenes[,2:5]
group <- c(1,1,2,2)
y <- DGEList(counts=x, group = group)
design <- model.matrix(~group)
y <- estimateDisp(y,design)

fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit)

X <-qlf$table
X$p.adjust <- -log10(p.adjust(X$PValue, method = "BH"))
X$PValue <- -log10(X$PValue)

X$HUR_won <- X$logFC <= -1 & X$p.adjust >= 1.30103
X$coRAP_won <- X$logFC >= 1 & X$p.adjust >= 1.30103
X$DiffBind <- paste(X$HUR_won,X$coRAP_won,sep = "_")
X$DiffBind <- gsub("FALSE_FALSE","N.S.",X$DiffBind)
X$DiffBind <- gsub("FALSE_TRUE","coRAP",X$DiffBind)
X$DiffBind <- gsub("TRUE_FALSE","HUR",X$DiffBind)



aa <- ggplot(data=X) +
  geom_point(aes(x=logFC,y=p.adjust,color=DiffBind, alpha = abs(logFC), size=abs(logFC)), pch=16) +
  scale_color_manual(values = c(acqua_blues[1],acqua_blues[6],"grey75")) +
  theme_classic(base_size = 20) +
  ylab("FDR (-log10)") +
  xlab("FC (log2)") +
  ylim(0,25) +
  xlim(-6,6)
aaa <- aa + geom_text(aes(x=5,y=25), label = table(X$DiffBind)["coRAP"], size=10, color=acqua_blues[1]) +
  geom_text(aes(x=-5,y=25), label = table(X$DiffBind)["HUR"], size=10, color=acqua_blues[6])

#pdf(file = "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_4/Gene_diff_bind_HUR.vs.coRAP.pdf",6,4.5)
aaa
#dev.off()

HUR_DB <- X
coRAP_HU_diffb <- HUR_coRAP_commongenes[HUR_DB$DiffBind == "coRAP","gene_name"]

```

Supplementary FIGURE S4.B right panel
```{r, fig.width=6, fig.height=4.5}

aa <- as.data.frame(intersect(PT, coRAP))
colnames(aa) <- "gene_name"

HU_C <- merge(PTBP1,aa,by="gene_name")
HU_C <- HU_C[,c("Rep1","Rep2","gene_name")]
BBB <- melt(HU_C) %>% group_by(gene_name,variable) %>% summarise(Gene_Counts = sum(value))
BBB <- as.data.frame(BBB)
HU_C <- reshape2::dcast(BBB, gene_name~variable)
colnames(HU_C) <- c("gene_name","PTBP1_1","PTBP1_2")

HUPT_C <- merge(HUR_PTBP1,aa,by="gene_name")
HUPT_C <- HUPT_C[,c("Rep1","Rep2","gene_name")]
BBB <- melt(HUPT_C) %>% group_by(gene_name,variable) %>% summarise(Gene_Counts = sum(value))
BBB <- as.data.frame(BBB)
HUPT_C <- reshape2::dcast(BBB, gene_name~variable)
colnames(HUPT_C) <- c("gene_name","coRAP_1","coRAP_2")

HUR_coRAP_commongenes <- merge(HU_C,HUPT_C,by="gene_name")



x <- HUR_coRAP_commongenes[,2:5]
group <- c(1,1,2,2)
y <- DGEList(counts=x, group = group)
design <- model.matrix(~group)
y <- estimateDisp(y,design)

fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit)

X <-qlf$table
X$p.adjust <- -log10(p.adjust(X$PValue, method = "BH"))
X$PValue <- -log10(X$PValue)

X$HUR_won <- X$logFC <= -1 & X$p.adjust >= 1.30103
X$coRAP_won <- X$logFC >= 1 & X$p.adjust >= 1.30103
X$DiffBind <- paste(X$HUR_won,X$coRAP_won,sep = "_")
X$DiffBind <- gsub("FALSE_FALSE","N.S.",X$DiffBind)
X$DiffBind <- gsub("FALSE_TRUE","coRAP",X$DiffBind)
X$DiffBind <- gsub("TRUE_FALSE","PTBP1",X$DiffBind)


aa <- ggplot(data=X) +
  geom_point(aes(x=logFC,y=p.adjust,color=DiffBind, alpha = abs(logFC), size=abs(logFC)), pch=16) +
  scale_color_manual(values = c(acqua_blues[1],"grey75",acqua_blues[10])) +
  theme_classic(base_size = 20) +
  ylab("FDR (-log10)") +
  xlab("FC (log2)") +
  ylim(0,25) +
  xlim(-10,10)
bbb <- aa + geom_text(aes(x=7.5,y=25), label = table(X$DiffBind)["coRAP"], size=10, color=acqua_blues[1]) +
  geom_text(aes(x=-7.5,y=25), label = table(X$DiffBind)["PTBP1"], size=10, color=acqua_blues[10])

#pdf(file = "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_4/Gene_diff_bind_PT.vs.coRAP.pdf",6,4.5)
bbb
#dev.off()

PTBP1_DB <- X
coRAP_PT_diffb <- HUR_coRAP_commongenes[PTBP1_DB$DiffBind == "coRAP","gene_name"]

```
Supplementary FIGURE S4.C
```{r, fig.width=7.5, fig.height=3}


a <- unique(c(coRAP_only,coRAP_HU_diffb,coRAP_PT_diffb))
entrez_IDs <- na.omit(as.data.frame(unlist(mapIds(org.Hs.eg.db, a, 'ENTREZID', 'SYMBOL')))[,1])
all <- unique(c(HUR_PTBP1$gene_name, HUR$gene_name,PTBP1$gene_name))
entrez_IDs_all <- na.omit(as.data.frame(unlist(mapIds(org.Hs.eg.db, all, 'ENTREZID', 'SYMBOL')))[,1])

BPs <- enrichGO(
                gene = entrez_IDs,
                universe = entrez_IDs_all, 
                keyType       = "ENTREZID",
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                minGSSize     = 10,
                readable      = TRUE
                )
BP <- as.data.frame(BPs)
BP$Description <- factor(BP$Description, levels = BP$Description)

CCs <- enrichGO(
                gene = entrez_IDs,
                universe = entrez_IDs_all, 
                keyType       = "ENTREZID",
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                minGSSize     = 10,
                readable      = TRUE
                )
CC <- as.data.frame(CCs)
CC$Description <- factor(CC$Description, levels = CC$Description)


GOs <- rbind(BP,CC)
aa <- ggplot2::ggplot(data=GOs, aes(x=Description, y=1)) + 
  geom_point(aes(color=p.adjust), size=sqrt(GOs$Count), alpha=0.9) + 
  coord_flip() + 
  theme_classic(base_size = 17.5) + 
  scale_color_gradient(low=acqua_blues[3],high=acqua_blues[8]) +
  ylab(NULL) + 
  xlab("GO:BP           GO:CC") +
  theme(axis.text.x = element_blank(), 
        axis.line = element_blank(), 
        axis.ticks.x = element_blank(),
        )

#pdf(file = "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_4/GOs.for.coRAP.pdf",7.5,3)
aa
#dev.off()


#gos <- list()
#for (i in 1:nrow(GOs)){
  
#  gos[[i]] <- unlist(str_split(GOs[i,"geneID"],"\\/"))
  
#}
#names(gos) <- as.character(GOs$Description)
#upset(fromList(gos), nsets = nrow(GOs), keep.order = TRUE)

```

Supplementary FIGURE S4.D 
```{r, fig.height=5, fig.width=9}

gene_name <- unique(unlist(str_split(GOs$geneID,"\\/")))
genes <- as.data.frame(gene_name)

features <- HUR_PTBP1[,c("gene_name","feature","BS","Mean_FCI","Mean_FCH")]
features <- merge(genes,features,by="gene_name")

#pdf(file = "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_4/GOs.for.coRAP_features_n_FCs_for_features.pdf",9,5)
par(bty="n",mfrow=c(1,2))
pie1(table(features$feature), percentage=T, col=acqua_blues[c(1,4,7,11)])
boxplot2(data=features,log2(Mean_FCI)~feature, outline=F, col=acqua_blues[c(1,4,7,11)], las=2, range=1, boxwex=0.5, ylab="log2(FC)")
#dev.off()

```





```{r}
sessionInfo()
```

