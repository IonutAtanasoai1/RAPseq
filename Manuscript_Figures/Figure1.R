---
title: "Figure 1"
output:
  html_document:
    df_print: paged
---

# Authors: Ionut Atanasoai

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



```{r}

AA <- rnorm(n=5000, mean=110, sd=20)
BB <- rnorm(n=10000, mean=1870, sd=180)
CC <- rnorm(n=20000, mean=5070, sd=360)

ABC <- c(AA,BB,CC)
DD <- rnorm(n=20000, mean=70, sd=15)



par(mfrow=c(1,2))
plot(density(log2(ABC), bw=0.1), xlim=log2(c(20,10000)))
plot(density(log2(DD), bw=0.1), xlim=log2(c(20,10000)))




```



```{r, fig.width=5, fig.height=10}

qPCR <- read.table(file = "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_1/HUR_RAP_qPCR_results.txt", header = T, stringsAsFactors = F)
qPCR_means <- qPCR %>% group_by(Substrate, Target) %>% summarise(Mean_Fold_Change = mean(Fold_Change))
qPCR_means <- as.data.frame(qPCR_means)

aa <- ggplot2::ggplot() +
  geom_bar(data=qPCR_means, aes(x=Target, y=Mean_Fold_Change, fill=Target), stat="identity", color="black", alpha=0.9) +
  geom_point(data=qPCR, aes(x=Target, y=Fold_Change, fill=Target), size=3.5, pch=21, color="black") +
  facet_wrap(~Substrate) + 
  theme_pubr() + 
  theme(axis.text.x = element_text(angle=90)) +
  ylim(0,15) +
  scale_fill_manual(values=c(greys[5], acqua_greens[c(8,5,2)])) +
  geom_hline(yintercept = c(0,1,2,3,4,5,10,15), lty=2, size=0.1, color="black")

#pdf(file = "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_1/HUR_RAPqPCR_barplot.pdf",2.5,5)
aa
#dev.off()

```

