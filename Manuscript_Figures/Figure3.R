---
title: "Figure 3"
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

Read in the metafile containing normalized counts and sequences

```{r}


HURs <- read.table(file = "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_3/Metafile/HuRs.metafile.ANNOTATED.final.bed", stringsAsFactors = F, header = T, sep="\t")
HURs$feature <- gsub("exon","ncRNA",HURs$feature)
AA <- HURs[HURs$phylo100_30mean != ".",]
BB <- HURs[HURs$phylo100_30mean == ".",]
BB$phylo100_30mean <- NA
HURs <- rbind(AA,BB)


AA <- HURs[HURs$phylo100_30mean_NEG != ".",]
BB <- HURs[HURs$phylo100_30mean_NEG == ".",]
BB$phylo100_30mean_NEG <- NA
HURs <- rbind(AA,BB)

HURs$phylo100_30mean <- as.numeric(HURs$phylo100_30mean)
HURs$phylo100_30mean_NEG <- as.numeric(HURs$phylo100_30mean_NEG)


```



write file for DREME de novo motif discovery
```{r}

orths <- c("hs","mm","md","gg","xt")
SEQ_RAPs <- list()

for (i in 1:length(orths)){
  
  SEQ_RAPs[[i]] <- HURs[grep(orths[i],HURs$RBP),c("positive_fa","negative_fa")]
  SEQ_RAPs[[i]]$positive_fa_sub <- str_sub(SEQ_RAPs[[i]]$positive_fa,75,125)
  SEQ_RAPs[[i]]$negative_fa_sub <- str_sub(SEQ_RAPs[[i]]$negative_fa,75,125)
  SEQ_RAPs[[i]]$header <- rep(">",nrow(SEQ_RAPs[[i]]))
  SEQ_RAPs[[i]]$seq_nr <- 1:nrow(SEQ_RAPs[[i]])
  SEQ_RAPs[[i]] <- SEQ_RAPs[[i]][,c("header","seq_nr","positive_fa_sub","negative_fa_sub")]
 
}



names(SEQ_RAPs) <- orths

#################################   write tables for dreme  ##############################


#PATH <- "/Users/ionutatanasoai/Documents/RAPseq_Manuscript/FIGURE_3/FOR_DREME/"
#sapply(names(SEQ_RAPs), function (x) write.table(SEQ_RAPs[[x]], file = paste(paste(PATH,x,sep = ""),".fastas.txt",sep = ""), row.names = F, col.names = F, sep = "\t") )  




## Command linde arguments used for DREME:
## hs : dreme -p positive_set_of_sequences_written_above.fa -n negative_set_of_sequences_written_above.fa -rna -m 1 -g 1
## mm : dreme -p positive_set_of_sequences_written_above.fa -n negative_set_of_sequences_written_above.fa -rna -m 1 -g 1
## md : dreme -p positive_set_of_sequences_written_above.fa -n negative_set_of_sequences_written_above.fa -rna -m 1 -g 1
## gg : dreme -p positive_set_of_sequences_written_above.fa -n negative_set_of_sequences_written_above.fa -rna -m 1 -g 1 -k 3
## xt : dreme -p positive_set_of_sequences_written_above.fa -n negative_set_of_sequences_written_above.fa -rna -m 1 -g 1


```



Figure 3D
```{r}

### check UUU centrality in all orthologs

plot( density(unlist(str_locate_all(HURs[grep("hs",HURs$RBP),"positive_fa"],"TTT"))-100,bw=2) , lwd=5, col=acqua_greens[1], bty="n", main=NA, xlab=NA, ylim=c(0.004,0.01), las=1)
points( density(unlist(str_locate_all(HURs[grep("mm",HURs$RBP),"positive_fa"],"TTT"))-100,bw=2), type="l", lwd=5, col=acqua_greens[3] )
points( density(unlist(str_locate_all(HURs[grep("md",HURs$RBP),"positive_fa"],"TTT"))-100,bw=2), type="l", lwd=5, col=acqua_greens[5] )
points( density(unlist(str_locate_all(HURs[grep("gg",HURs$RBP),"positive_fa"],"TTT"))-100,bw=2), type="l", lwd=5, col=acqua_greens[7] )
points( density(unlist(str_locate_all(HURs[grep("xt",HURs$RBP),"positive_fa"],"TTT"))-100,bw=2), type="l", lwd=5, col=acqua_greens[9] )

points( density(unlist(str_locate_all(HURs[grep("hs",HURs$RBP),"negative_fa"],"TTT"))-100,bw=2), type="l", lwd=5, col=greys[1] )
points( density(unlist(str_locate_all(HURs[grep("mm",HURs$RBP),"negative_fa"],"TTT"))-100,bw=2), type="l", lwd=5, col=greys[2] )
points( density(unlist(str_locate_all(HURs[grep("md",HURs$RBP),"negative_fa"],"TTT"))-100,bw=2), type="l", lwd=5, col=greys[3] )
points( density(unlist(str_locate_all(HURs[grep("gg",HURs$RBP),"negative_fa"],"TTT"))-100,bw=2), type="l", lwd=5, col=greys[4] )
points( density(unlist(str_locate_all(HURs[grep("xt",HURs$RBP),"negative_fa"],"TTT"))-100,bw=2), type="l", lwd=5, col=greys[5] )


```

Figure 3E

```{r}


HURs$UUUs <- str_count(str_sub(HURs$positive_fa,50,150),"TTT")
HURs$UUUs[HURs$UUUs > 7] <- 7
RAP <- HURs
a <- ggplot2::ggplot(data=RAP, aes(x=UUUs, y=Ortholog_Average_FCI)) + 
  stat_smooth(method="loess", color=acqua_blues[3], fill=acqua_blues[9], span=0.8, level=0.99) + 
  theme_classic2(base_size = 5) + 
  ylab("Fold Change") + 
  xlab("No of Uracil triplets") +
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7)) +
  coord_cartesian(ylim=c(3.4,6.5))
b <- ggplot2::ggplot(data=RAP, aes(x=as.character(UUUs), y=Ortholog_Average_FCI)) +
  geom_boxplot(outlier.shape = NA, fill=NA, color=acqua_blues[1], size=1) +
  theme_classic2(base_size = 5) + 
  ylab("Fold Change") + 
  xlab("No of Uracil triplets") +
  coord_cartesian(ylim=c(0,summary(RAP[RAP$UUUs == 6,"Ortholog_Average_FCI"])[[5]]*2))
AB <- ggarrange(a,b,NULL,ncol=3)






A <- ggplot2::ggplot(data=RAP, aes(x=UUUs, y=hs_Mean_FCI)) + 
  stat_smooth(method="loess", color=acqua_greens[1], fill=acqua_greens[2], span=0.8, level=0.99) + 
  theme_classic2(base_size = 5) + 
  ylab("Fold Change") + 
  xlab("No of Uracil triplets") +
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7)) +
  coord_cartesian(ylim = c(3,7))
B <- ggplot2::ggplot(data=RAP, aes(x=UUUs, y=mm_Mean_FCI)) + 
  stat_smooth(method="loess", color=acqua_greens[3], fill=acqua_greens[4], span=0.8, level=0.99) + 
  theme_classic2(base_size = 5) + 
  ylab("Fold Change") + 
  xlab("No of Uracil triplets") +
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7)) +
  coord_cartesian(ylim = c(3,7))
C <- ggplot2::ggplot(data=RAP, aes(x=UUUs, y=md_Mean_FCI)) + 
  stat_smooth(method="loess", color=acqua_greens[5], fill=acqua_greens[6], span=0.8, level=0.99) + 
  theme_classic2(base_size = 5) + 
  ylab("Fold Change") + 
  xlab("No of Uracil triplets") +
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7)) +
  coord_cartesian(ylim = c(3,7))
D <- ggplot2::ggplot(data=RAP, aes(x=UUUs, y=gg_Mean_FCI)) + 
  stat_smooth(method="loess", color=acqua_greens[7], fill=acqua_greens[8], span=0.8, level=0.99) + 
  theme_classic2(base_size = 5) + 
  ylab("Fold Change") + 
  xlab("No of Uracil triplets") +
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7)) +
  coord_cartesian(ylim = c(3,7))
E <- ggplot2::ggplot(data=RAP, aes(x=UUUs, y=xt_Mean_FCI)) + 
  stat_smooth(method="loess", color=acqua_greens[9], fill=acqua_greens[10], span=0.8, level=0.99) + 
  theme_classic2(base_size = 5) + 
  ylab("Fold Change") + 
  xlab("No of Uracil triplets") +
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7)) +
  coord_cartesian(ylim = c(3,7))

ABCDE <- ggarrange(A,B,C,D,E, ncol=5)
A <- ggarrange(AB,ABCDE,nrow=2)

#pdf(file = "/Users/ionutatanasoai/Documents/RAPseq_Manuscript/FIGURE_3/HURs_all_orths_FC~UUUs.pdf",9,4)
A
#dev.off()





```

Figure 3F

```{r,fig.width=6,fig.height=7}

### boxplots for seeing if the orthologs have putatively different binding affinitites

orths <- c("hs","mm","md","gg","xt")
FCIs <- list()
for (i in 1:5){
  aa <- paste(orths[i],"_Mean_FCI",sep="")
  bb <- paste(orths[i],"_Mean_FCH",sep="")
  FCIs[[i]] <- log2(HURs[HURs[,bb]>=2,aa]+1)
  
}
names(FCIs) <- orths


#pdf(file = "/Users/ionutatanasoai/Documents/RAPseq_Manuscript/FIGURE_3/HURs_boxplots_FCs_changes_in_affinity.pdf",5,6)
par(bty="n")
boxplot2(FCIs,outline=F, las=2, lty=1, range=1,  col=NA, ylab="log2(FC+1)", ylim=c(0,5), boxwex=0.9)
#dev.off()




```

Figure 3G
```{r, fig.width=6, fig.height=5}


HURs$N_orthologs <- as.character(HURs$N_orthologs)
HURs$N_orthologs <- factor(HURs$N_orthologs, levels = c("1","2","3","4","5"))


AAA <- HURs[1:5,]
set.seed(1)
for (i in c("1","2","3","4","5")){
  AA <- HURs[HURs$N_orthologs == i,]
  
  #### if you want to remove upper outliers use the next commented lines
  #a <- summary(AA$Ortholog_Average_BS)
  #b <- a[5] + (1.5*(a[5] - a[2]))
  #AA <- AA[AA$Ortholog_Average_BS <= b,]
  
  AA <- AA[sample( 1:nrow(AA),min(table(HURs$N_orthologs)) ),]
  AAA <- rbind(AAA,AA)
}
AAA <- AAA[-c(1:5),]


coordinates <- function(x) {
  data.frame(y = median(x, na.rm=T),
             ymin = quantile(x, na.rm=T)[[2]], 
             ymax = quantile(x, na.rm=T)[[4]]) 
}


aa <- ggplot(data=AAA, aes(x=N_orthologs, y=log2(Ortholog_Average_BS))) +
  geom_violin(aes(fill=N_orthologs), bw=0.75) +
  scale_fill_manual(values=acqua_greens[c(9,7,5,3,1)]) +
  stat_summary(fun.data = coordinates, geom = "pointrange", size = 1.5, na.rm=T) +
  theme_pubr() 

aa

```



Figure 3H
```{r, fig.width=6, fig.height=4}


AA <- HURs[1:2,]
set.seed(1)
for (i in 1:5){
  BB <- HURs[HURs$N_orthologs == i,]
  BB <- BB[sample( rownames(BB),min(table(HURs$N_orthologs)) ),]
  AA <- rbind(AA,BB)
}
AA <- AA[-c(1:2),]


coordinates <- function(x) {
  data.frame(y = median(x, na.rm=T),
             ymin = quantile(x, na.rm=T)[[2]], 
             ymax = quantile(x, na.rm=T)[[4]]) 
}

bb <- ggplot(data = AA, aes(x=as.character(N_orthologs) ,y=phylo100_30mean, color=as.character(N_orthologs))) +
  stat_summary(fun.data = coordinates, geom = "pointrange", size = 1.5, na.rm=T, position = position_nudge(x=-0.1, y=0))  +
  theme_classic() +
  scale_color_manual(values = acqua_greens[c(9,7,5,3,1)]) +
  ylab("Mean phyloP 100 vertebrates scores") + 
  guides(color=FALSE)

bb + stat_summary(data=AA , aes(x=as.character(N_orthologs),y=phylo100_30mean_NEG), color=greys[2], fun.data = coordinates, geom = "pointrange", size = 1.5, na.rm=T, position = position_nudge(x=0.1, y=0)) 


```


Figure 3I
```{r, fig.width=4,fig.height=4}

N_orths <- AA$N_orthologs
UUUs <- str_count(str_sub(AA$positive_fa,85,115),"TTT")
UUUs[UUUs>5] <- 5
UUUs <- factor(as.character(UUUs),levels = c("5","4","3","2","1","0"))
aa <- table(data.frame(N_orths,UUUs))
aa <- aa/rowSums(aa)

barplot( t(aa), col=c(yellows[c(1,3,5,7,9)],greys[7]), space=0.025, las=1 )

```




Figure 3J

```{r,fig.width=8,fig.height=8}

### pie charts for mrNA features 
HURs_FISH <- read.table(file = "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_3/Metafile/HuRs_FISH.metafile.ANNOTATED.final.bed", stringsAsFactors = F, header = T, sep="\t")
HURs_FISH$feature <- gsub("exon","ncRNA",HURs_FISH$feature)
AA <- HURs
AA$feature <- factor(AA$feature, levels=c("3UTR","5UTR","CDS","ncRNA","intron"))
BB <- HURs_FISH
BB$feature <- factor(BB$feature, levels=c("3UTR","5UTR","CDS","ncRNA","intron"))


#pdf(file = "/Users/ionutatanasoai/Documents/RAPseq_Manuscript/FIGURE_3/HURs_pie_Charts_Features.pdf",5,5)
par(mfrow=c(2,2), mar=c(2,2,2,2))
pie(x=100)
pie1(table(AA[AA$hs_Mean_FCH >= 2,"feature"]), main="hs in HUMAN", col=c("white",acqua_greens[c(4,7,11,1)]), percentage = T)
pie1(table(BB[BB$dr_Mean_FCH >= 2,"feature"]), main="dr in FISH", col=c("white",acqua_greens[c(4,7,11,1)]), percentage = T)
pie1(table(BB[BB$hs_Mean_FCH >= 2,"feature"]), main="hs in FISH", col=c("white",acqua_greens[c(4,7,11,1)]), percentage = T)
#dev.off()

```


hs_transcriptionally_active_genome_size <- 1987059556
dr_transcriptionally_active_genome_size <- 850149907






Write fastas for DREME for fish substrate assays

```{r}


#orths <- c("dr","hs")
#SEQ_RAPs <- list()

#for (i in 1:length(orths)){
  
#  SEQ_RAPs[[i]] <- HURs_FISH[grep(orths[i],HURs_FISH$RBP),c("positive_fa","negative_fa")]
#  SEQ_RAPs[[i]]$positive_fa_sub <- str_sub(SEQ_RAPs[[i]]$positive_fa,75,125)
#  SEQ_RAPs[[i]]$negative_fa_sub <- str_sub(SEQ_RAPs[[i]]$negative_fa,75,125)
#  SEQ_RAPs[[i]]$header <- rep(">",nrow(SEQ_RAPs[[i]]))
#  SEQ_RAPs[[i]]$seq_nr <- 1:nrow(SEQ_RAPs[[i]])
#  SEQ_RAPs[[i]] <- SEQ_RAPs[[i]][,c("header","seq_nr","positive_fa_sub","negative_fa_sub")]
 
#}



#names(SEQ_RAPs) <- orths

#################################   write tables for dreme  ##############################


#PATH <- "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_3/FOR_DREME/FISH/"
#sapply(names(SEQ_RAPs), function (x) write.table(SEQ_RAPs[[x]], file = paste(paste(PATH,x,sep = ""),".fastas.txt",sep = ""), row.names = F, #col.names = F, sep = "\t") )  



## Command linde arguments used for DREME:
## hs : dreme -p positive_set_of_sequences_written_above.fa -rna -g 1 -k 3 -s 2 -m 2
## dr : dreme -p positive_set_of_sequences_written_above.fa -rna -g 1 -k 3 -s 2 -m 2

```


Supplementary Figure 3E

```{r}



plot( density(unlist(str_locate_all(HURs_FISH[grep("hs",HURs_FISH$RBP),"positive_fa"],"TTT"))-100,bw=3) , lwd=5, col=acqua_blues[2], bty="n", main=NA, xlab=NA, ylim=c(0.003,0.01), las=1)
points( density(unlist(str_locate_all(HURs_FISH[grep("dr",HURs_FISH$RBP),"positive_fa"],"TTT"))-100,bw=3), type="l", lwd=5, col=acqua_blues[6] )

points( density(unlist(str_locate_all(HURs_FISH[grep("hs",HURs_FISH$RBP),"negative_fa"],"TTT"))-100,bw=3), type="l", lwd=5, col=greys[3] )
points( density(unlist(str_locate_all(HURs_FISH[grep("dr",HURs_FISH$RBP),"negative_fa"],"TTT"))-100,bw=3), type="l", lwd=5, col=greys[5] )




```



Supplementary FigureF

```{r}



HURs_FISH$UUUs <- str_count(str_sub(HURs_FISH$positive_fa,50,150),"TTT")
HURs_FISH$UUUs[HURs_FISH$UUUs > 7] <- 7
RAP <- HURs_FISH


par(bty="n", mfrow=c(1,2))
boxplot(data=RAP[RAP$hs_Mean_FCH>=2,], hs_Mean_FCI~UUUs, outline=F, range=1, col=acqua_blues[2], ylab="hs Fold Change", boxwex=0.9, lty=1)
boxplot(data=RAP[RAP$dr_Mean_FCH>=2,], dr_Mean_FCI~UUUs, outline=F, range=1, col=acqua_blues[6], ylab="dr Fold Change", boxwex=0.9, lty=1)


```





Supplementary Figure 3G

```{r,fig.height=7,fig.width=3}



orths <- c("hs","dr")
FCIs <- list()
for (i in 1:2){
  aa <- paste(orths[i],"_Mean_FCI",sep="")
  bb <- paste(orths[i],"_Mean_FCH",sep="")
  FCIs[[i]] <- log2(HURs_FISH[HURs_FISH[,bb]>=2,aa]+1)
  
}
names(FCIs) <- orths


#pdf(file = "/Users/ionutatanasoai/Documents/RAPseq_Manuscript/FIGURE_3/HURs_FISH_boxplots_FCs_changes_in_affinity.pdf",3,7.5)
par(bty="n")
boxplot2(FCIs,outline=F, las=2, lty=1, range=1,  col=NA, ylab="log2(FC+1)", ylim=c(0,7), boxwex=0.9)
#dev.off()


```




FISH GOs

```{r}

a <- unique(HURs_FISH[HURs_FISH$dr_Mean_FCH >= 2 ,"gene_name"])
entrez_IDs <- na.omit(as.data.frame(unlist(mapIds(org.Dr.eg.db, a, 'ENTREZID', 'SYMBOL')))[,1])


ALL_GOs_dr <- enrichGO(gene = entrez_IDs,
                keyType       = "ENTREZID",
                OrgDb         = org.Dr.eg.db,
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 1,
                minGSSize     = 10,
                readable      = TRUE)
ALL_GOs_dr <- as.data.frame(ALL_GOs_dr)


```

HUMAN GOs

```{r}


a <- unique(HURs[HURs$hs_Mean_FCH >= 2 ,"gene_name"])
entrez_IDs <- na.omit(as.data.frame(unlist(mapIds(org.Hs.eg.db, a, 'ENTREZID', 'SYMBOL')))[,1])


ALL_GOs_hs <- enrichGO(gene = entrez_IDs,
                keyType       = "ENTREZID",
                OrgDb         = org.Hs.eg.db,
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 1,
                minGSSize     = 10,
                readable      = TRUE)
ALL_GOs_hs <- as.data.frame(ALL_GOs_hs)


```


Supplementary Figure 3H upper panel
```{r}


ABC <- list(unique(ALL_GOs_hs$Description), unique(ALL_GOs_dr$Description))
names(ABC) <- c("Human", "Zebrafish")


v <- euler(ABC, shape="ellipse")
plot(v, fills=c(alpha(acqua_blues[2],0.7), alpha(acqua_blues[8],0.7), yellows[2]), quantities=TRUE, edges=F, main="GO child terms overlaps")


```

Considering only 3'UTR and intronic binding sites, what is their relative rations?
Supplementary Figure 3H lower panel
```{r}

common <- paste(intersect(unique(ALL_GOs_hs$ID), unique(ALL_GOs_dr$ID)),collapse = "|")

gene_name <- unique(unlist(str_split(paste(ALL_GOs_hs[grep(common,ALL_GOs_hs$ID),"geneID"],collapse = "/"),"\\/")))
a <- as.data.frame(gene_name)
AA <- HURs[HURs$hs_Mean_FCH >= 2,]
a <- merge(AA[grep("intron|3UTR",AA$feature), c("gene_name","feature")],a,by="gene_name")
hs_common_features <- table(a$feature)

gene_name <- unique(unlist(str_split(paste(ALL_GOs_dr[grep(common,ALL_GOs_dr$ID),"geneID"],collapse = "/"),"\\/")))
a <- as.data.frame(gene_name)
AA <- HURs_FISH[HURs_FISH$dr_Mean_FCH >= 2,]
a <- merge(AA[grep("intron|3UTR",AA$feature),c("gene_name","feature")],a,by="gene_name")
dr_common_features <- table(a$feature)



par(mfrow=c(1,2))
pie1(hs_common_features, col = c("white",acqua_blues[2]), main="human")
pie1(dr_common_features, col = c("white",acqua_blues[8]), main="zebrafish")


```

Supplementary Figure 3I

```{r, fig.height=5, fig.width=5}

common <- paste(intersect(unique(ALL_GOs_hs$ID), unique(ALL_GOs_dr$ID)),collapse = "|")





Pathways <- ALL_GOs_dr[grep(common,ALL_GOs_dr$ID),c("Description","geneID")]
colnames(Pathways) <- c("Pathway","geneID")
TTT <- data.frame()
for (i in 1:nrow(Pathways)){
  
  path <- Pathways[i,1]
  genes <- str_split(Pathways[i,2],"\\/")
  path <- rep(path,length(genes))
  DT <- data.frame(path,genes)
  colnames(DT) <- c("Pathway","geneID")
  TTT <- rbind(TTT,DT)
  
}
colnames(TTT) <- c("Description","gene_name")
TTT <- unique(TTT)

AA <- HURs_FISH[HURs_FISH$dr_Mean_FCH >= 2,]
AA <- AA[grep("intron|3UTR", AA$feature),c("gene_name","feature")] 
dr <- merge(TTT,AA,by="gene_name")
dr <- dr[,c("Description","feature")]
fish_features <- table(dr)/rowSums(table(dr))
fish_features <- fish_features[order(-fish_features[,1]),]

Pathways <- ALL_GOs_hs[grep(common,ALL_GOs_hs$ID),c("Description","geneID")]
colnames(Pathways) <- c("Pathway","geneID")
TTT <- data.frame()
for (i in 1:nrow(Pathways)){
  
  path <- Pathways[i,1]
  genes <- str_split(Pathways[i,2],"\\/")
  path <- rep(path,length(genes))
  DT <- data.frame(path,genes)
  colnames(DT) <- c("Pathway","geneID")
  TTT <- rbind(TTT,DT)
  
}
colnames(TTT) <- c("Description","gene_name")
TTT <- unique(TTT)

AA <- HURs[HURs$hs_Mean_FCH >= 2,]
AA <- AA[grep("intron|3UTR", AA$feature),c("gene_name","feature")] 
hs <- merge(TTT,AA,by="gene_name")
hs <- hs[,c("Description","feature")]
human_features <- table(hs)/rowSums(table(hs))
human_features <- human_features[rownames(fish_features),]

intron_ratios <- fish_features[,2]/human_features[,2]



par(bty="n")
plot(log2(intron_ratios[intron_ratios<=2^1.75]), fish_features[intron_ratios<=2^1.75, 2]*100, 
     ylim=c(0,100), 
     xlim=c(-4.5,4.5), 
     cex=3, 
     pch=16, 
     col=alpha(yellows[1],0.5), 
     las=1,
     xlab="Pathway Intron Ratios log2(Dr/Hs)",
     ylab="Pathway Intron % (Dr)")
points(log2(intron_ratios[intron_ratios>2^1.75]), fish_features[intron_ratios>2^1.75, 2]*100,  
       cex=3, 
       pch=16, 
       col=alpha(acqua_blues[5],0.7), 
       xaxt="n", 
       yaxt="n")

```

Supplementary Figure 3J

```{r, fig.height=5, fig.width=14}

par(mar=c(2,20,2,2), mfrow=c(1,2))
barplot(t(human_features[intron_ratios>2^1.75,])*100, horiz=T, las=1, col=c("white",acqua_blues[2]), space=0.05)
barplot(t(fish_features[intron_ratios>2^1.75,])*100, horiz=T, las=1, col=c("white",acqua_blues[8]), space=0.05)

```

```{r}
sessionInfo()
```
