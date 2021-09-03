---
title: "Figure 6"
output:
  html_document:
    df_print: paged
---

# Author: Ionut Atanasoai

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

Figure 6A
```{r, fig.width=4, fig.height=4}
  
RIC <- read.csv(file = "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_6/RBPomes/RIC_table.csv", stringsAsFactors = F, header = T, sep=";")
RIC <- RIC[RIC$Hs_HEK293.RIC == "TRUE" | RIC$Hs_HuH7.RIC == "TRUE" | RIC$Hs_HeLa.RIC == "TRUE" | RIC$Hs_HeLa.RBDmap == "TRUE" | RIC$Hs_JURKAT.2018.RIC == "TRUE" | RIC$Hs_JURKAT.2018.eRIC == "TRUE" | RIC$Hs_Cyto.eRIC == "TRUE" | RIC$Hs_Nuc.eRIC == "TRUE",]
RIC <- unique(RIC$UNIQUE)

XRNAX <- read.csv(file = "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_6/RBPomes/XRNAX_table.csv", stringsAsFactors = F, header = T, sep=";")
XRNAX <- unique(XRNAX$Gene.name)

OOPS <- read.csv(file = "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_6/RBPomes/OOPS.csv", stringsAsFactors = F, header = T, sep=";")
OOPS <- unique(AnnotationDbi::select(org.Hs.eg.db, keys=OOPS$Uniprot_ID, columns = c('UNIPROT', 'SYMBOL'), keytype="UNIPROT")[,2])

PTex <- read.csv(file = "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_6/RBPomes/PTex.csv", stringsAsFactors = F, header = T, sep=";")
PTex <- unique(PTex$Gene.name)

AB <- list(PTex, XRNAX, RIC, OOPS)
names(AB) <- c("PTex", "XRNAX", "RIC", "OOPS")
v <- venn(AB)

#pdf(file = "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_6/venn_RPMomes.pdf",3,3)
plot(v, fills=c(alpha(acqua_greens[c(1,4,7,10)],0.5)), quantities=TRUE, edges=F)
#dev.off()


A <- intersect(RIC,XRNAX)
AA <- intersect(A,PTex)
AAA <- intersect(AA,OOPS)






RBDs <- read.table(file = "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_6/RBPomes/RBDs_in_RBP_Census.txt", stringsAsFactors = F, header = T, sep="\t")
RBDs <- paste(RBDs$Pfam.RNA.binding.domains, collapse = "|")
RBP_Census <- read.csv(file = "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_6/RBPomes/Canonical_RBPs_Gerstberger_et_al.csv", stringsAsFactors = F, header = T, sep=";")
ZZ <- RBP_Census[grep(RBDs,RBP_Census$domains.count.),]
RBP_Census <- RBP_Census[grep("established",RBP_Census$supporting.evidence....pubmed.ID.),]
RBP_Census <- unique(RBP_Census$gene.name)
RBP_Census <- unique(c(RBP_Census,ZZ$gene.name))
RBPDB <- read.table(file = "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_6/RBPomes/RBPDB.utoronto.txt", stringsAsFactors = F, header = T, sep="\t")
RBPDB <- RBPDB$Annotation.ID
RBPDB <- unique(AnnotationDbi::select(org.Hs.eg.db, keys=RBPDB, columns = c('ENSEMBL','SYMBOL'), keytype="ENSEMBL")[,2])
canonical_RBPs <- unique(c(RBP_Census, RBPDB))
AAA <- setdiff(AAA,canonical_RBPs)





UNIPROT_IDs <- AnnotationDbi::select(org.Hs.eg.db, keys=AAA, columns = c('SYMBOL','UNIPROT','GENENAME'), keytype="SYMBOL")
colnames(UNIPROT_IDs) <- c("gene_name","uniprot","Description")
BBB <- data.frame(c("NONE","NONE"), c("NONE","NONE"), c("NONE","NONE"), c("NONE","NONE"), c("NONE","NONE"), c("NONE","NONE"), c("NONE","NONE"))
colnames(BBB) <- c("uniprot", "length", "hmm.acc",  "hmm.name", "start", "end", "type")
for (i in UNIPROT_IDs$uniprot){
  B <- uniprot2pfam(i)
  BBB <- rbind(BBB,B)
}
BBB <- BBB[BBB$uniprot!="NONE",]
BBB <- merge(BBB,UNIPROT_IDs,by="uniprot")
gene_name <- paste(unique(BBB[grep(RBDs,BBB$hmm.name),"gene_name"]),collapse = "|")
BBB <- BBB[grep(gene_name,BBB$gene_name, invert=T),]

## at this stage there are 88 proteins left after removing high confident RBPs with RBDs based on published annotations. 
## After manual curating of each of the 88 entries, by removing splicing factors, elongation factors etc, there are 

manual_curation <- c("TCERG1|AQR|PRPF40A|PSIP1|GTF2F1|ERH|TCOF1|PRPF4B|TMA16|SCAF11|EIF3|ERH|NVL|RNH1|RPRD2|TCEA1|ZNF598")
BBB <- BBB[grep(manual_curation,BBB$gene_name, invert=T),]
ncRBPs <- unique(BBB$gene_name)

A <- intersect(RIC,XRNAX)
AA <- intersect(A,PTex)
AAA <- intersect(AA,OOPS)


cRBPs <- setdiff(AAA,ncRBPs)
intersection <- c(length(cRBPs),length(ncRBPs))
names(intersection) <- c("cRBPs","ncRBPs")


barplot(intersection, las=1, ylim=c(0,550))
text(x=0.7,y=length(cRBPs)+25,label=length(cRBPs), cex=2)
text(x=1.9,y=length(ncRBPs)+25,label=length(ncRBPs), cex=2)

```





Supplementary Figure S6 A
```{r, fig.width=8, fig.height=4}


PATH <- "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/RAP_Annotated_PEAKS/P102/"
PEAKS <- list.files(path = PATH)
Names <- unlist(str_split(PEAKS,"\\.peak"))[grep("txt",unlist(str_split(PEAKS,"\\.peak")),invert = T)]
Tables <- list()

for (i in 1:length(Names)){
  Tables[[i]] <- as.data.frame(read.table(paste(PATH,PEAKS[i],sep = ""), stringsAsFactors = F, header = T))
}

names(Tables) <- Names
Tables <- Tables[-47]

AAA <- Tables[["YBX3"]][1:3,]

for (i in names(Tables)){
  AAA <- rbind(AAA,Tables[[i]])
}

Tables <- AAA[-c(1:3),]
rm(AAA)

Tables <- Tables[grep("RBFOX|HuR|PTBP|IGF2BP|HNRNP|YBX",Tables$RBP),]
Tables <- Tables[grep("HuRPTBP1|IGF2BP1R|IGF2BP2b|IGF2BP3R|IGF2BP3I",Tables$RBP,invert=T),]
Tables$RBP <- factor(Tables$RBP, levels = names(table(Tables[,c("RBP")]))[order(-t(table(Tables[,c("RBP")])))])



cRBPs <- as.vector(table(Tables$RBP))
cRBPs_BSs <- Tables$Mean_FCH

can_BS_per_Gene <- c()
for (i in unique(Tables$RBP)){
  can_BS_per_Gene <- c(can_BS_per_Gene,mean(as.vector(table(Tables[Tables$RBP == i,"gene_name"]))))
}
names(can_BS_per_Gene) <- unique(Tables$RBP)




NonCan_RBPs <- read.table("/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/RAP_Annotated_PEAKS/ncRBPs/ncRBPs_all_peaks.txt", header = T, stringsAsFactors = F, sep="\t")



NonCan_RBPs$RBP <- factor(NonCan_RBPs$RBP, levels = names(table(NonCan_RBPs[,c("RBP")]))[order(-t(table(NonCan_RBPs[,c("RBP")])))])

AA <- unique(NonCan_RBPs[,c("RBP","gene_name")])

BS_per_Gene <- c()
for (i in unique(NonCan_RBPs$RBP)){
  BS_per_Gene <- c(BS_per_Gene,mean(as.vector(table(NonCan_RBPs[NonCan_RBPs$RBP == i,"gene_name"]))))
}
names(BS_per_Gene) <- unique(NonCan_RBPs$RBP)



#pdf(file = "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_6/barplot_n_sites.pdf",8,4)
barplot(table(NonCan_RBPs[,"RBP"]), las=2, ylim=c(0,max(table(NonCan_RBPs$RBP))+500), ylab="Number of Binding Sites")
abline(h=1000)
#dev.off()
#pdf(file = "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_6/barplot_n_sites_bellow200.pdf",4,3)
barplot(table(NonCan_RBPs$RBP)[table(NonCan_RBPs$RBP)<200], las=2)
#dev.off()

```




Supplementary Fig S6 B
```{r, fig.width=5, fig.height=10}

a <- table(NonCan_RBPs$RBP)[table(NonCan_RBPs$RBP)>1000]
AA <- NonCan_RBPs[grep(paste(names(a),collapse="|"),NonCan_RBPs$RBP),]



BB <- AA[1:2,]
medians <- c()
set.seed(1)
for (i in names(a)){
  B <- AA[AA$RBP == i,]
  B <- B[sample(1:nrow(B),1000),]
  medians <- c(medians, median(B$Mean_FCH))
  BB <- rbind(BB,B)
}
names(medians) <- names(a)
medians <- medians[order(-medians)]

AA$RBP <- factor(AA$RBP, levels=names(medians))
BB$RBP <- factor(BB$RBP, levels=names(medians))





#pdf(file = "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_6/boxplots_FCs_ncRBPs_above_1k.pdf",5,10)
par(mfrow=c(2,1))
boxplot(data=AA, log2(Mean_FCH)~RBP, las=2, outline=F, xlab=NA, ylab="log2(FC)  All Peaks")
abline(h=c(2,3,4,5), lty=2, lwd=0.5)
boxplot(data=BB, log2(Mean_FCH)~RBP, las=2, outline=F, xlab=NA, ylab="log2(FC)  Subsampled")
abline(h=c(2,3,4,5), lty=2, lwd=0.5)
#dev.off()

```



Figure 6B,C,D
```{r, fig.width=12, fig.height=4.5}

a <- table(NonCan_RBPs$RBP)[table(NonCan_RBPs$RBP)>1000]
b <- table(NonCan_RBPs$RBP)[table(NonCan_RBPs$RBP)<=1000]
AB <- list(cRBPs,as.vector(a),as.vector(b))
names(AB) <- c("canonical_RBPs","above_1k", "bellow_1k")


ncRBPs_BSs_high <- NonCan_RBPs[grep(paste(names(a),collapse="|"),NonCan_RBPs$RBP),"Mean_FCH"]
ncRBPs_BSs_low <- NonCan_RBPs[grep(paste(names(b),collapse="|"),NonCan_RBPs$RBP),"Mean_FCH"]
CD <- list(log2(cRBPs_BSs), log2(ncRBPs_BSs_high), log2(ncRBPs_BSs_low))
names(CD) <- c("canonical_RBPs","above_1k", "bellow_1k")

set.seed(1)
cRBPs_BSs_subs <- sample(cRBPs_BSs,1810)
ncRBPs_BSs_high_subs <- sample(ncRBPs_BSs_high, 1810)
CD_sub <- list(log2(cRBPs_BSs_subs), log2(ncRBPs_BSs_high_subs), log2(ncRBPs_BSs_low))
names(CD_sub) <- c("canonical_RBPs","above_1k", "bellow_1k")



above_1k <- names(table(NonCan_RBPs[,c("RBP")]))[order(-t(table(NonCan_RBPs[,c("RBP")])))][1:10]

above <- as.vector(BS_per_Gene[grepl(paste(above_1k,collapse = "|"),names(BS_per_Gene))])
bellow <- as.vector(BS_per_Gene[grepl(paste(above_1k,collapse = "|"),names(BS_per_Gene)) == "FALSE"])
can_BS_per_Gene <- as.vector(can_BS_per_Gene)

BSperGENE <- list(can_BS_per_Gene, above, bellow)
names(BSperGENE) <- c("cRBPs","above_1k", "bellow_1k")





#pdf(file = "/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_6/boxplots_fig6B.pdf",9,4)
par(mfrow=c(1,4))
boxplot2(AB, las=2, ylab="Number of RAPseq Peaks", lty=1)
abline(h=1000)
boxplot2(BSperGENE, las=2, ylab="Binding Sites per Gene (mean)", lty=1)

boxplot2(CD, las=2, ylab="log2(RBP vs Halo) - All Peaks", outline=F, lty=1, ylim=c(0,8.25))
boxplot2(CD_sub, las=2, ylab="log2(RBP vs Halo) - Downsampled", outline=F, lty=1, ylim=c(0,8.25))
#dev.off()



```


```{r, fig.width=10, fig.height=7}

AAA <- NonCan_RBPs
AAA$Biotype <- AAA$gene_type
AAA$Biotype <- paste(AAA$Biotype,AAA$feature, sep = "ZZZZ")
AAA <- AAA[grep("protein_codingZZZZexon",AAA$Biotype,invert=T),]
AAA <- AAA[grep("TEC",AAA$Biotype,invert=T),]
AAA$Biotype <- gsub("protein_codingZZZZ","",AAA$Biotype)
AAA$Biotype <- gsub("ZZZZexon","",AAA$Biotype)
AAA$Biotype <- gsub("ZZZZintron","",AAA$Biotype)
AAA$Biotype <- factor(AAA$Biotype, levels=rev(unique(AAA$Biotype)[order(unique(AAA$Biotype))]))
AAA$RBP <- factor(AAA$RBP, levels = names(table(NonCan_RBPs[,c("RBP")]))[order(-t(table(NonCan_RBPs[,c("RBP")])))])

aaa <- ggplot(data=AAA) + 
  geom_bar(aes(y=RBP, fill=Biotype), position="fill") +
  coord_flip() +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90),
    text = element_text(color="black")
    ) +
  xlab("%") +
  scale_fill_manual(values = c(pinks[c(1,6)], oranges[c(1,5)], yellows[c(2,6)], greys[5], acqua_blues[c(10,7,4,1)]))


BBB <- AAA
BBB$Biotype <- as.character(BBB$Biotype)
BBB <- BBB[grep("3UTR|5UTR|CDS|intron", BBB$Biotype, invert=T),]
BBB$Biotype <- factor(BBB$Biotype, levels=rev(unique(BBB$Biotype)[order(unique(BBB$Biotype))]))
BBB$RBP <- factor(BBB$RBP, levels = names(table(NonCan_RBPs[,c("RBP")]))[order(-t(table(NonCan_RBPs[,c("RBP")])))])

bbb <- ggplot(data=BBB) + 
  geom_bar(aes(y=RBP, fill=Biotype), position="fill") +
  coord_flip() +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90),
    text = element_text(color="black")
    ) +
  xlab("%") +
  scale_fill_manual(values = c(pinks[c(1,6)], oranges[c(1,5)], yellows[c(2,6)], greys[5]) )



ccc <- ggplot(data=BBB[BBB$Biotype != "lncRNA",]) + 
  geom_bar(aes(y=RBP, fill=Biotype)) +
  coord_flip() +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90),
    text = element_text(color="black")
    ) +
  scale_fill_manual(values = c(pinks[c(1,6)], oranges[c(1,5)], yellows[c(2,6)], greys[5]) )


#pdf(file="/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_6/RBPvsFeature_N_biotype_barplot.pdf", width=7, height=5)
ggarrange(aaa,bbb,ccc, nrow=3, common.legend = T, legend="right")
#dev.off()

```



```{r, fig.width=10, fig.height=8}


Inputs <- read.table(file="/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/INPUTs.TPMs.txt", stringsAsFactors = F, header=T)
Input_genes <- unlist(str_split(Inputs[Inputs$TPM_RAP >= 1,"gene_ID"],"\\."))[seq(1,length(unlist(str_split(Inputs[Inputs$TPM_RAP >= 1,"gene_ID"],"\\."))),2)]
ncRBP_bound_genes <- unique(unlist(str_split(NonCan_RBPs$gene_ID,"\\."))[seq(1,length(unlist(str_split(NonCan_RBPs$gene_ID,"\\."))),2)])
ALL_Genes <- unique(c(Input_genes,ncRBP_bound_genes))







Names <- unique(NonCan_RBPs$RBP)
Table <- data.frame(c(0,0), c(0,0), c(0,0), c(0,0), c(0,0))
npath <- c()
for (i in Names){
gene_IDs <- unlist(str_split(NonCan_RBPs[NonCan_RBPs$RBP == i,"gene_ID"],"\\."))[seq(1,length(unlist(str_split(NonCan_RBPs[NonCan_RBPs$RBP == i,"gene_ID"],"\\."))),2)]
BP <- enrichGO(gene = gene_IDs,
                universe = ALL_Genes,
                keyType       = "ENSEMBL",
                OrgDb         = org.Hs.eg.db,
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                minGSSize     = 10,
                maxGSSize     = 100,
                readable      = TRUE)

npath <- c(npath,nrow(as.data.frame(BP)))
}
counts <- data.frame(Names,npath)
counts <- counts[counts$npath>0,]
Names <- as.character(counts$Names)





Table <- data.frame(c("NONE","NONE"),c(0,0), c(0,0), c(0,0), c(0,0),c(0,0),c(0,0),c(0,0),c(0,0),c(0,0),c(0,0))
colnames(Table) <- c("ONTOLOGY","ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","geneID","Count","RBP")

for (i in Names){
gene_IDs <- unlist(str_split(NonCan_RBPs[NonCan_RBPs$RBP == i,"gene_ID"],"\\."))[seq(1,length(unlist(str_split(NonCan_RBPs[NonCan_RBPs$RBP == i,"gene_ID"],"\\."))),2)]
BP <- enrichGO(gene = gene_IDs,
                universe = ALL_Genes,
                keyType       = "ENSEMBL",
                OrgDb         = org.Hs.eg.db,
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                minGSSize     = 10,
                maxGSSize     = 100,
                readable      = TRUE)
BP1 <- as.data.frame(BP)
BP1$RBP <- rep(i,nrow(BP1))
Table <- rbind(Table,BP1)
print(nrow(Table))
}
Table <- Table[Table$ONTOLOGY != "NONE",]
Table$Count <- as.numeric(Table$Count)
Table <- Table[Table$Count>=5,]
Table$p.adjust <- as.numeric(Table$p.adjust)
A <- as.numeric(unlist(str_split(Table$GeneRatio,"\\/"))[seq(1,length( unlist(str_split(Table$GeneRatio,"\\/"))),2)])
B <- as.numeric(unlist(str_split(Table$GeneRatio,"\\/"))[seq(2,length( unlist(str_split(Table$GeneRatio,"\\/"))),2)])
Table$Percent_genes <- (A/B) * 100
Table$Description <- as.character(Table$Description)


BP_Table <- Table[Table$ONTOLOGY == "BP",]
MF_Table <- Table[Table$ONTOLOGY == "MF",]
CC_Table <- Table[Table$ONTOLOGY == "CC",]


##### collapse GO:BP ########
AA <- data.frame(c(0,0), c(0,0), c(0,0), c(0,0), c(0,0), c(0,0), c(0,0), c(0,0))
colnames(AA) <- c("Parent_ID","ParentTerm","N_Childs","ONTOLOGY","RBP","FDR","Gene_Count","Gene_names")
for(i in unique(BP_Table$RBP)){
  
  BP <- BP_Table[BP_Table$RBP == i,]
  BP$P_ID <- paste(BP$ID,BP$RBP,sep = "_")
  simMatrix <- calculateSimMatrix(BP$ID,
                                 orgdb="org.Hs.eg.db",
                                 ont="BP",
                                 method="Rel")
  
  scores <- setNames(-log10(BP$p.adjust), BP$ID)
  
  reducedTerms_BP <- reduceSimMatrix(simMatrix,
                                 scores,
                                 threshold=0.8,
                                orgdb="org.Hs.eg.db")
  
  reducedTerms_BP$RBP <- rep(i, nrow(reducedTerms_BP))
  reducedTerms_BP$Ontology <- rep("BP",nrow(reducedTerms_BP))
  reducedTerms_BP$P_ID <- paste(reducedTerms_BP$go,reducedTerms_BP$RBP,sep = "_")
  reducedTerms_BP <- reducedTerms_BP[,c("P_ID","Ontology","parent","parentTerm")]
  BP <- merge(BP,reducedTerms_BP,by="P_ID")

    
BB <- data.frame(c(0,0), c(0,0), c(0,0), c(0,0), c(0,0), c(0,0), c(0,0), c(0,0))
colnames(BB) <- c("Parent_ID","ParentTerm","N_Childs","ONTOLOGY","RBP","FDR","Gene_Count","Gene_names")
  for (j in unique(BP$parent)){
    
    FDR <- min(BP[BP$parent == j,"p.adjust"])
    A <- paste(BP[BP$parent == j,"geneID"],collapse = "/")
    Gene_names <- paste((unique(unlist(str_split(A,"\\/")))),collapse = "/")
    Gene_Count <- length(unique(unlist(str_split(A,"\\/"))))
    RBP <- unique(BP$RBP)
    ONTOLOGY <- unique(BP$ONTOLOGY)
    ParentTerm <- unique(BP[BP$parent == j,"parentTerm"])
    Parent_ID <- j
    N_Childs <- nrow(BP[BP$parent == j,])
    
    A <- data.frame(Parent_ID, ParentTerm, N_Childs, ONTOLOGY, RBP, FDR, Gene_Count, Gene_names)
    BB <- rbind(BB,A)
    
  }
BB <- BB[BB$Parent_ID != "0",]
AA <- rbind(AA,BB)

}
AA_BP <- AA[AA$Parent_ID != "0",]


##### collapse GO:MF ########
AA <- data.frame(c(0,0), c(0,0), c(0,0), c(0,0), c(0,0), c(0,0), c(0,0), c(0,0))
colnames(AA) <- c("Parent_ID","ParentTerm","N_Childs","ONTOLOGY","RBP","FDR","Gene_Count","Gene_names")
names <- names(table(MF_Table$RBP))[table(MF_Table$RBP)>1]
one <- paste(names(table(MF_Table$RBP))[table(MF_Table$RBP)==1],collapse = "|")
ones <- MF_Table[grep(one,MF_Table$RBP),]
ones$N_Childs <- rep(1,nrow(ones))
ones <- ones[,c("ID","Description","N_Childs","ONTOLOGY","RBP","p.adjust","Count","geneID")]
colnames(ones) <- c("Parent_ID","ParentTerm","N_Childs","ONTOLOGY","RBP","FDR","Gene_Count","Gene_names")
for(i in names){
  
  MF <- MF_Table[MF_Table$RBP == i,]
  MF$P_ID <- paste(MF$ID,MF$RBP,sep = "_")
  simMatrix <- calculateSimMatrix(MF$ID,
                                 orgdb="org.Hs.eg.db",
                                 ont="MF",
                                 method="Rel")
  
  scores <- setNames(-log10(MF$p.adjust), MF$ID)
  
  reducedTerms_MF <- reduceSimMatrix(simMatrix,
                                 scores,
                                 threshold=0.8,
                                 orgdb="org.Hs.eg.db")
  
  reducedTerms_MF$RBP <- rep(i, nrow(reducedTerms_MF))
  reducedTerms_MF$Ontology <- rep("MF",nrow(reducedTerms_MF))
  reducedTerms_MF$P_ID <- paste(reducedTerms_MF$go,reducedTerms_MF$RBP,sep = "_")
  reducedTerms_MF <- reducedTerms_MF[,c("P_ID","Ontology","parent","parentTerm")]
  MF <- merge(MF,reducedTerms_MF,by="P_ID")

    
BB <- data.frame(c(0,0), c(0,0), c(0,0), c(0,0), c(0,0), c(0,0), c(0,0), c(0,0))
colnames(BB) <- c("Parent_ID","ParentTerm","N_Childs","ONTOLOGY","RBP","FDR","Gene_Count","Gene_names")
  for (j in unique(MF$parent)){
    
    FDR <- min(MF[MF$parent == j,"p.adjust"])
    A <- paste(MF[MF$parent == j,"geneID"],collapse = "/")
    Gene_names <- paste((unique(unlist(str_split(A,"\\/")))),collapse = "/")
    Gene_Count <- length(unique(unlist(str_split(A,"\\/"))))
    RBP <- unique(MF$RBP)
    ONTOLOGY <- unique(MF$ONTOLOGY)
    ParentTerm <- unique(MF[MF$parent == j,"parentTerm"])
    Parent_ID <- j
    N_Childs <- nrow(MF[MF$parent == j,])
    
    A <- data.frame(Parent_ID, ParentTerm, N_Childs, ONTOLOGY, RBP, FDR, Gene_Count, Gene_names)
    BB <- rbind(BB,A)
    
  }
BB <- BB[BB$Parent_ID != "0",]
AA <- rbind(AA,BB)

}
AA_MF <- AA[AA$Parent_ID != "0",]
AA_MF <- rbind(AA_MF,ones)

##### collapse GO:CC ########
AA <- data.frame(c(0,0), c(0,0), c(0,0), c(0,0), c(0,0), c(0,0), c(0,0), c(0,0))
colnames(AA) <- c("Parent_ID","ParentTerm","N_Childs","ONTOLOGY","RBP","FDR","Gene_Count","Gene_names")
names <- names(table(CC_Table$RBP))[table(CC_Table$RBP)>1]
one <- paste(names(table(CC_Table$RBP))[table(CC_Table$RBP)==1],collapse = "|")
ones <- CC_Table[grep(one,CC_Table$RBP),]
ones$N_Childs <- rep(1,nrow(ones))
ones <- ones[,c("ID","Description","N_Childs","ONTOLOGY","RBP","p.adjust","Count","geneID")]
colnames(ones) <- c("Parent_ID","ParentTerm","N_Childs","ONTOLOGY","RBP","FDR","Gene_Count","Gene_names")
for(i in names){
  
  CC <- CC_Table[CC_Table$RBP == i,]
  CC$P_ID <- paste(CC$ID,CC$RBP,sep = "_")
  simMatrix <- calculateSimMatrix(CC$ID,
                                 orgdb="org.Hs.eg.db",
                                 ont="CC",
                                 method="Rel")
  
  scores <- setNames(-log10(CC$p.adjust), CC$ID)
  
  reducedTerms_CC <- reduceSimMatrix(simMatrix,
                                 scores,
                                 threshold=0.8,
                                orgdb="org.Hs.eg.db")
  
  reducedTerms_CC$RBP <- rep(i, nrow(reducedTerms_CC))
  reducedTerms_CC$Ontology <- rep("CC",nrow(reducedTerms_CC))
  reducedTerms_CC$P_ID <- paste(reducedTerms_CC$go,reducedTerms_CC$RBP,sep = "_")
  reducedTerms_CC <- reducedTerms_CC[,c("P_ID","Ontology","parent","parentTerm")]
  CC <- merge(CC,reducedTerms_CC,by="P_ID")

    
BB <- data.frame(c(0,0), c(0,0), c(0,0), c(0,0), c(0,0), c(0,0), c(0,0), c(0,0))
colnames(BB) <- c("Parent_ID","ParentTerm","N_Childs","ONTOLOGY","RBP","FDR","Gene_Count","Gene_names")
  for (j in unique(CC$parent)){
    
    FDR <- min(CC[CC$parent == j,"p.adjust"])
    A <- paste(CC[CC$parent == j,"geneID"],collapse = "/")
    Gene_names <- paste((unique(unlist(str_split(A,"\\/")))),collapse = "/")
    Gene_Count <- length(unique(unlist(str_split(A,"\\/"))))
    RBP <- unique(CC$RBP)
    ONTOLOGY <- unique(CC$ONTOLOGY)
    ParentTerm <- unique(CC[CC$parent == j,"parentTerm"])
    Parent_ID <- j
    N_Childs <- nrow(CC[CC$parent == j,])
    
    A <- data.frame(Parent_ID, ParentTerm, N_Childs, ONTOLOGY, RBP, FDR, Gene_Count, Gene_names)
    BB <- rbind(BB,A)
    
  }
BB <- BB[BB$Parent_ID != "0",]
AA <- rbind(AA,BB)

}
AA_CC <- AA[AA$Parent_ID != "0",]
AA_CC <- rbind(AA_CC,ones)

ALL_collapsed_GOs <- rbind(AA_BP,AA_MF,AA_CC)
colnames(ALL_collapsed_GOs) <- c("Parent_ID","GO_Term","N_Childs","ONTOLOGY","RBP","FDR","Gene_Count","Gene_names")
######   remove redundant terms that are found in more then one ncRBP    ##################
ALL_collapsed_GOs <- ALL_collapsed_GOs[grep("0071013|0003730|0010494|0072562|0035459|0061980|0030864",ALL_collapsed_GOs$Parent_ID,invert = T),]

ALL_collapsed_GOs$GO_Term <- factor(ALL_collapsed_GOs$GO_Term, levels = unique(ALL_collapsed_GOs[rev(order(ALL_collapsed_GOs$RBP)),"GO_Term"]))




GO_features <- data.frame(c(0,0), c(0,0))
colnames(GO_features) <- c("GO_Term","feature")
for (i in unique(ALL_collapsed_GOs$RBP))
  {
  
  AA <- ALL_collapsed_GOs[ALL_collapsed_GOs$RBP == i,]
  BB <- data.frame(c(0,0), c(0,0), c(0,0), c(0,0))
  colnames(BB) <- c("RBP_GENE","GO_Term","RBP","feature")
  for (j in AA$GO_Term)
    {
    ZZ <- AA[AA$GO_Term == j,"Gene_names"]
    gene_name <- unlist(str_split(ZZ,"\\/"))
    GO_Term <- rep(j,length(ZZ))
    RBP <- rep(i,length(ZZ))
    RBP_GENE <- paste(RBP,gene_name,sep="_")
    ZZ <- data.frame(GO_Term,RBP,RBP_GENE)
    
    YY <- NonCan_RBPs[NonCan_RBPs$RBP == i,]
    YY$RBP_GENE <-paste(YY$RBP,YY$gene_name,sep="_")
    YY <- YY[,c("RBP_GENE","feature")]
    
    ZZ <- merge(ZZ,YY,by="RBP_GENE")
    BB <- rbind(BB,ZZ)
  }
  AA <- merge(AA,BB,by="GO_Term")
  AA <- AA[,c("GO_Term","feature")]
  GO_features <- rbind(GO_features,AA)
  
}
GO_features$feature <- as.character(GO_features$feature)
GO_features <- GO_features[GO_features$GO_Term != "0",]
GO_features <- GO_features[GO_features$feature != "exon",]
GO_features$GO_Term <- factor(GO_features$GO_Term, levels=levels(ALL_collapsed_GOs$GO_Term))
GO_features <- as.data.frame(table(GO_features))
#colnames(GO_features) <- c("GO_Term","variable","value")
#GO_features <- cast(GO_features)

AA <- ALL_collapsed_GOs[,c("GO_Term","RBP")]
GO_features <- merge(GO_features,AA,by="GO_Term")



#pdf(file="/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_6/GOvsRBP_dotplot.pdf", width=8.5, height=7.5)
aaa <- ggplot2::ggplot(data=ALL_collapsed_GOs) + 
  geom_point( aes(x=GO_Term, y=RBP, alpha=FDR, color=ONTOLOGY, size=Gene_Count) ) + 
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_line(color="grey75", size=0.1),
    axis.text.x = element_text(angle = 90)
    ) +
  scale_color_manual(values = c(acqua_greens[5],"grey25",yellows[3])) +
  coord_flip() +
  scale_alpha_continuous(range = c(1, 0.4)) 
#aaa
#dev.off()



#pdf(file="/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_6/GOvsFeature_barplot.pdf", width=7, height=7.5)
bbb <- ggplot() + geom_bar(data=GO_features, aes(y=Freq, x=GO_Term, fill=feature), stat="identity", position="fill") + 
  coord_flip() +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90)
    ) +
  scale_fill_manual(values = c( oranges[4], pinks[3], acqua_blues[7], greys[4] ))  
#bbb
#dev.off()

aaa
bbb



```



```{r}



Inputs <- read.table(file="/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/INPUTs.TPMs.txt", stringsAsFactors = F, header=T)
Input_genes <- unlist(str_split(Inputs[Inputs$TPM_RAP >= 1,"gene_ID"],"\\."))[seq(1,length(unlist(str_split(Inputs[Inputs$TPM_RAP >= 1,"gene_ID"],"\\."))),2)]
ncRBP_bound_genes <- unique(unlist(str_split(NonCan_RBPs[,"gene_ID"],"\\."))[seq(1,length(unlist(str_split(NonCan_RBPs$gene_ID,"\\."))),2)])
ALL_Genes <- unique(c(Input_genes,ncRBP_bound_genes))



gene_IDs <- unlist(str_split(NonCan_RBPs[NonCan_RBPs$RBP == "MAPRE1","gene_ID"],"\\."))[seq(1,length(unlist(str_split(NonCan_RBPs[NonCan_RBPs$RBP == "MAPRE1","gene_ID"],"\\."))),2)]

BP <- enrichGO(gene = gene_IDs,
                universe = ALL_Genes,
                keyType       = "ENSEMBL",
                OrgDb         = org.Hs.eg.db,
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                minGSSize     = 10,
                maxGSSize     = 100,
                readable      = TRUE)



AA <- as.data.frame(BP)[,"geneID"]
gene_name <- unique(unlist(str_split(AA,"\\/")))
AA <- as.data.frame(gene_name)
MAPRE1 <- NonCan_RBPs[NonCan_RBPs$RBP == "MAPRE1",]
MAPRE1 <- unique(MAPRE1[,c("gene_name","Gene_BS")])
AA <- merge(AA,MAPRE1,by="gene_name")
BB <- AA$Gene_BS
names(BB) <- AA$gene_name

BP@result <- BP@result[grep("cortical cytoskeleton",BP@result$Description),]

aa <- cnetplot(BP, node_label="all",showCategory = 5, foldChange = log2(BB)) +
  scale_color_gradient(low = acqua_blues[9], high=acqua_blues[2])

pdf(file="/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_6/MAPRE1_RNAnetwork.pdf", width=5.25, height=3.5)
aa
dev.off()


aa


```


```{r}


Inputs <- read.table(file="/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/INPUTs.TPMs.txt", stringsAsFactors = F, header=T)
Input_genes <- unlist(str_split(Inputs[Inputs$TPM_RAP >= 1,"gene_ID"],"\\."))[seq(1,length(unlist(str_split(Inputs[Inputs$TPM_RAP >= 1,"gene_ID"],"\\."))),2)]
ncRBP_bound_genes <- unique(unlist(str_split(NonCan_RBPs[,"gene_ID"],"\\."))[seq(1,length(unlist(str_split(NonCan_RBPs$gene_ID,"\\."))),2)])
ALL_Genes <- unique(c(Input_genes,ncRBP_bound_genes))



gene_IDs <- unlist(str_split(NonCan_RBPs[NonCan_RBPs$RBP == "STMN1","gene_ID"],"\\."))[seq(1,length(unlist(str_split(NonCan_RBPs[NonCan_RBPs$RBP == "STMN1","gene_ID"],"\\."))),2)]

BP <- enrichGO(gene = gene_IDs,
                universe = ALL_Genes,
                keyType       = "ENSEMBL",
                OrgDb         = org.Hs.eg.db,
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                minGSSize     = 10,
                maxGSSize     = 100,
                readable      = TRUE)



AA <- as.data.frame(BP)[,"geneID"]
gene_name <- unique(unlist(str_split(AA,"\\/")))
AA <- as.data.frame(gene_name)
STMN1 <- NonCan_RBPs[NonCan_RBPs$RBP == "STMN1",]
STMN1 <- unique(STMN1[,c("gene_name","Gene_BS")])
AA <- merge(AA,STMN1,by="gene_name")
BB <- AA$Gene_BS
names(BB) <- AA$gene_name

BP@result <- BP@result[grep("migration",BP@result$Description),]

bb <- cnetplot(BP, node_label="all",showCategory = 5, foldChange = log2(BB)) +
  scale_color_gradient(low = acqua_blues[9], high=acqua_blues[2])


#pdf(file="/Users/ionut.atanasoai/Documents/RAPseq_Manuscript/FIGURE_6/STMN1_RNAnetwork.pdf", width=5.25, height=3.5)
bb
#dev.off()

bb

```

