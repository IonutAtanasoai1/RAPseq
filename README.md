# RNA Affinity Purification and Sequencing (RAP-seq)
RAP-seq is an in vitro binding assay between a HaloTag fused recombinant RNA Binding Protein (RBP) of interest and a pool of cellular extracted native total RNA. The pool of total RNA, previously to being incubated with the RBP is fragment to a normal distribution described by a median of ~40 nucleotides. After the binding assay is performed, the bound RNA molecules are recovered and cloned into an Illumina compatibile sequencing RNA library. Illumina NGS is used to sequence the bound molecules and deconvolute their identity by alignment to the respective reference genome of the species from which the RNA was first extracted. The following lines describe what each script used for processing the NGS data and deconvoluting the identity of the bound RNA molecules does.

## Read Processing, Alignment and Peak Calling
RAP-seq NGS data processing involves 3 major steps (FASTQ File processing, Alignment to a reference genome and Peak calling) that takes as input FASTQ Files and provides as output a BED Files representing genomic coordinates of candidate RBP binding sites and their relative signals and significance in each of the RBP replicate binding assay and the two controls, namely the Input and HaloTag only libraries. These candidate binding sites need to be further annotated and filtered for multiple testing corrected significance and fold change thresholds.

## Read Processing
#### 1. Trim Adapters -> A.trim_PE.sh
The first step regards removal of the NEXTFLEX v3 small RNA library preparation kit adapters using A.trim_PE.sh. This bash script makes use of cutadapt to trim the adapters, perform a base quality filter and output all reads that are 20 or more nucleotides long. The length filter accounts for the presence of the library prep kit inserted UMIs at the 5' and 3' end of the RNA molecules.

#### 2. UMI Extraction -> B.extract_UMIs.sh
B.extract_UMIs.sh makes use of awk to remove the first 4 bases of both read 1 and read 2, paste them together and append the 8 nucleotide UMI to both read 1 and read 2 read names.

#### 3. Removal of rRNA and tRNA Mapping Reads -> C.hisat2_remov_rRNA_tRNA.sh and D.select_read2.sh
The C.hisat2_remov_rRNA_tRNA.sh bash script uses hisat2 to align read 1 against a collection of rRNA (UCSC table browser, genome version hg38) and tRNA sequences (GtRNAdb). The collection of sequences is available in the RAP-seq manuscript supplementary files. The script will output all reads that do not aligh to rRNAs and tRNAs.
After the read 1 FASTQ File has beed depleted of rRNA and tRNA reads, the read 2 FASTQ File is parsed to select the mates that correspond to the ones filtered in the read 1 File using unix awk and join commands highlithed in the D.select_read2.sh script.

## Alignment
#### 4. Reference Genome Alignment -> E.hisat2_genome_alignment_uniq.mapped_rmdup.sh
The E.hisat2_genome_alignment_uniq.mapped_rmdup.sh script makes use of hisat2, samtools and umi_tools to align the rRNA and tRNA depleted mates to a reference genome, retrieve uniquely mapping reads and remove PCR duplicates based on the UMIs. The alignment considers default hisat2 parameters and outputs bam formatted files. Samtools is used to convert the bam to sam and the "NH:i:1" flag reported by hisat2 is used to retrieve uniquely mapping reads. umi_tools dedup is then used to remove PCR duplicates.

## Peak Calling
#### 5. Identification Of Candidate Genomic Coordinates -> F.make_bedgraphs_and_find_candidate_coordinates.sh
F.make_bedgraphs_and_find_candidate_coordinates.sh first computes pileups, subtracts signals from Input and HaloTag control assays and then builds a bed formatted file of candidate regions that need to be further subjected to filtering procedures. The bash script makes use of bedtools (genomecov and intersect commands are used) and macs2 (bdgcmp and bdgpeakcall) softwares.
In order to use this script, one needs to provide a 2 column file, (or even 1 column, but modify the awk line of code at the beggining of the script) containing the scales in the second column that are going to be used by bedtools genomecov to compute bedgraphs. In the RAP-seq manuscript the scales used are present in the supplementary files. If the user is analyzing one or generally few assays she/he can skip the first lines and provide directly to bedtools the scale (this number will be used by bedtools to multiply the pileups by) for normalizing the coverage signals. In RAP-seq all coverage profiles are normalized to the sequencing depth of the HaloTag control library. 
After having built the bedgraphs, macs2 bgdcmp will take each file and subtract the signal present in Halo. The resulting subtraction bedgraphs will again be input of macs2 bdgcmp to also subtract the remaining signal from the Input library. After both subtractions, macs2 bdgpeakcall will output a large set of candidate genomic coordinates with a minimum length requirement of 20 nucleotides. 
These candidate regions are kept for future filtering if they are present in both replicates and their summits reside at most 10 nucleotides from each other. 
Considering that the replicates are used to identify candidate regions that will be further filtered, the user should always have ideally 2 replicates per RBP, otherwise the thresholds for filtering will need to be optimised in a case by case basis. The output of this script will be a 7 column bed file with the following columns: chromosome, start, end, summit replicate 1, summit replicate 2, width, distance between summits.

#### 6. Compute P-values -> G.get_pvals_from_bedgraphs.sh
G.get_pvals_from_bedgraphs.sh takes as input the pileups computed by bedtools genomecov in section 5 of this readme file and subtracts the Input signal from both RBP and HaloTag libraries and then used macs2 bdgcmp with the argument -m ppois to compute the significance of signal enrichments in the RBP libraries over the HaloTag control library.

#### 7. Assign Strand And Count Mates -> H.1.run_mean_summits.sh, H.2.mean_summits.R, I.strand_read_counts.sh
The first 2 scripts (H.1.run_mean_summits.sh and H.2.mean_summits.R, both written by Sofia Papavasileiou) are used to call and run R from bash to compute the average location of the summit based on the two replicate summits. It outputs a 5 column bed file with the 4th and the 5th column being the original start and end of the candidate coordinates and the 2nd and 3rd column the locations of the mean summit between the 2 replicates +/- 4 nucleotides. 
I.strand_read_counts.sh takes the bed file from above and its corresponding bam file and uses bedtools bamtobed to output one bed file for all read 1 (read 1 is the sense of the RNA molecule) mapping to the minus strand and one to the plus strand. bedtools map is then used to count the reads at each summit for both strands and if the reads on one strand are 10 times more than on the other then that strand is assigned to the respective summit and this needs to hold true for both replicates. The output is a bed file with the first 5 columns from the fist mean_summits.R and the 6th column representing the strand, the 7th the number of read 1 reads counted for replicate one and the 8th column those for replicate two. 

#### 8. Map Signals To Candidate Coordinates -> J.scoring_candidate_coordinates.sh
J.scoring_candidate_coordinates.sh uses bedtools map to score the summits with the scaled pileups from section 5 for both RBP replicates, HaloTag and Input controls and the respective replicate p-values computed in section 6. It adds 5 extra columns after the 8th column of the bed file output from section 7. 9th: Replicate 1 p-values (-log10 scale), 10th: Replicate 2 p-value (-log10 scale), 11th: Replicate 1 Normalized pileup signal, 12th: Replicate 2 Normalized pileup signal, 13th: HaloTag Normalized pileup signal and 14th: Input Normalized pileup signal.

#### 9. Get The Sequence Of Candidate Coordinates -> K.get_sequences.sh
To the 14 column bed file produced in section 8, two more columns will be added by K.get_sequences.sh representing the sequence underlying the candidate coordinate (column 15) and the sequence of control sites taken 1000 nucleotides 3' from the summit of the candidate region (column 16). Both sequence column represent 200 nucleotide bins centered around the peak summit (column 15) or 1000 nucleotides 3' of it (column 16). To retrieve the sequences bedtools getfasta was used and the strandness option -s was used to retrive strand specific sequences.

## Peak Annotation and Filtering
#### 10. Filter and Annotate Peaks -> L.Peak_annotation_n_filtering.R
In order to annotate and filter the final candidate coordinates from section 9, R is used to compute enrichments and correct for multiple testing the p-values computed by macs2 (section 6). The R package GenomicFeatures is used to parse a gencode GTF annotation file and annotate the filtered candidate coordinates with gene names, gene IDs and relative RNA features (CDS, ncRNA, introns, etc.). In order to ensure unique annotations, a feature priority is given in the following order: CDS > 3'UTR > 5'UTR > non coding RNA > intron. The L.Peak_annotation_n_filtering.R is an R script that contains all the code necessary to load all R packages used in the RAP-seq manuscript and to read in the candidate coordinates from section 9 and output a final filtered and annotated 32 column tabular file containing RBP binding sites and the final output of the RAP-seq pipeline. The tabular file has the following columns:
01. Chromosome 
02. Start 
03. End 
04. Peak ID 
05. Name of the RNA Binding Protein profiled 
06. Strand 
07. Summit start 
08. Summit end
09. Replicate 1 normalized pileup signal
10. Replicate 2 normalized pileup signal
11. Halo normalized pileup signal
12. Input normalized pileup signal
13. -log10pval_rep1 
14. -log10pval_rep2 
15. -log10FDR_rep1 BH multiple testing correction of the macs2 computed p values
16. -log10FDR_rep2 BH multiple testing correction of the macs2 computed p values
17. Fold Change of the normalized pileup signals over Halo for Replicate 1
18. Fold Change of the normalized pileup signals over Halo for Replicate 2
19. Fold Change of the normalized pileup signals over Input for Replicate 1
20. Fold Change of the normalized pileup signals over Input for Replicate 2
21. The average Fold Change over Halo 
22. The average Fold Change over Input
23. Binding Score = The average Fold Change considering both Halo and Input controls weighed by the log2 normailzed average pileup signals of both replicates
24. Gene ENSEMBL ID 
25. Gene Name
26. Gene Type 
27. RNA feature (intron, CDS, UTR, etc) 
28. Gene Binding Score = The sum of all Binding Scores of all Binding Sites mapping to one gene 
29. Local IDR computed using the replicate -log10pvalues 
30. Global IDR computed using the replicate -log10pvalues
31. The sequence surrounding the RBP binsing site (200 nucleotide bin centered around the peak summit) 
33. The sequence of adjacent control unbound sites (200 nucleotide bin taken 1000 nucleotides 3' of the binding site)
