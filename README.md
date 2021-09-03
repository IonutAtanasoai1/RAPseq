# RNA Affinity Purification and Sequencing (RAP-seq)
RAP-seq is an in vitro binding assay between a recombinant RNA Binding Protein (RBP) of interest and a pool of cellular extracted native total RNA. The pool of total RNA, previously to being incubated with the RBP is fragment to a normal distribution described by a median of ~40 nucleotides. After the binding assay is performed, the bound RNA molecules are recovered and cloned into an Illumina compatibile sequencing RNA library. Illumina NGS is used to sequence the bound molecules and deconvolute their identity by alignment to the respective reference genome from which the RNA was first extracted. The following lines describe what each script used for processing the NGS data and deconvoluting the identity of the bound RNA molecules does.

## Read Processing, Alignment and Peak Calling
RAP-seq NGS data processing involves 6 major steps that takes as input FASTQ Files and provides as output a BED File representing genomic coordinates of candidate RBP binding sites and their relative signals and significance in each of the RBP replicate binding assay and the two controls, namely the Input and Tag only libraries. These candidate binding sites need to be further filtered for corrected significance and fold change thresholds.

### 1. Read Processing (Trim Adapters) -> A.trim_PE.sh
The first step regards removal of the NEXTFLEX v3 small RNA library preparation kit adapters using A.trim_PE.sh. This bash script makes use of cutadapt to trim the adapters, perform a base quality filter and output all reads that are 20 or more nucleotides long. The length filter accounts for the presence of the library prep kit inserted UMIs at the 5' and 3' end of the reads.

### 2. Read Processing (UMI extraction) -> B.extract_UMIs.sh
B.extract_UMIs.sh makes use of awk to remove the first 4 bases of both read 1 and read 2, paste them together and append the 8 nucleotide UMI to both read 1 and read 2 read names.

### 3. Read Processing (Removal of rRNA and tRNA mapping reads) -> C.hisat2_remov_rRNA_tRNA.sh and D.select_read2.sh
The C.hisat2_remov_rRNA_tRNA.sh bash script uses hisat2 to align read 1 against a collection of rRNA (UCSC table browser, genome version hg38) and tRNA sequences (GtRNAdb). The collection of sequences is available in the RAP-seq manuscript supplementary files. The script will output all reads that do not aligh to rRNAs and tRNAs.
After the read 1 FASTQ File has beed depleted of rRNA and tRNA reads, the read 2 FASTQ File is parsed to select the mates that correspond to the ones filtered in the read 1 File using unix awk and join commands highlithed in the D.select_read2.sh script.

### 4. Alignment (Reference Genome Alignment) -> E.hisat2_genome_alignment_uniq.mapped_rmdup.sh
The E.hisat2_genome_alignment_uniq.mapped_rmdup.sh script makes use of hisat2, samtools and umi_tools to align the rRNA and tRNA depleted mates to a reference genome, retrieve uniquely mapping reads and remove PCR duplicates based on the UMIs. The alignment considers default hisat2 parameters and outputs bam formatted files. Samtools is used to convert the bam to sam and the "NH:i:1" flag reported by hisat2 is used to retrieve uniquely mapping reads. umi_tools dedup is then used to remove PCR duplicates.

### 5. Peak Calling (Identification Of Candidate Genomic Coordinates) -> F.make_bedgraphs_and_find_candidate_coordinates.sh
F.make_bedgraphs_and_find_candidate_coordinates.sh first computes pileups, subtracts signals from Input and HaloTag control assays and then builds a bed formatted file of candidate regions that need to be further subjected to filtering procedures. The bash script makes use of bedtools (genomecov and intersect commands are used) and macs2 (bdgcmp and bdgpeakcall) softwares.
In order to use this script, one need to provide a 2 column file, (or even 1 column, but modify the awk line of code at the beggining of the script) containing the scales in the second column that are going to be used by bedtools genomecov to compute bedgraphs. In the RAP-seq manuscript the scales used are present in the supplementary files. If the user is analyzing one or generally few assays she/he can skip the first lines and provide directly to bedtools the scale (this number will be used by bedtools to multiply the pileups by) for normalizing the coverage signals. In RAP-seq all coverage profiles are normalized to the sequencing depth of the HaloTag control library. 
After having built the bedgraphs, macs2 bgdcmp will take each file and subtract the signal present in Halo. The resulting subtraction bedgraphs will again be input of macs2 bdgcmp to also subtract the remaining signal from the Input library. After both subtractions, macs2 bdgpeakcall will output a large set of candidate genomic coordinates with a minimum length requirement of 20 nucleotides. 
These candidate regions are kept for future filtering if they are present in both replicates and their summits reside at most 10 nucleotides from each other. 
Considering that the replicates are used to identify candidate regions that will be further filtered, the user should always have ideally 2 replicates per RBP, otherwise the thresholds for filtering will need to be optimised in a case by case basis. The output of this script will be a 7 column bed file with the following columns: chromosome, start, end, summit replicate 1, summit replicate 2, width, distance between summits.

### 6. Peak Calling (Compute p-values) -> G.get_pvals_from_bedgraphs.sh


