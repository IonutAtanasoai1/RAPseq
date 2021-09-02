# RAP-seq Read Processing, Alignment and Peak Calling
RAP-seq NGS data processing involves 6 major steps that takes as input FASTQ Files and provides as output a BED File representing genomic coordinates of candidate RBP binding sites and their relative signals and significance in each of the RBP replicate binding assay and the two controls, namely the Input and Tag only libraries. These candidate binding sites need to be further filtered for corrected significance and fold change thresholds.

## 1. Read Processing -> A.trim_PE.sh
The first step regards removal of the NEXTFLEX v3 small RNA library preparation kit adapters using A.trim_PE.sh. This bash script makes use of cutadapt to trim the adapters, perform a base quality filter and output all reads that are 20 or more nucleotides long. The length filter accounts for the presence of the library prep kit inserted UMIs at the 5' and 3' end of the reads.

## 2. UMI extraction -> B.extract_UMIs.sh
B.extract_UMIs.sh makes use of awk to remove the first 4 bases of both read 1 and read 2, paste them together and append the 8 nucleotide UMI to both read 1 and read 2 read names.

## 3. Removal of rRNA and tRNA mapping reads -> C.hisat2_remov_rRNA_tRNA.sh and D.select_read2.sh
The C.hisat2_remov_rRNA_tRNA.sh bash script uses hisat2 to align read 1 against a collection of rRNA (UCSC table browser, genome version hg38) and tRNA sequences (GtRNAdb). The collection of sequences is available in the RAP-seq manuscript supplementary files. The script will output all reads that do not aligh to rRNAs and tRNAs.

After the read 1 FASTQ File has beed depleted of rRNA and tRNA reads, the read 2 FASTQ File is parsed to select the mates that correspond to the ones filtered in the read 1 File using unix awk and join commands highlithed in the D.select_read2.sh script.

## 4. 

