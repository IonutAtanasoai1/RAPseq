# RAP-seq Read Processing, Alignment and Peak Calling
RAP-seq NGS data processing involves 7 major steps that takes as input FASTQ Files and provides as output a BED File representing genomic coordinates of candidate RBP binding sites and their relative signals and significance in each of the RBP replicate binding assay and the two controls, namely the Input and Tag only libraries. These candidate binding sites need to be further filtered for corrected significance and fold change thresholds.

## 1. Read Processing
The first step regards removal of the NEXTFLEX v3 small RNA library preparation kit adapters using A.trim_PE.sh. This bash script makes use of cutadapt to trim the adapters, perform a base quality filter and output all reads that are 20 or more nucleotides long. The length filter accounts for the presence of the library prep kit inserted UMIs at the 5' and 3' end of the reads.

## 2. UMI extraction 

