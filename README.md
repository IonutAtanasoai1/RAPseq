# RAP-seq Read Processing, Alignment and Peak Calling and Threshold Filtering
RAP-seq NGS data processing involves 7 major steps that takes as input FASTQ Files and provides as output a BED File representing genomic coordinates of candidate RBP binding sites and their relative signals and significance in each of the RBP replicate binding assay and the two controls, namely the Input and Tag only libraries. These candidate binding sites need to be further filtered for corrected significance and fold change thresholds.

## Read Processing
