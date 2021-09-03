#!/bin/bash
#Made by Ionut Atanasoai
#Made on 200513
#UPPMAX commands
#SBATCH -A snic2020-15-293 
#SBATCH -p core
#SBATCH -n 16




module load bioinfo-tools
module load cutadapt/2.3


NAMES=($(ls ./READ_1/*R1*.fastq.gz | sed 's/\.R1/\t/g' | awk '{print $1}' | sort -u | sed 's/\// /g' | awk '{print $3}'))


mkdir TRIMMED

for i in $(seq 0 $(ls ./READ_1/*R1*.fastq.gz | wc -l | awk '{print $1-1}'))
do
	cutadapt -j 16 -m 20 --trim-n --nextseq-trim=20 -a TGGAATTCTCGGGTGCCAAGG -A GATCGTCGGACTGTAGAACTCTGAAC -o ./TRIMMED/$(echo ${NAMES[$i]}).read1_trimmed.fastq -p ./TRIMMED/$(echo ${NAMES[$i]}).read2_trimmed.fastq ./READ_1/$(echo ${NAMES[$i]})*R1*fastq.gz ./READ_2/$(echo ${NAMES[$i]})*R2*fastq.gz
done

