#!/bin/bash
#Made by Ionut Atanasoai
#Made on 210102
#UPPMAX commands
#SBATCH -A snic2020-15-293
#SBATCH -p core
#SBATCH -n 12


module load bioinfo-tools
module load HISAT2/2.2.1



NAMES=($(ls *.R1*fastq | sed 's/\.R1/\t/g' | awk '{print $1}' | sort -u))

Read1=($(ls *.R1*fastq))


mkdir NO_rRNA_FASTQs_R1

for i in $(seq 0 $(ls *.R1*fastq | wc -l | awk '{print $1-1}'))


do

        hisat2 -p 12 --un NO_rRNA_FASTQs_R1/$(echo ${NAMES[$i]}).R1.no_rRNA_tRNA.fastq -x rRNAs_tRNAs -U $(echo ${Read1[$i]}) -S $(echo ${NAMES[$i]}).sam
        rm *.sam

done
