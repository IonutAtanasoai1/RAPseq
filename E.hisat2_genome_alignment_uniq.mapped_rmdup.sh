#!/bin/bash
#Made by Ionut Atanasoai
#Made on 210106
#UPPMAX commands
#SBATCH -A snic2020-15-293
#SBATCH -p core
#SBATCH -n 12


module load bioinfo-tools
module load HISAT2/2.1.1
module load samtools/1.9
module load umi_tools/1.0.0

NAMES=($(ls *.R1.*fastq | sed 's/\./\t/g' | awk '{print $1}' | sort -u))

Read1=($(ls *.R1.*fastq))
Read2=($(ls *.R2.*fastq))


touch unique_mapping_reads.txt
touch unique_mapping_reads_without_PCRdup.txt


for i in $(seq 0 $(ls *.R1.*fastq | wc -l | awk '{print $1-1}'))


do



        hisat2 -p 12 --no-unal -x GRCh38.primary_assembly.genome -1 $(echo ${Read1[$i]}) -2 $(echo ${Read2[$i]}) -S sam
        samtools view -H sam > header
        samtools view -@ 12 sam | grep 'NH:i:1' > unique
        cat header unique > $(echo ${NAMES[$i]}).sam
        rm sam
        rm header
        rm unique

        samtools view -@ 12 -b $(echo ${NAMES[$i]}).sam -o $(echo ${NAMES[$i]}).bam
        samtools sort -@ 12 $(echo ${NAMES[$i]}).bam -o $(echo ${NAMES[$i]}).sorted.bam
        samtools index $(echo ${NAMES[$i]}).sorted.bam
        samtools view -@ 12 -c $(echo ${NAMES[$i]}).sorted.bam &>> unique_mapping_reads.txt
        rm $(echo ${NAMES[$i]}).sam
        rm $(echo ${NAMES[$i]}).bam

        umi_tools dedup -I $(echo ${NAMES[$i]}).sorted.bam --paired -S $(echo ${NAMES[$i]}).umied_dedup.bam
        samtools index $(echo ${NAMES[$i]}).umied_dedup.bam
        samtools view -@ 12 -c $(echo ${NAMES[$i]}).umied_dedup.bam &>> unique_mapping_reads_without_PCRdup.txt
        rm  $(echo ${NAMES[$i]}).sorted.bam
        rm  $(echo ${NAMES[$i]}).sorted*bai



done

mkdir ALIGNMENTS
mv *.umied_dedup.ba* ALIGNMENTS
mv unique_mapping_reads.txt ALIGNMENTS
mv unique_mapping_reads_without_PCRdup.txt ALIGNMENTS
