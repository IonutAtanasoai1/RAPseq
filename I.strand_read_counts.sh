#!/bin/bash
#Made by Ionut Atanasoai
#Made on 210208
#UPPMAX commands
#SBATCH -A snic2020-15-293
#SBATCH -p core
#SBATCH -n 2



module load bioinfo-tools
module load BEDTools/2.29.2

NAMES=($(ls *9nt_window.common.bed | sed 's/\_/ /g' | awk '{print $1}'))

mv /proj/snic2020-16-228/P102/TRIMMED/UMIED/NO_rRNA_FASTQs_R1/SUBSETTED/ALIGNMENTS/*.bam ./
mv /proj/snic2020-16-228/P102/TRIMMED/UMIED/NO_rRNA_FASTQs_R1/SUBSETTED/ALIGNMENTS/*.bai ./


for i in $(seq 0 $(ls *9nt_window.common.bed | wc -l | awk '{print $1-1}'))

do

	cat $(echo ${NAMES[$i]})*9nt_window.common.bed > bed.txt
	bedtools bamtobed -split -i $(echo ${NAMES[$i]})*rep1*bam | grep '/1' | grep 'chr' | awk '{if($6=="-"){print $0}}' | bedtools sort > bam.minus.txt
	bedtools bamtobed -split -i $(echo ${NAMES[$i]})*rep1*bam | grep '/1' |	grep 'chr' | awk '{if($6=="+"){print $0}}' | bedtools sort > bam.plus.txt

	bedtools map -a bed.txt -b bam.plus.txt -o count -c 1 > A.txt
	bedtools map -a A.txt -b bam.minus.txt -o count -c 1 | awk '{if($6>=5 || $7>=5){print $1,$2,$3,$4,$5,$6+0.1,$7+0.1}}' | tr " " "\t" > B.txt
	awk '{if($6/$7>=10){print $1,$2,$3,$4,$5,"+",$6-0.1}}' B.txt | tr " " "\t" > B.plus.txt
	awk '{if($7/$6>=10){print $1,$2,$3,$4,$5,"-",$7-0.1}}' B.txt | tr " " "\t" > B.minus.txt
	cat B.plus.txt B.minus.txt | bedtools sort > B.rep1.stranded.counted.txt

	bedtools bamtobed -split -i $(echo ${NAMES[$i]})*rep2*bam | grep '/1' | grep 'chr' | bedtools sort > bam.2.txt
	bedtools map -a B.rep1.stranded.counted.txt -b bam.2.txt -s -o count -c 1 | awk '{if($8>=5){print $0}}' > $(echo ${NAMES[$i]})_stranded.counted.bed

	rm *txt


done

mkdir STRANDED_COUNTED
mv *stranded.counted.bed ./STRANDED_COUNTED/

mv *.bam /proj/snic2020-16-228/P102/TRIMMED/UMIED/NO_rRNA_FASTQs_R1/SUBSETTED/ALIGNMENTS/
mv *.bai /proj/snic2020-16-228/P102/TRIMMED/UMIED/NO_rRNA_FASTQs_R1/SUBSETTED/ALIGNMENTS/

