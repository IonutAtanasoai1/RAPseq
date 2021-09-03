#!/bin/bash
#Made by Ionut Atanasoai
#Made on 20210218
#UPPMAX commands
#SBATCH -A snic2020-15-293
#SBATCH -p core
#SBATCH -n 2



module load bioinfo-tools
module load BEDTools/2.29.2


NAMES=($(ls *scored.bed | sed 's/\./ /g' | awk '{print $1}'))
for i in $(seq 0 $(ls *scored.bed | wc -l | awk '{print $1-1}'))


do

	grep -v 'chrM' $(echo ${NAMES[$i]})*scored.bed | awk '{if($2>96){print $0}}' > no_chrM.txt
	grep 'chrM' $(echo ${NAMES[$i]})*scored.bed | awk '{if($2>100 && $3<15460){print $0}}' > chrM.txt
	cat no_chrM.txt chrM.txt | bedtools sort > bed.txt

	awk '{print $1,$2-96,$3+96,$4,$5,$6}' bed.txt | tr " " "\t" > pos.txt
	awk '{print $1,$2+1000,$3+1000,$4,$5,$6}' pos.txt | tr " " "\t" > neg.txt
	bedtools getfasta -s -fi /proj/snic2020-16-228/P102/TRIMMED/UMIED/NO_rRNA_FASTQs_R1/SUBSETTED/GRCh38.primary_assembly.genome.fa -bed pos.txt | grep -v '>' > pos.seq
	bedtools getfasta -s -fi /proj/snic2020-16-228/P102/TRIMMED/UMIED/NO_rRNA_FASTQs_R1/SUBSETTED/GRCh38.primary_assembly.genome.fa	-bed neg.txt | grep -v '>' > neg.seq
	paste bed.txt pos.seq neg.seq > $(echo ${NAMES[$i]})_final.bed

	rm *.txt
	rm *.seq

done



mkdir FINAL
mv *final.bed FINAL/
