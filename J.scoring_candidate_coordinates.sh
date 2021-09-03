#!/bin/bash
#Made by Ionut Atanasoai
#Made on 210208
#UPPMAX commands
#SBATCH -A snic2020-15-293
#SBATCH -p core
#SBATCH -n 2



module load bioinfo-tools
module load BEDTools/2.29.2



mv /proj/snic2020-16-228/P102/TRIMMED/UMIED/NO_rRNA_FASTQs_R1/SUBSETTED/ALIGNMENTS/SCALED_BEDGRAPHs_CORRECT/MINUS_INPUT/PVALS_BIS/*pvalues*.bdg ./
mv /proj/snic2020-16-228/P102/TRIMMED/UMIED/NO_rRNA_FASTQs_R1/SUBSETTED/ALIGNMENTS/SCALED_BEDGRAPHs_CORRECT/*scaled.bedgraph ./
mv /proj/snic2020-16-228/P102/TRIMMED/UMIED/NO_rRNA_FASTQs_R1/SUBSETTED/ALIGNMENTS/SCALED_BEDGRAPHs_CORRECT/MINUS_INPUT/*bdg ./


tail -n+2 Halo*scaled.bedgraph | grep 'chr' | bedtools sort > Halo.bedgraph
tail -n+2 Input*scaled.bedgraph | grep 'chr' | bedtools sort > Input.bedgraph
tail -n+2 Halo*input.bdg | grep 'chr' | bedtools sort > Halo.sub.bedgraph


NAMES=($(ls *stranded.counted.bed | sed 's/\_/ /g' | awk '{print $1}' | sort -u))
for i in $(seq 0 $(ls *stranded.counted.bed | sed 's/\_/ /g' | awk '{print $1}' | sort -u | wc -l | awk '{print $1-1}'))


do
	
	
	
	tail -n+2 $(echo ${NAMES[$i]})*rep1*scaled.bedgraph | grep 'chr' | bedtools sort > rep1.bedgraph.txt
	tail -n+2 $(echo ${NAMES[$i]})*rep2*scaled.bedgraph | grep 'chr' | bedtools sort > rep2.bedgraph.txt
	tail -n+2 $(echo ${NAMES[$i]})*rep1*input.bdg | grep 'chr' | bedtools sort > rep1.subed.txt
    tail -n+2 $(echo ${NAMES[$i]})*rep2*input.bdg | grep 'chr' | bedtools sort > rep2.subed.txt


	bedtools sort -i $(echo ${NAMES[$i]})*rep1*pvalues*.bdg | grep 'chr' > sorted.pvalues_rep1.txt
	bedtools sort -i $(echo ${NAMES[$i]})*rep2*pvalues*.bdg | grep 'chr' > sorted.pvalues_rep2.txt

	bedtools sort -i $(echo ${NAMES[$i]})*stranded.counted.bed | grep 'chr' > sorted.txt 
	
	bedtools map -a sorted.txt -b sorted.pvalues_rep1.txt -c 4 -o max  > a.txt
	bedtools map -a a.txt -b sorted.pvalues_rep2.txt -c 4 -o max > b.txt

	bedtools map -a b.txt -b rep1.bedgraph.txt -c 4 -o mean > c.txt
	bedtools map -a c.txt -b rep2.bedgraph.txt -c 4 -o mean > d.txt
	bedtools map -a d.txt -b Halo.bedgraph -c 4 -o mean > e.txt
	bedtools map -a e.txt -b Input.bedgraph -c 4 -o mean > f.txt
	bedtools map -a f.txt -b rep1.subed.txt -c 4 -o mean > g.txt
	bedtools map -a	g.txt -b rep2.subed.txt	-c 4 -o mean > h.txt
	bedtools map -a h.txt -b Halo.sub.bedgraph -c 4 -o mean > i.txt


	cat i.txt > $(echo ${NAMES[$i]}).scored.bed

	rm *txt

done

rm Halo.bedgraph
rm Input.bedgraph
rm Halo.sub.bedgraph

mkdir SCORED
mv *scored.bed ./SCORED




mv *pvalues*.bdg /proj/snic2020-16-228/P102/TRIMMED/UMIED/NO_rRNA_FASTQs_R1/SUBSETTED/ALIGNMENTS/SCALED_BEDGRAPHs_CORRECT/MINUS_INPUT/PVALS_BIS/
mv *scaled.bedgraph /proj/snic2020-16-228/P102/TRIMMED/UMIED/NO_rRNA_FASTQs_R1/SUBSETTED/ALIGNMENTS/SCALED_BEDGRAPHs_CORRECT/
mv *input.bdg /proj/snic2020-16-228/P102/TRIMMED/UMIED/NO_rRNA_FASTQs_R1/SUBSETTED/ALIGNMENTS/SCALED_BEDGRAPHs_CORRECT/MINUS_INPUT/
