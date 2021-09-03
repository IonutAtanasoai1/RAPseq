#!/bin/bash
#Made by Ionut Atanasoai
#Made on 210112
#UPPMAX commands
#SBATCH -A snic2020-15-293
#SBATCH -p core
#SBATCH -n 2




module load bioinfo-tools
module load samtools/1.9
module load BEDTools/2.29.2


BAMs=($(ls *.bam))
rm scales.txt
touch scales.txt



cat /proj/snic2020-16-228/SCALES_p102.txt | awk '{print $2}' &>> scales.txt




RPMs=($(cat scales.txt))



for i in $(seq 0 $(ls *bam | wc -l | awk '{print $1-1}'))

do

        bedtools genomecov -bga -split -trackline -scale $(echo ${RPMs[$i]}) -ibam $(echo ${BAMs[$i]}) > $(echo ${BAMs[$i]}).scaled.bedgraph

done

mkdir SCALED_BEDGRAPHs_CORRECT
mv *scaled.bedgraph SCALED_BEDGRAPHs_CORRECT/






cd SCALED_BEDGRAPHs_CORRECT


module load bioinfo-tools
module load MACS/2.2.6
module load BEDTools/2.29.2



## Substracting Halo

NAMES=($(ls *.bedgraph | grep -v 'Halo' | sed 's/\./ /g' | awk '{print $1}'))
for i in $(seq 0 $(ls *.bedgraph | grep -v 'Halo' | wc -l | awk '{print $1-1}'))


do


	tail -n+2 $(echo ${NAMES[$i]})*bedgraph > aaa.bedgraph
	tail -n+2 Halo*bedgraph > bbb.bedgraph

	macs2 bdgcmp -t aaa.bedgraph -c bbb.bedgraph -m subtract -o ccc.bedgraph
	grep 'chr' ccc.bedgraph > ddd.bedgraph
	awk '{if($4>=0){print $0}}' ddd.bedgraph | tr " " "\t" > eee
	awk '{if($4<0){print $1,$2,$3,"0"}}' ddd.bedgraph | tr " " "\t" > fff
	cat eee fff > ggg
	bedtools sort -i ggg > $(echo ${NAMES[$i]}).minus_Halo.bdg
	rm aaa*
	rm bbb*
	rm ccc*
	rm ddd*
	rm eee
  rm fff
	rm ggg

done
rm -rf NAMES

## Substracting the input
NAMES=($(ls *minus_Halo.bdg | grep -v 'Input' | sed 's/\./ /g' | awk '{print $1}'))
for i in $(seq 0 $(ls *minus_Halo.bdg | grep -v 'Input' | wc -l | awk '{print $1-1}'))


do


        cat $(echo ${NAMES[$i]})*minus_Halo.bdg > aaa.bedgraph
        cat Input*minus_Halo.bdg > bbb.bedgraph

        macs2 bdgcmp -t aaa.bedgraph -c bbb.bedgraph -m subtract -o ccc.bedgraph
        grep 'chr' ccc.bedgraph > ddd.bedgraph
        awk '{if($4>=0){print $0}}' ddd.bedgraph | tr " " "\t" > eee
        awk '{if($4<0){print $1,$2,$3,"0"}}' ddd.bedgraph | tr " " "\t" > fff   # If the substraction results in a negative number (Halo more signal than the protein then I print 0)
        cat eee fff > ggg
        bedtools sort -i ggg > $(echo ${NAMES[$i]}).minus_Halo.minus_input.bdg
        rm aaa*
        rm bbb*
        rm ccc*
        rm ddd*
        rm eee
	rm fff
	rm ggg

done
rm -rf NAMES

mkdir MINUS_HALO_INPUT
mv *minus_Halo.minus_input.bdg MINUS_HALO_INPUT

cd MINUS_HALO_INPUT


module load bioinfo-tools
module load MACS/2.2.6
module load BEDTools



NAMES=($(ls *.minus_input.bdg | sed 's/\./ /g' | awk '{print $1}'))
for i in $(seq 0 $(ls *.minus_input.bdg | wc -l | awk '{print $1-1}'))


do
	
	grep 'chr' $(echo ${NAMES[$i]})*.minus_input.bdg > fixed.bed
	bedtools sort -i fixed.bed > aaa.bdg
	macs2 bdgpeakcall -i aaa.bdg -c 0.5 -l 20 -g 5 -o $(echo ${NAMES[$i]}).peaks.bed
	rm fixed.bed
	rm aaa.bdg

done
rm -rf NAMES

mkdir PEAKS
mv *.peaks.bed ./PEAKS
cd ./PEAKS





NAMES=($(ls *peaks.bed | sed 's/\./ /g' | awk '{print $1}'))
for i in $(seq 0 $(ls *peaks.bed | wc -l | awk '{print $1-1}'))

do

	awk '{print $1,$2,$3,$2+$10}' $(echo ${NAMES[$i]})*peaks.bed | grep 'chr' | tr " " "\t" > $(echo ${NAMES[$i]}).peaks.fixed.bed

done
rm -rf NAMES


mkdir FIXED_BEDs
mv *fixed.bed ./FIXED_BEDs
cd ./FIXED_BEDs





NAMES=($(ls *fixed.bed | sed 's/\_/ /g' | awk '{print $1}' | sort -u))
for i in $(seq 0 $(ls *rep1*fixed.bed | wc -l | awk '{print $1-1}'))

do

        bedtools intersect -a $(echo ${NAMES[$i]})*rep1*fixed.bed -b $(echo ${NAMES[$i]})*rep2*fixed.bed -wb | awk '{print $1,$2,$3,$4,$8,$3-$2,$8-$4}' > aaa.txt
	awk '{if($6>=20){print $0}}' aaa.txt | awk '{if($7<11 && $7>-11){print $0}}' | awk '{if($4>$2 && $4<$3){print $0}}' | awk '{if($5>$2 && $5<$3){print $0}}' | tr " " "\t" > $(echo ${NAMES[$i]}).common.peaks.bed
	rm aaa.txt	
done

mkdir COMMON
mv *common.peaks.bed ./COMMON

