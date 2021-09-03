#!/bin/bash
#Made by Ionut Atanasoai
#Made on 210102
#UPPMAX commands
#SBATCH -A snic2020-15-293
#SBATCH -p core
#SBATCH -n 2



NAMES=($(ls *.R1.*fastq | sed 's/\.R1/\t/g' | awk '{print $1}' | sort -u))

Read1=($(ls *.R1.*fastq))


mkdir SUBSETTED


for i in $(seq 0 $(ls *.R1.*fastq | wc -l | awk '{print $1-1}'))


do

	cat $(echo ${NAMES[$i]})*.R1.*fastq | paste - - - - | sort -k1,1 > AAA.txt
	cat ../$(echo ${NAMES[$i]})*.R2.*fastq | paste - - -	- | sort -k1,1 > BBB.txt
	join -1 1 -2 1 AAA.txt BBB.txt > CCC.txt

	awk '{print $1,$2"\n"$3"\n"$4"\n"$5}' CCC.txt > ./SUBSETTED/$(echo ${NAMES[$i]}).R1.for_genome_alignment.fastq
	awk '{print $1,$6"\n"$7"\n"$8"\n"$9}' CCC.txt > ./SUBSETTED/$(echo ${NAMES[$i]}).R2.for_genome_alignment.fastq

	rm *.txt

done

