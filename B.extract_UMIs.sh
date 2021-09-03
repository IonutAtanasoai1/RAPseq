#!/bin/bash
#Made by Ionut Atanasoai
#Made on 200513
#UPPMAX commands
#SBATCH -A snic2020-15-293 
#SBATCH -p core
#SBATCH -n 4


 

NAMES=($(ls *read1*.fastq | sed 's/\.read1/\t/g' | awk '{print $1}' | sort -u))


mkdir UMIED

for i in $(seq 0 $(ls *read1*.fastq | wc -l | awk '{print $1-1}'))
do

	cat $(echo ${NAMES[$i]})*read1*.fastq | paste - - - - > ALL_R1.txt
	cat $(echo ${NAMES[$i]})*read2*.fastq |	paste -	- - - >	ALL_R2.txt
	
	awk '{print $3}' ALL_R1.txt | sed 's/....$//g' | sed 's/^....//g' > fixed.3.R1.txt
	awk '{print $5}' ALL_R1.txt | sed 's/....$//g' | sed 's/^....//g' > fixed.5.R1.txt

    awk '{print $3}' ALL_R2.txt | sed 's/....$//g' | sed 's/^....//g' > fixed.3.R2.txt
    awk '{print $5}' ALL_R2.txt | sed 's/....$//g' | sed 's/^....//g' > fixed.5.R2.txt	

	awk '{print $3}' ALL_R1.txt | grep -o '^....' > umi.R1.txt
	awk '{print $3}' ALL_R2.txt | grep -o '^....' >	umi.R2.txt
	paste -d '' umi.R1.txt umi.R2.txt > UMI.txt

	awk '{print $1}' ALL_R1.txt > ONE.R1.txt
	paste -d '_' ONE.R1.txt UMI.txt > FIELD_1.R1.txt

	awk '{print $1}' ALL_R2.txt > ONE.R2.txt
	paste -d '_' ONE.R2.txt	UMI.txt	> FIELD_1.R2.txt
	
	awk '{print $2}' ALL_R1.txt > TWO.R1.txt
	awk '{print $4}' ALL_R1.txt > FOUR.R1.txt
        awk '{print $2}' ALL_R2.txt > TWO.R2.txt
        awk '{print $4}' ALL_R2.txt > FOUR.R2.txt	

	paste FIELD_1.R1.txt TWO.R1.txt | awk '{print $1,$2}' > FIELDS_1_2.R1.txt
	paste FIELDS_1_2.R1.txt fixed.3.R1.txt FOUR.R1.txt fixed.5.R1.txt | tr "\t" "\n" > ./UMIED/$(echo ${NAMES[$i]}).R1.trimmed_umied.fastq

        paste FIELD_1.R2.txt TWO.R2.txt | awk '{print $1,$2}' > FIELDS_1_2.R2.txt
        paste FIELDS_1_2.R2.txt fixed.3.R2.txt FOUR.R2.txt fixed.5.R2.txt | tr "\t" "\n" > ./UMIED/$(echo ${NAMES[$i]}).R2.trimmed_umied.fastq		

	rm *.txt

done

