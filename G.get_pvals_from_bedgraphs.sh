#!/bin/bash
#Made by Ionut Atanasoai
#Made on 200525
#UPPMAX commands
#SBATCH -A snic2020-15-293
#SBATCH -p core
#SBATCH -n 2




module load bioinfo-tools
module load MACS/2.2.6
module load BEDTools/2.29.2


## Substracting the input
NAMES=($(ls *.bam.scaled.bedgraph | grep -v 'Input' | sed 's/\./ /g' | awk '{print $1}'))
for i in $(seq 0 $(ls *.bam.scaled.bedgraph | grep -v 'Input' | wc -l | awk '{print $1-1}'))


do


        cat $(echo ${NAMES[$i]})*.bam.scaled.bedgraph > aaa.bedgraph
        cat Input*.bam.scaled.bedgraph > bbb.bedgraph

        macs2 bdgcmp -t aaa.bedgraph -c bbb.bedgraph -m subtract -o ccc.bedgraph
        grep 'chr' ccc.bedgraph > ddd.bedgraph
        awk '{if($4>=0){print $0}}' ddd.bedgraph | tr " " "\t" > eee
        awk '{if($4<0){print $1,$2,$3,"0"}}' ddd.bedgraph | tr " " "\t" > fff   # If the substraction results in a negative number (Input more signal than the protein then print 0)
        cat eee fff > ggg
        bedtools sort -i ggg > $(echo ${NAMES[$i]}).minus_input.bdg
        rm aaa*
        rm bbb*
        rm ccc*
        rm ddd*
        rm eee
	rm fff
	rm ggg

done
rm -rf NAMES

mkdir MINUS_INPUT
mv *minus_input.bdg MINUS_INPUT

cd MINUS_INPUT/


## Computing Pvalues in -log10 scale
NAMES=($(ls *.bdg | grep -v 'Halo' | sed 's/\./ /g' | awk '{print $1}'))
for i in $(seq 0 $(ls *.bdg | grep -v 'Halo' | wc -l | awk '{print $1-1}'))


do


	tail -n+2 $(echo ${NAMES[$i]})*bdg > aaa.bedgraph
	tail -n+2 Halo*bdg > bbb.bedgraph

	macs2 bdgcmp -t aaa.bedgraph -c bbb.bedgraph -p 1 -m ppois -o ccc.bedgraph
	grep 'chr' ccc.bedgraph > ddd.bedgraph
	bedtools sort -i ddd.bedgraph > $(echo ${NAMES[$i]}).pvalues_from_minus_Inputs.bdg
	rm aaa*
	rm bbb*
	rm ccc*
	rm ddd*


done


mkdir PVALS_BIS
mv *pvalues_from_minus_Inputs.bdg PVALS_BIS/
