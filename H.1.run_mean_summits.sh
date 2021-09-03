#bash
#Made by Sofia Papavasileiou
NAMES=($(ls *bed | sed 's/\./ /g' | awk '{print $1}' | sort -u))
FILES=($(ls *.bed))

for i in $(seq 0 $(ls *.bed | sed 's/\./ /g' | awk '{print $1}' | sort -u | wc -l | awk '{print $1-1}')) 
do
	Rscript mean_summits.R  $(echo ${FILES[$i]})  table.txt	
	sed 's/"//g' table.txt >  $(echo ${NAMES[$i]})_9nt_window.common.bed
	rm table.txt
done
