#!/bin/bash

all_bins=$(ls ARS_bins_v2/bin_summary/bin_by_bin/*/*-gene_calls.txt)
out="bin_COG_diversity.csv"
rm $out
touch $out

get_shannon_weaver_index()
{	
	file=$1
	bin_name=$(echo $file | awk -F'/' '{print $5}' | sed 's/-gene_calls.txt//g')
	tmp_var=$(awk -F'\t' '/meta/{ print $11}' $file | sed 's/!!!/\n/g' | awk '/[PCGEFHI]/' | sort | uniq -c | awk '{ print $1 }')
	total=0
	for each in $tmp_var; do total=$(( $total + $each )); done
	shannon_weaver_index=0
	for each in $tmp_var; do
		prop=$(echo "$each / $total" | bc -l)
		logNum=$(echo "l($prop)" | bc -l)
		shannon_weaver_index=$(echo "$shannon_weaver_index - $logNum * $prop" | bc -l)
	done
	echo $bin_name,$shannon_weaver_index >> $out
}

for file in $all_bins; do 
	echo doing $file 
	get_shannon_weaver_index $file;
done


# One-letter abbreviations for the functional categories: 

# P, inorganic ion transport and metabolism; 
# C, energy production and conversion; 
# G, carbohydrate metabolism and transport; 
# E, amino acid metabolism and transport; 
# F, nucleotide metabolism and transport; 
# H, coenzyme metabolism; 
# I, lipid metabolism; 

# J, translation, including ribosome structure and biogenesis; 
# L, replication, recombination and repair; 
# K, transcription; 
# O, molecular chaperones and related functions; 
# M, cell wall structure and biogenesis and outer membrane; 
# N, secretion, motility and chemotaxis; 
# T, signal transduction; 
# P, inorganic ion transport and metabolism; 
# C, energy production and conversion; 
# G, carbohydrate metabolism and transport; 
# E, amino acid metabolism and transport; 
# F, nucleotide metabolism and transport; 
# H, coenzyme metabolism; 
# I, lipid metabolism; 
# D, cell division and chromosome partitioning; 
# R, general functional prediction only; 
# S, no functional prediction.




