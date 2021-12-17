#!/bin/bash

transeed="../BLAST_tara_OM-RGC_v2/OM-RGC_v2_Reference_Pieces/Transposase_Database.faa"
firstarg=$1

echo "search in Tara${firstarg}..."

for db in $(ls Tara${firstarg}*/*.db); do
	basename=${db%.db}
	mv ${basename}.fasta.n* ../trash/. 
	mv ${basename}_transposase.blastp ../trash/.
	fastaname=${basename}_anvi-genecall.fasta
	outname=${basename}_anvi-genecall_transposase.blastp
	foldername=$(echo $basename | awk -F'/' '{print $1}' )
	echo will move $foldername to f_$foldername

	anvi-get-sequences-for-gene-calls -c $db -o $fastaname
	makeblastdb -in $fastaname -dbtype nucl
	tblastn -query $transeed -db $fastaname -outfmt 6 -evalue 1e-05 -out $outname -num_threads 5
	mv $foldername f_${foldername}
	echo "done with $basename"
done
