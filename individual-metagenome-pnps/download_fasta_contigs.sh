#!/bin/bash

for i in `seq 001 999`; do 
	var=$(printf "%03d\n" $i)
	echo $var
	downlink="ftp://ftp.sra.ebi.ac.uk/vol1/ERZ842/ERZ842$var/"
	result=$(curl $downlink| grep ".fasta.gz")
	echo $result
	errnum=$(echo $result | awk -F " " '/ERR/ { print $9 }' | sed s/.fasta.gz//)
	echo $errnum
	cat new_download.csv | grep $errnum
	echo test
	if [[ $? == 0 ]]
	then
		echo "I found $errnum"
		line=$(cat new_download.csv | grep $errnum)
		newname=$(echo $line | awk -F "," '{ print $5 }' | tr -d '\r')
		echo ${downlink}
		wget ${downlink}${errnum}.fasta.gz
		gunzip ${errnum}.fasta.gz
		mv ${errnum}.fasta ${newname}.fasta
		anvi-gen-contigs-database --contigs-fasta ${newname}.fasta -T 10 --output ${newname}.db
	fi
done