#!/bin/bash

for i in `seq 001 999`; do 
	var=$(printf "%03d\n" $i)
	downlink="ftp://ftp.sra.ebi.ac.uk/vol1/ERZ843/ERZ843$var/"
	result=$(curl $downlink| grep -e "bam$")
	errnum=$(echo $result | awk -F " " '/ERR/ { print $9 }' | sed s/.bam//)
	cat filter_big_download.csv | grep $errnum
	if [[ $? == 0 ]]
	then
		echo "I found $errnum"
		line=$(cat filter_big_download.csv | grep $errnum)
		newname=$(echo $line | awk -F "," '{ print $5 }' | tr -d '\r')
		echo ${downlink}
		wget ${downlink}${errnum}.bam
		wget ${downlink}${errnum}.bam.bai
		mv ${errnum}.bam ${newname}.bam
		mv ${errnum}.bam.bai ${newname}.bam.bai
	fi
	echo $result
	echo $errnum
	echo $var
done

