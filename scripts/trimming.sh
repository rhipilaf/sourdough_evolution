#!/bin/bash

path_input="/media/data/Bakery/Master_Theo/other_strains/fq_files"
path_output="00_trimming"
samples=$(ls -R $path_input/ | awk -F "_" '{print $1}' - | uniq)
fastp="/usr/local/bin/fastp"
encoding="phred33"
force="n" #To force the trimming to be done, even if the output file already exists
ext="fastq" #can be either fastq or fq.gz

$fastp -v


if [[ "$encoding" == "phred64" ]]; then
	phred="--phred64"
else
	phred=""
fi


for s in $samples; do

	fqi1=${s}_1.$ext
	fqi2=${s}_2.$ext
	fqf1=${s}_1.$ext
	fqf2=${s}_2.$ext

	if [[ -f $path_output/${s}_1.fq.gz && -f $path_output/${s}_2.fq.gz ]]; then

		if [[ "$force" != "y" ]]; then
			echo "$path_output/${s}_*.fastq already created from $path_input/${s}_*.fq.gz." 
			continue

		else
			echo "IMPORTANT: $path_output/${s}_*.fq.gz already exist and is being overwritten."

		fi

	fi

	echo "Processing and writing $path_input/${s}_*.fq..."


	$fastp -l 50 -w 8 $phred -i $path_input/$fqi1 -I $path_input/$fqi2 -o $path_output/$fqf1 -O $path_output/$fqf2 -h $path_output/$s.html -j $path_output/$s.json

	if [[ "$ext" == "fastq" ]]; then
		gzip -c $path_output/$fqf1 > $path_output/${s}_1.fq.gz
		gzip -c $path_output/$fqf2 > $path_output/${s}_2.fq.gz

		rm $path_output/$fqf1
		rm $path_output/$fqf2
	fi

done
