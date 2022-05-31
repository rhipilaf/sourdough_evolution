#!/bin/bash

path_bam=02_mapping/
path_depth=02_mapping/
samtools=/opt/samtools-1.10/samtools

for b in $(ls -p $path_bam*.bam | awk -F '.bam' '{print $1}' - | awk -F '/' '{print $2}' -  )
do
	bam=$b.bam
	depth=$b.bam.depth

	 if [[ -f $path_depth$depth ]]; then
		echo "$path_depth$depth already created from $path_bam$bam."
		continue
	fi

	echo "Proccessing and writing $path_depth$depth..."

	$samtools depth -aa $path_bam$bam > $path_depth$depth

done
