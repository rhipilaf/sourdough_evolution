#!/bin/bash

path_sra=/media/data/Bakery/Master_Theo/other_strains/sra_files
path_output=/media/data/Bakery/Master_Theo/other_strains/fq_files
fasterqdump=/opt/sratoolkit.3.0.0-ubuntu64/bin/fastq-dump
log_file=$path_output/sra_to_fq_$(date '+%Y-%m-%d_%H-%M-%S').log
touch $log_file

for s in $(ls $path_sra | awk -F '.sra' '{print $1}' -  )
do
    sra=${s}.sra
    fq1=${s}_1.fastq
    fq2=${s}_2.fastq

    if [[ -f $path_output/$fq1 && -f $path_output/$fq2 ]]; then
        echo "$s already dumped into fastq."
        continue
    fi

    $fasterqdump --split-files -O $path_output $path_sra/$sra 1>>$log_file 2>&1

done
