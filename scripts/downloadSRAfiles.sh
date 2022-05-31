#!/bin/bash

path_output=/media/data/Bakery/Master_Theo/other_strains/sra_files
log_file=$path_output/sra_dl_$(date '+%Y-%m-%d_%H-%M-%S').log
prefetch=/opt/sratoolkit.3.0.0-ubuntu64/bin/prefetch
touch $log_file

for s in $(cut -f8 99_metadata/other_strains.tsv | tail -n +2)
do
    sra=${s}.sra

    if test -f "$path_output/$sra"; then
        echo "$sra already downloaded."
        continue
    fi

    $prefetch -o $path_output/$sra $s 1>>$log_file 2>&1

done

