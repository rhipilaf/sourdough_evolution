#!/bin/bash

# ARGS
# $1 = path and name of the reference genome
# $2 = a text file containing a list of *.g.vcf files. The file can be created via piped (|) ls.
# $3 = path and name of the output file

refgenome=$1
list_gvcf=$2
path_output=$3
format_output=".g.vcf"
log_file=$(echo $path_output | awk -F $format_output '{print $1}' - )_$(date '+%Y-%m-%d_%H-%M-%S')$format_output.log
gatk=/opt/gatk-4.2.5.0/gatk


touch $log_file
echo $(date '+%Y-%m-%d_%H-%M-%S') >> $log_file
echo "" >> $log_file
cat $list_gvcf >> $log_file
echo "" >> $log_file

$gatk CombineGVCFs -R $refgenome --variant $list_gvcf -O $path_output 1>>$log_file 2>&1
