#!/bin/bash

# ARGS
# First argument ($1) must be the number of threads to use on the processor
# Second argument ($2) must be the path to the reference genome
# Third argument ($3) must be the path to trimmed .fq.gz files (input folder)
# Fourth argument ($4) must be the path of the output folder

path_trim=$3
path_output=$4
samples_fqgz=$(ls $path_trim | grep ".fq.gz$" | awk -F ".fq.gz" '{print $1}' - | uniq) #All samples that have .fq.gz file in the trimming folder
#samples_done=$(ls $path_output | grep ".bam$" | awk -F ".bam" '{print $1}' - | uniq) #All samples that already have a .bam file in the mapping folder
#samples_todo=$(setdiff between samples_bam & samples_done)
threads=$1 # Number of threads to use simultaneously
nb_samples=${#samples_bam[@]}
refgenome=$2
samtools=/opt/samtools-1.10/samtools

#TO WRITE PROPERLY : if sample_todo == NULL, echo "All fastq files were trimmed. So nothing to do for me. If you want to rerun mapping, remove corresponding .fq.gz files in the trimming folder." ; break


function mapping {
  
  local s=$1 # The argument of the function must be the id of the currently analyzed sample
  
  # for s in $($data_path | awk -F "_" '{print $1}' - | uniq;
  bwa mem -t 8 -M -R "@RG\tID:${s}\tSM:${s}" -o $path_output/$s.sam $refgenome $path_trim/${s}_1.fq.gz $path_trim/${s}_2.fq.gz

  # Calcul d'une stat du fichier sam
  $samtools stat $path_output/$s.sam > $path_output/$s.stat

  # Suite d'instruction pour nettoyer, ordonner et compresser le fichier sam en bam
  $samtools fixmate -r -m -@ 6 $path_output/$s.sam - | samtools sort -@ 6 - | samtools markdup -r -@ 6 -O BAM - $path_output/$s.bam

  # Indexer le fichier bam
  $samtools index $path_output/$s.bam

  # Supprimer le sam (si on est sûr)
  rm $path_output/$s.sam ;
  
}


echo $samples_fqgz | xargs -d' ' -P$threads -I % -n1 bash -c 'echo "Genome mapping for %.." ; echo % 1>'$path_output'/%.mapping.log 2>&1 ; sleep 2' {}
echo "-- JOB COMPLETED --"
echo ""

