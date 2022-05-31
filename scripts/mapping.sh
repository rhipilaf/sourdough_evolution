#!/bin/bash

# ARGS:
# First argument ($1) must be the number of threads to use on the processor
# Second argument ($2) must be the path to the reference genome
# Third argument ($3) must be the path to trimmed .fq.gz files (input folder)
# Fourth argument ($4) must be the path of the output folder

export readonly path_trim=$3
export readonly path_output=$4
export readonly samples_fqgz=$(ls $path_trim | grep ".fq.gz$" | awk -F "_[1-2].fq.gz" '{print $1}' - | uniq) #All samples that have .fq.gz file in the trimming folder
export readonly samples_done=$(ls $path_output | grep ".bam$" | awk -F ".bam" '{print $1}' - | uniq) #All samples that already have a .bam file in the mapping folder
export readonly samples_todo=$(comm -23 <(printf "%s\n" $samples_fqgz | sort) <(printf "%s\n" $samples_done | sort) | tr '\n' ' ')
export readonly threads=$1 # Number of threads to use simultaneously
export readonly refgenome=$2
export readonly samtools=/opt/samtools-1.10/samtools

#TO WRITE PROPERLY : if sample_todo == NULL, echo "All fastq files were trimmed. So nothing to do for me. If you want to rerun mapping, remove corresponding .fq.gz files in the trimming folder." ; break

echo "SAMPLES TO DO:"
echo $samples_todo
echo "---"

# Function that does the mapping with the sample name as only argument
mapping(){


  echo $path_output "in"
  local s=$1 # The argument of the function must be the id of the currently analyzed sample

  # Mise en oeuvre du mapping. Creation du fichier .sam
  bwa mem -t 12 -M -R "@RG\tID:${s}\tSM:${s}" -o $path_output/$s.sam $refgenome $path_trim/${s}_1.fq.gz $path_trim/${s}_2.fq.gz

  # Calcul d'une stat du fichier sam
  $samtools stat $path_output/$s.sam > $path_output/$s.stat

  # Suite d'instruction pour nettoyer, ordonner et compresser le fichier sam en bam
  $samtools fixmate -r -m -@ 6 $path_output/$s.sam - | samtools sort -@ 12 - | samtools markdup -r -@ 12 -O BAM - $path_output/$s.bam

  # Indexer le fichier bam
  $samtools index $path_output/$s.bam

  # Supprimer le sam (si on est sûr)
  rm $path_output/$s.sam ;

}


# Run mapping for each sample name in parallel according to the number of threads defines at the script calling (first argument)
echo $samples_todo | xargs -d' ' -P$threads -I % -n1 bash -c "$(declare -f mapping) ; echo 'Genome mapping for %' ; mapping % 1>$path_output/%.mapping.log 2>&1" {}

echo "-- JOB COMPLETED --"
echo ""

