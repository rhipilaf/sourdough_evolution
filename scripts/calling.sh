#!/bin/bash

# ARGS
# First argument ($1) must be the number of threads to use on the processor
# Second argument ($2) must be the path to the reference genome
# Third argument ($3) must be the path to bam files (input folder)
# Fourth argument ($4) must be the path to the output folder

export readonly path_bam=$3
export readonly path_output=$4
export readonly gatk=/opt/gatk-4.2.5.0/gatk
export readonly samples_bam=$(ls $path_bam | grep ".bam$" | awk -F ".bam" '{print $1}' - | uniq) #All samples that have .bam file in the mapping folder
export readonly samples_done=$(ls $path_output | grep ".g.vcf$" | awk -F ".g.vcf" '{print $1}' - | uniq) #All samples that already have .g.vcf file in the calling folder
export readonly samples_todo=$(comm -23 <(printf "%s\n" $samples_bam | sort) <(printf "%s\n" $samples_done | sort) | tr '\n' ' ') #Samples that only are in $samples_bam, not in $samples_done
export readonly threads=$1 # Number of threads to use simultaneously
export readonly refgenome=$2
export readonly ploidy=2


function calling {

	local s=$1
	$gatk --java-options "-Xmx16g" HaplotypeCaller -R $refgenome -I $path_bam/$s.bam -O $path_output/$s.g.vcf --sample-ploidy $ploidy -ERC GVCF
}


# Runs SNP calling parallely with the number of thread specified (first argument at the script calling)
echo $samples_todo | xargs -d' ' -P$threads -I % -n1 bash -c "$(declare -f calling) ; echo 'SNP calling for %' ; calling % 1>$path_output/%.g.vcf.log 2>&1" {}


# Notes :
# - the $(declare -f calling) and the 'export readonly ...'s in this script are needed to allow the access of the calling function to variables when called with xargs (details here : https://confluence.acfr.usyd.edu.au/pages/viewpage.action?pageId=51839827)
# - in xargs, -d' ' is needed to cut elements of the list retreived by echo.
# - in xargs, -n1 is for the number of elements of the list to be passed to the call as arguments.
# - in xargs, -I % specifies the character (%) to be use as the position where piped argument has to be inserted in the call.
# - at the end of the xargs line, don't know what {} is useful for.
