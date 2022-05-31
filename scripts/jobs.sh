#!/bin/bash

#./trimming.sh

#./mapping.sh 1 01_refgenome/S288c_ABC.fasta 00_trimming 02_mapping

#./calling.sh 3 01_refgenome/S288c_ABC.fasta 02_mapping 03_calling



## Building the all-strains cohort

#./combining.sh 01_refgenome/S288c_ABC.fasta 99_metadata/cohort_lists/gvcfs_all.list 04_SNPfiltering/allGenotypes_combined.g.vcf

#./buildCohorteBialSNP.sh 04_SNPfiltering/allGenotypes_combined.g.vcf 04_SNPfiltering/allGenotypes_filtered 01_refgenome/S288c_ABC.fasta

#sequenceLength.pl -i 01_refgenome/S288c_ABC.fasta > 01_refgenome/S288c_ABC.length
#./runGATKtoMatrix.R 04_SNPfiltering/allGenotypes_filtered.table 01_refgenome/S288c_ABC.length 05_foranalysis/allGenotypes_AlleleFreq.csv



## Building bakery strains cohort

#./combining.sh 01_refgenome/S288c_ABC.fasta 99_metadata/cohort_lists/gvcfs_bak.list 04_SNPfiltering/bakGenotypes_combined.g.vcf

#./buildCohorteBialSNP.sh 04_SNPfiltering/bakGenotypes_combined.g.vcf 04_SNPfiltering/bakGenotypes_filtered 01_refgenome/S288c_ABC.fasta

#sequenceLength.pl -i 01_refgenome/S288c_ABC.fasta > 01_refgenome/S288c_ABC.length
#./runGATKtoMatrix.R 04_SNPfiltering/bakGenotypes_filtered.table 01_refgenome/S288c_ABC.length 05_foranalysis/bakGenotypes_AlleleFreq.csv



## Building

./combining.sh 01_refgenome/S288c_ABC.fasta 99_metadata/cohort_lists/gvcfs_sub.list 04_SNPfiltering/subGenotypes_combined.g.vcf

./buildCohorteBialSNP.sh 04_SNPfiltering/subGenotypes_combined.g.vcf 04_SNPfiltering/subGenotypes_filtered 01_refgenome/S288c_ABC.fasta

sequenceLength.pl -i 01_refgenome/S288c_ABC.fasta > 01_refgenome/S288c_ABC.length
./runGATKtoMatrix.R 04_SNPfiltering/subGenotypes_filtered.table 01_refgenome/S288c_ABC.length 05_foranalysis/subGenotypes_AlleleFreq.csv



#./gettingdepth.sh

#for f in $(ls 02_mapping/*.bam.depth); do ./depth_compression.R $f 05_foranalysis; done
