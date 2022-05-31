#!/bin/bash

# Argument list
# 1) => g.vcf cohorte (input)
# 2) => output base name
# 3) => reference genomes

# Default path/values
gatk=/opt/gatk-4.2.5.0/gatk

# First genotype the g.vcf
$gatk --java-options "-Xmx16g" GenotypeGVCFs -R $3 -V $1 -O tmp01.vcf

# Extract biallelic SNP
$gatk --java-options "-Xmx16g" SelectVariants -R $3 -V tmp01.vcf -O tmp02.vcf --select-type-to-include SNP --restrict-alleles-to BIALLELIC

# Clean-up
rm -f tmp01.vcf*

# Apply hard-filtering
$gatk --java-options "-Xmx16g" VariantFiltration -R $3 -V tmp02.vcf -O tmp03.vcf -filter "QD < 2.0" --filter-name Low_QD -filter "FS > 60.0" --filter-name High_FS -filter "SOR > 3.0" --filter-name High_SOR -filter "MQ < 40.0" --filter-name Low_MQ -filter "MQRankSum < -2.5 || MQRankSum > 2.5" --filter-name "Bad_MQRS" -filter "ReadPosRankSum < -8.0 || ReadPosRankSum > 8.0" --filter-name Bad_RPRS

# Clean-up
rm -f tmp02.vcf*

# Trim filtered
$gatk --java-options "-Xmx16g" SelectVariants -R $3 -V tmp03.vcf -O $2.vcf --exclude-filtered

# Clean-up
rm -f tmp03.vcf*

# Extract table
$gatk --java-options "-Xmx16g" VariantsToTable -R $3 -V $2.vcf -O $2.table -F CHROM -F POS -F REF -F ALT -GF AD

