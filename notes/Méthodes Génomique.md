Données de séquences

Une ligne de fichier fastq consiste en un *read* et des adaptateurs.

Chaque read apparatient à une paire qui ont été lu sur le même brin. Dans les deux fichiers .fastq, les paires de reads sont dans le même ordre.

### Couverture

Taille génome S. cerev. = 12 Mb

**Longueur read x 2 x nb_paires_read / taille génome = couverture pour chaque base.**

Pour des individus, une couverture de 100 suffit. Pour des pops, une couverture de 200 c'est pas mal. En dessous de 50, c'est chaud :/

Avec des reads de 150b, la probabilité de chevauchement entre reads et plus probable qu'en 75

### Controle qualité

FastQC reports

### Nettoyage des données

```bash
nohup ./trimming.sh 1>./00_trimming/trimming_$(date '+%Y-%m-%d_%H-%M-%S').log 2>&1 &
```

### Indexation du génome de référence (S288C)

`bwa` = donne les fichiers `.amb`, `.ann` , `.bwt`, `.pac`, `.sa`

`gatk` = donne un fichier `.dict`

`samtools` = donne un fichier `.fai`

### Mapping (`bwa mem`)

Ne pas générer tous les fichiers `.sam` un par un car ils sont très gros. Il faut le nettoyer avant de lancer un second.

### Nettoyage des reads mappés

`samtools stat`

`samtools flimate`

`-r` = supprime les flags 4 et 256 ()
`-m` = calcule un score d'appariement
`-@` # = nombre de threads de mémoire utilisés. Ne pas trop en utiliser car cela peut prendre plus de temps

`samtools markdup`

`-r` = remove duplicated

```bash
nohup ./mapping.sh 1 01_refgenome/S288c_ABC.fasta 00_trimming 02_mapping 1>02_mapping/mapping_$(date '+%Y-%m-%d_%H-%M-%S').log 2>&1 & 
```

### SNP calling avec `gatk`

`-ERC GVCF` = histoire d'optimisation. Redemander à Hugo à quoi ça sert

```bash
nohup ./calling.sh 3 01_refgenome/S288c_ABC.fasta 02_mapping 03_calling 1>03_calling/calling_$(date '+%Y-%m-%d_%H-%M-%S').log 2>&1 &
```

### For PoolSeq analysis

#### SNP filtering avec `gatk`

On parle ici de *hard filtering method* de `gatk`

Sélection des SNP pour lesquels on est sûr qu'ils existent

```bash
/opt/gatk-4.2.5.0/gatk --java-options "-Xmx16g" SelectVariants -R 01_refgenome/S288c_ABC.fasta -V 03_callpool/2-EC1118_4n.vcf -O 03_callpool/2-EC1118_4n_rawSNP.vcf --select-type-to-include SNP
```

##### `gatk` SNP statistics

###### Extraction des statistiques

```bash
/opt/gatk-4.2.5.0/gatk --java-options "-Xmx16g" VariantsToTable -R 01_refgenome/S288c_ABC.fasta -V 03_callpool/2-EC1118_4n_rawSNP.vcf -O 03_callpool/2-EC1118_4n_rawSNP_stat.txt -F CHROM -F POS -F REF -F ALT -F AF -F QUAL -F DP -F QD -F FS -F SOR -F MQ -F MQRankSum -F ReadPosRankSum 
```

###### Visualisation des SNPs filtrés selon les critères suivants

```
QD < 25
Test exact de Fischer : FS > 10
SOR > 2
MQ < 50
MQRankSum < -2.5 ou > 2.5
ReadPosRankSum < -3 ou > 3
```

```bash
runPlotGATKStat.R 03_callpool/2-EC1118_4n_rawSNP_stat.txt ouput.pdf
```

##### Ajout des flags nécessaire au filtre

```bash
/opt/gatk-4.2.5.0/gatk --java-options "-Xmx16g"  VariantFiltration -R 01_refgenome/S288c_ABC.fasta -V 03_callpool/2-EC1118_4n_rawSNP.vcf  -O 03_callpool/2-EC1118_4n_HQSNP_tagged.vcf -filter "QD < 25.0" --filter-name Low_QD -filter "FS > 10.0" --filter-name High_FS -filter "SOR > 2.0" --filter-name High_SOR -filter "MQ < 50.0" --filter-name Low_MQ -filter "MQRankSum < -2.5 || MQRankSum > 2.5" --filter-name "Bad_MQRS" -filter "ReadPosRankSum < -3.0 || ReadPosRankSum > 3.0" --filter-name Bad_RPRS
```

`gatk [...] VariantFiltration [...]` functioning :

- Spécifier les règles d'exclusion et non d'inclusion !
  
- Les valeurs doivent être des chiffres à virgule
  
- Lors du filtrage, il ajoute en flag le nom du filtre si le SNP correspond à la règle de filtre. Si elle ne correspond pas, le flag devient `PASS`
  

##### Effective filtering

`gatk [...] SelectVariants [...]` gives a file with the SNPs that only have the SNP with the `PASS` flag

```bash
/opt/gatk-4.2.5.0/gatk --java-options "-Xmx16g"  SelectVariants -R 01_refgenome/S288c_ABC.fasta -V 03_callpool/2-EC1118_4n_HQSNP_tagged.vcf -O 03_callpool/2-EC1118_4n_HQSNP.vcf --exclude-filtered
```

#### Recalibration pour refaire un SNP calling sachant les SNPs de haute qualité

L'idée est de donner du poids à de potentiels SNPs ayant le même profil que les SNPs sélectionnés plus haut. Mais le process utilisé reste un mystère.

##### Génération d'un fichier de recalibration

- le fichier bam utilisé est le fichier .bam de base lors du premier SNP calling fait avec `gatk`.
  
- `BQSR_table.tmp` est le fichier de recalibration qui est temporaire. Dans le cas ou je parallélise, il faut bien que ce soit un nom propre à l'échantillon.
  

```bash
/opt/gatk-4.2.5.0/gatk --java-options "-Xmx16g" BaseRecalibrator -R 01_refgenome/S288c_ABC.fasta -I /media/workspace/cbecerra/ALE2022/02_mappool/fixed_platform/2-EC1118.bam -O BQSR_table.tmp --known-sites 03_callpool/2-EC1118_4n_HQSNP.vcf 
```

##### Effective recalibration using the calibration file

```bash
/opt/gatk-4.2.5.0/gatk --java-options "-Xmx16g" ApplyBQSR -R 01_refgenome/S288c_ABC.fasta -I /media/workspace/cbecerra/ALE2022/02_mappool/fixed_platform/2-EC1118.bam -O 02_mappool/2-EC1118_recal_4n.bam --bqsr-recal-file BQSR_table.tmp
```

Et maintenant il faut relancer le SNP calling sur le nouveau fichier bam recalibré `*_recal.bam`

### Analyse des données génomiques

#### Distance génétique

Dans Gutaker et al 2019, = 1- IBS (identity by state ; [Identity by type - Wikipedia](https://en.wikipedia.org/wiki/Identity_by_type))