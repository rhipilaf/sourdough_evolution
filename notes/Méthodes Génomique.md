# Données de séquences

## Répertoires de stockage des données brutes

```bash
/media/data/Bakery/These_Lucas/BGI_JAN2022/ #Souches Lucas
/media/data/Bakery/These_Lauriane/BGI_JAN2022/ #Souches Lauriane
/media/data/Bakery/Master_Theo/other_strains/fq_files/ #Souches en plus
/media/data/Bakery/BGI_MAR2022/ #Souches Lauriane, Lucas et commerciales
```

### Files list generation

```bash
 for f in $(cat 99_metadata/data_repos.txt); do ls $f*_[1-2].f*; done | grep "^/*" 1> 99_metadata/raw_files.txt
```

## Processing

Une ligne de fichier fastq consiste en un *read* et des adaptateurs.

Chaque read apparatient à une paire qui ont été lu sur le même brin. Dans les deux fichiers .fastq, les paires de reads sont dans le même ordre.

### Fonction pour avoir la longueur d'une sequence

`sequenceLength.pl` directement accessible sur le serveur

Génère un tableau nécessaire à mettre dans l'argument 2 du script `runGATKtoMatrix.R`

### Couverture

Taille génome S. cerev. = 12 Mb

**Longueur read x 2 x nb_paires_read / taille génome = couverture pour chaque base.**

Pour des individus, une couverture de 100 suffit. Pour des pops, une couverture de 200 c'est pas mal. En dessous de 50, c'est chaud :/

Avec des reads de 150b, la probabilité de chevauchement entre reads et plus probable qu'en 75

### Encoding detection

```bash
vsearch="/media/workspace/guillert/sdgh_evol/98_softwares/vsearch-2.21.1-linux-x86_64/bin/vsearch"
for f in $(cat 99_metadata/raw_files.txt | grep "_1.f")
do
sample=$(echo $f | awk -F "/" '{print $NF}' - | awk -F "." '{print $1}' - )
$vsearch --fastq_chars $f > 99_metadata/encoding/$sample.chars
done
```

```bash
touch 99_metadata/encoding/phred_summary.txt
for f in $(ls 99_metadata/encoding/*.chars)
do
sample=$(echo $f | awk -F "/" '{print $NF}' - | awk -F "." '{print $1}' - )
echo $sample " " $(grep "Guess:" $f | grep "(phred") >> 99_metadata/encoding/phred_summary.txt
done
```

BGI = phred64

Peter et al = phred33

Evolya = phred33

Bigey = phred33

Gallone = phred33

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

### Mapping/alignement des reads (`bwa mem`)

Ne pas générer tous les fichiers `.sam` un par un car ils sont très gros. Il faut le nettoyer avant de lancer un second.

### Nettoyage des reads mappés

`samtools stat`

`samtools fixmate`

`-r` = supprime les flags 4 et 256
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

### Extraction des SNPs

#### Combiner les GVCFs

Permet d'avoir l'ensemble des souches dans un seul fichier avec l'ensemble des SNPs étant variable au moins pour une des souches :

```bash
ls 03_calling/*.g.vcf > gvcf.list
nohup ./combining.sh 01_refgenome/S288c_ABC.fasta gvcf.list 04_SNPfiltering/allGenotypes_combined.g.vcf 1>04_SNPfiltering/combining_$(date '+%Y-%m-%d_%H-%M-%S').log 2>&1 &
```

#### Récupérer les SNP bialléliques

```bash
nohup ./buildCohorteBialSNP.sh 04_SNPfiltering/allGenotypes_combined.g.vcf 04_SNPfiltering/allGenotypes_filtered 01_refgenome/S288c_ABC.fasta 1>>04_SNPfiltering/buildCohorteBialSNP.sh.log 2>&1 &
```

Le hard filtering est paramétré selon les recommandations

#### Extraction des fréquences alléliques

Extraires les données de fréquences alléliques

```bash
sequenceLength.pl -i 01_refgenome/S288c_ABC.fasta > 01_refgenome/S288c_ABC.length
./runGATKtoMatrix.R 04_SNPfiltering/allGenotypes_filtered.table 01_refgenome/S288c_ABC.length 05_foranalysis/allGenotypes_AlleleFreq.csv
```

#### Extraire la couverture

```bash
nohup ./gettingdepth.sh 1>04_SNPfiltering/combining_$(date '+%Y-%m-%d_%H-%M-%S').log 2>&1 &
```

##### Compresser la couverture

Fenetre de 1000, pas de 500

```bash
for f in $(ls 02_mapping/*.bam.depth); do ./depth_compression.R $f 06_depth; done
```

Warning :

```R
Processing 02_mapping/ERR1308783.bam.depth -DONE
Warning message: In scan(file = file, what = what, sep = sep, quote = quote, dec = dec,  :
le nombre d'objets lus n'est pas un multiple du nombre de colonnes
```

#### Extraire les séquences consensus

```bash
gatk=/opt/gatk-4.2.5.0/gatk
$gatk FastaAlternateReferenceMaker -R 01_refgenome/S288c_ABC.fasta -O 05_foranalysis/allGenotypes_consensus.fasta -V 04_SNPfiltering/allGenotypes_filtered.vcf
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

Et maintenant il faut relancer le SNP calling sur le nouveau fichier bam recalibré `*_recal.bam`.

## Analyse des données génomiques

### Distance génétique

Dans Gutaker et al 2019, = 1- IBS (identity by state ; [Identity by type - Wikipedia](https://en.wikipedia.org/wiki/Identity_by_type))

### Fréquences alléliques

Supprimer les lignes non informatives qui ont les même freq allélique

Suppr les lignes ayant au moins un -1 sauf peut-être pour l'analyse de la couverture

#### Plot fréquence allélique vs position dans le génome

Mettre en évidence les limites des chromosomes avec des lignes verticales

### Couverture

`samtools depth -aa fichier.bam > fichier.bam.depth`

`-aa` = toutes les positions même celles qui ne sont pas couvertes (affiche les 0)

nombre de reads moyens le long du génome : fenetre de 1000 avec un pas de 500
Fonction de fenetre glissante peut être codée en C (utilisable sous R une fois compilé)
`mcapply` ? = multicore apply
Profondeur en log2 ?
Afficher les bordures de chromosomes !
Voir ce qu'il se passe au niveau des régions ABC qui ont une couverture très faible

### D de Tajima

Doit être fait sur une population

*"**Tajima’s D**. We calculated nucleotide diversity (mean of pairwise sequence differences within population) and number of segregating sites for each population in each photoperiod response gene (Supplementary Methods 2.18). We then calculated Tajima’s D by applying equations implemented in `tajima.test` function of the `pegas` package to quad-allelic individuals. Finally, we visualized the empirical distribution for Tajima’s D using the `yarrr` package in `R` v.3.4.1 (ref. 55) and highlighted genes in the 99th percentile of this empirical distribution."* (Gutaker et al. 2019)

## Collecte données NCBI

Répertoire des données : `/media/data/Bakery/Master_Theo`

### Téléchargement des fichiers `.sra`

`downloadingSRAfiles.sh` script (*à fignoler!*)

```bash
#!/bin/bash

path_ouput=/media/data/Bakery/Master_Theo/other_strains/sra_files
log_file=$path_ouput/sra_dl_$(date '+%Y-%m-%d_%H-%M-%S').log
prefetch=/opt/sratoolkit.3.0.0-ubuntu64/bin/prefetch
touch $log_file

for s in $(cut -f8 99_metadata/other_strains.tsv | tail -n +2)
do
    sra=${x}.sra

    if test -f "$path_ouput/$sra"; then
        echo "$sra already downloaded."
        continue
    fi

    $prefetch -o $path_output/$sra $s 1>>$log_file 2>&1
done
```

To backgroundly start the downloading

```bash
nohup ./downloadSRAfiles.sh 1>>downloadSRAfiles.sh.log 2>&1 &
```

### Extraction des fichiers` .fastq` à partir d'un `.sra`

`SRAtoFq.sh` (a rendre plus polyvalent, car cette version prend tout ce qu'il y a dans le répertoire à sra et )

```bash
#!/bin/bash

path_sra=/media/data/Bakery/Master_Theo/other_strains/sra_files
path_output=/media/data/Bakery/Master_Theo/other_strains/fq_files
fasterqdump=/opt/sratoolkit.3.0.0-ubuntu64/bin/fastq-dump
log_file=$path_ouput/sra_to_fq_$(date '+%Y-%m-%d_%H-%M-%S').log
touch $log_file

for s in $(ls $path_sra)
do
    sra=${s}.sra
    fq1=${s}_1.fastq
    fq2=${s}_2.fastq

    if [test -f "$path_ouput/$fq1"] && [test -f "$path_ouput/$fq2"]; then
        echo "$s already dumped into fastq."
        continue
    fi

    $fasterqdump --split-files -O $path_output $path_sra/$sra 1>>$log_file 2>&1

done
```

To backgroundly start the dumping :

```bash
nohup ./SRAtoFQ.sh 1>>SRAtoFQ.sh.log 2>&1 &
```

To compare files sizes between origin and destination :

```bash
source=/media/data/Bakery/Master_Theo/other_strains/sra_files
target=/media/workspace/guillert/sdgh_evol


for i in "$source"/*.sra
do
 f1=`stat -c%s $i`
 f2=`stat -c%s $target/${i##*/}`
  if [ "$f1" = "$f2" ]; then
        echo "$i" "$f1" VS "$target/${i##*/}" "$f2" "====>>>" "OK"
  else
        echo "$i" "$f1" VS "$target/${i##*/}" "$f2" "====>>>" "BAD"
  fi
done
```

```bash
grep "====>>> BAD" Zmachin.txt | awk -F '/' '{print $8}' - | awk -F ' ' '{print $1}' - > failed_sra.list
```

pangenome = union

core genome = intersection

regarder les stats de reads non mappés pour voir si certains échantillons

S288C + ABC = chimérique => S288C le plus complet / ABC (issues de EC1118) régions d'intérêt pour la fermentation (régions ABC poster à côté du bureau de Virginie)

```bash
nohup ./jobs.sh 1>jobs_sh_$(date '+%Y-%m-%d_%H-%M-%S').log 2>&1 &
```

## Divers

### Pour faire un assemblage rapide

```bash
/opt/SPAdes-3.15.3-Linux/bin/spades.py --pe-1 /media/data/Bakery/BGI_25   
```

$$
V_t = \frac{p_{t-1} + p_t + p_{t+1}}{xxxxxx} \approx \frac{d{p_t}}{dt}
$$

V_t = \frac{p_{t-1} + p_t + p_{t+1}}{xxxxxx} \approx \frac{d{p_t}}{dt}

$$
V_t = \frac{\sum_{t = -h}^{h}(p_t)}{\sum_{t = -h}^{h}(t)}
$$

### Renommage des fichiers de Fred Bigey

```bash
for f in $(find /media/data/Bakery/papier_Bigey2021/ | grep fastq); do 
sample=$(echo $f | awk -F "OSW" '{print $1}' - | awk -F "_" '{print $3}' - ) 
b=$(echo $f | awk -F "_" '{print $5}' - )
ln -sT $f /media/data/Bakery/Master_Theo/papier_Bigey2021_symlinks/${sample}_${b}.fq.gz 
done
```

(Jombart, 2009) = review de méthodo

PCoA et pas PCA (Higgins 1992)

Regroupement hiérarchique

# 
