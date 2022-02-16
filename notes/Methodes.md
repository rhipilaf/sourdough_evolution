# Méthodes

## Echantillonnage des levains

### Codage

- **TEMOIN** = levain du boulanger (*i.e.* levures déjà adaptées ?)

- **MAISON** = levain expérimental réalisé avec la farine utilisée habituellement par le boulanger (*i.e.* farine à laquelle les levures indigènes sont déjà adaptées?)

## Phénotypage au robot

Comme dans Bigey et al. (2021)

Mesure au cours du temps (entre 26 et 28 heures) toutes les 51 minutes, de la masse d'un milieu synthétique contenant *Saccharomyces cerevisiae* qui y fermente. La transformation des sucres en CO<sub>2</sub> entraine une perte de poids qui est un proxy de la production de CO<sub>2</sub> cumulée à chaque temps $t$ (avec $t_{max} \in [32,34]$).

Barboteur sur le dessus, donc pas de contamination possible, et pas d'entrée de matière, uniquement une sortie. Un tube ne peut que perdre de la matière.

Milieu liquide, composition cf. Bigey et al. (2021) <u>et Diego pour confirmation de la composition</u>.

A $t_0$ et à $t_{27}$, comptage au cytomètre de flux le nombre de cellules. L'effectif de la population à $t_{27}$ constitue un proxy de la fitness.

### Recette du milieu (levain synthétique)

= W-SSM (Wheat Sourdough Simulation Medium) modifié, d’après Vrancken, 2008

- Composition (par L) :
  
  - Wheat peptone : 24 g
  
  - MgSO4.7H2O : 0.2 g
  
  - MnSO4.H2O : 0.05 g
  
  - KH2PO4 : 4 g
  
  - K2HPO4 : 4 g
  
  - Tween 80 : 1 ml
  
  - Glucose : 15 g
  
  - Maltose : 35 g
  
  - Solution de vitamines 1000x : 1 ml
    
    - Composition de la solution de vitamines 1000x (par L) :
      
      - Cobalamine : 0.2 g
      
      - Acide folique : 0.2 g
      
      - Nicotinamide : 0.2 g
      
      - Acide pantothénique : 0.2 g
      
      - Pyridoxal-phosphate : 0.2 g
      
      - Thiamine : 0.2 g
      
      - La solution est ensuite filtrée sur 0.22µm puis stockée à -20°C

- pH ajusté à 4.5 avec de l’acide citrique. 60.7g dans 100ml d’eau distillée donne une solution 10x (environ 3M), stockée à -20°C

- Stérilisation 15 min en vapeur fluente

### Limites du protocole

- N'y a t il pas une production d'autres gazs (*e.g.* O<sub>2</sub>) qui entrainent une perte de masse mais qui n'est pas liée à la production de CO<sub>2</sub> ? Est-ce que l'hypothèse forte consistant à admettre que l'essentiel des gazs produits est du CO<sub>2</sub> est valable ?

- Milieu différent du levain réel, mais comme le montre Bigey et al (2021) *(et d'autres papiers ? Michel et al. 2022 ?)*, bien que synthétique, ce milieu peut mettre en évidence des phénotypes variables entre souches, notamment entre souches commerciales et sauvages.

- **Mesure de fitness** : c'est une approche prospective car on mesure des traits d'histoire de vie des souches isolées sans tenir compte d'une potentielle compétition directe entre souches. Une hypothèse forte permettant d'inférer une fitness à partir de ces mesures est que la compétition entre souches n'est qu'indirecte par une efficacité variable dans la consommation des ressources.

### Analyse

#### Extraction des paramètres

Les variables phénotypiques associées à chaque souche sont des statistiques décrivant à la fois :

##### La démographie de la population finale de levures

Ces mesures sont directement mesurées au du cytomètre de flux.

- $death_\%$ = pourcentage de cellules morte à $t_f$

- $pop_{size}$ = population vivante à $t_f$

##### La dynamique de fermentation alcoolique des sucres par la population de levures

Cette dynamique est mesurée par perte de poids du système (voir paragraphe 'phénotypage au robot'), proxy direct de la quantié de CO<sub>2</sub> produite par les levures : 

$$
C_6H_{12}O_6 \rightarrow 2 \cdot C_2H_5OH + 2 \cdot CO_2
$$

Cette dynamique peut-être décrite par des statistiques extraites à partir des données brutes de perte de poids ($m_t$) produites à l'aide du robot directement transformées en une valeur production cumulée de CO<sub>2</sub> $p_t$ en grammes.

$$
p_t = \frac{m_{t=2} - m_t}{volume} \times 1000
$$

et un débit de CO<sub>2</sub> $V_t$ en g/h habituellement mesuré par une pente locale mesurée sur trois points.

$$
V_t = \frac{p_{t-1} + p_t + p_{t+1}}{} \approx \frac{d{p_t}}{dt}
$$

###### Estimation des paramètres classique (méthode utilisée par le labo)

A partir de ces valeurs ont été mesurées 4 statistiques descriptive de la dynamique fermentaire : latence $\lambda$ (en h), CO<sub>2</sub> produit cumulé final $p_{max}$ (en g), débit maximal de CO<sub>2</sub> $V_{max}$(en g/h) et le temps où le débit maximal $t_{V_{max}}$.

- $\lambda = min(t) \thickspace \forall p_t > 1$

- $p_{max} = max(p_t)$

- $V_{max} = max(V_t)$

- $t_{V_{max}}=min(t)\thickspace \forall V_t = V_{max}$

###### Estimation des paramètre par ajustement d'un modèle

La dynamique de production de CO<sub>2</sub> par fermentation peut aussi être modélisée par un modèle de croissance de Gompertz modifié (formule ci-dessous) comme proposée par Mohammadi *et al.* (2014) pour la production d'éthanol, qui est stoechiométriquement égale à la production de CO<sub>2</sub>. 

$$
p_t=p_{max} \cdot \exp \bigg( -\exp \bigg[ \frac{V_{max}\times e}{p_{max}}(\lambda-t)+1\bigg]\bigg)
$$

L'ajustement de ce modèle de croissance à la production de CO<sub>2</sub> cumulée $p_t$ a été réalisé pour chacun des réplicats par la méthode des moindres carrés avec départs de chaines multiples avec le la fonction R `nls.multstart::nls_multstart()`.

Le $t_{V_{max}}$ a été calculé à posteriori de l'estimation des autres paramètres est cherchant le $t$ pour lequel la dérivée partielle en $t$ de la courbe ajustée aux données est maximale. Cette dérivée partielle s'exprime par la formule suivante :

$$

$$

En perspective, j'aimerai réaliser cet ajustement dans un cadre bayésien ou fréquentiste. Quelques pistes pour faire ça en bayésien : [Fitting differential equation models with Stan](https://shug3502.github.io/blog/DifferentialEqnsStan), [Sensitivity analysis : adjustr](https://corymccartan.github.io/adjustr/), [rgiordan/StanSensitivity](https://github.com/rgiordan/StanSensitivity), [Hierarchical Gompertz model avec *brms*](https://discourse.mc-stan.org/t/hierarchical-gompertz-model/13724/13), [Stan coding question for Gompertz curve fitting](https://discourse.mc-stan.org/t/stan-coding-question-for-gompertz-curve-fitting/14618), [Estimating Non-Linear Models with brms](https://cran.r-project.org/web/packages/brms/vignettes/brms_nonlinear.html).

#### Analyse de variance des paramètres

- Inoculum size as a random effect : (Davis et al., 2005; Liu et al., 2021)
* Unbalanced incomplete block designs
  
  * [Testing of hypotheses in unbalanced incomplete block designs | Request PDF](https://www.researchgate.net/publication/288692473_Testing_of_hypotheses_in_unbalanced_incomplete_block_designs)
  
  * [Chapter 8 Incomplete Block Designs | ANOVA and Mixed Models](https://stat.ethz.ch/~meier/teaching/anova/incomplete-block-designs.html)
  
  * [4.7 - Incomplete Block Designs | STAT 503](https://online.stat.psu.edu/stat503/lesson/4/4.7)
  
  * https://stats.oarc.ucla.edu/stata/faq/how-can-i-analyze-an-unbalanced-randomized-block-design/
  
  * https://ecommons.cornell.edu/bitstream/handle/1813/32993/BU-885-M.pdf;jsessionid=7F71B991E1D0A5C20EA531BA019B7818?sequence=1
  
  * [A thorough randomization in R of a Partially Balanced Incomplete Block Design.](http://rstudio-pubs-static.s3.amazonaws.com/312_5cf6573d98d247d0a23c1d36ed7e7485.html)

* AIC based comparison
  
  * Not a problem if it is not nested models. But this nestedness is a condition for LRT. [r - Comparing non nested models with AIC - Cross Validated](https://stats.stackexchange.com/questions/116935/comparing-non-nested-models-with-aic)

* Rank deficiency
  
  * [r - What is rank deficiency, and how to deal with it? - Cross Validated](https://stats.stackexchange.com/questions/35071/what-is-rank-deficiency-and-how-to-deal-with-it)

## Données de séquences

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

### Analyse de structure de population
