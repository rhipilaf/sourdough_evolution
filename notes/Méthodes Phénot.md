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
