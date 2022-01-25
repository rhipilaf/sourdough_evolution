# Méthodes

## Phénotypage au robot

Comme dans Bigey et al. (2021)

Mesure au cours du temps (entre 26 et 28 heures) toutes les 51 minutes, de la masse d'un milieu synthétique contenant *Saccharomyces cerevisiae* qui y fermente. La transformation des sucres en CO<sub>2</sub> entraine une perte de poids qui est un proxy de la production de CO<sub>2</sub> cumulée à chaque temps $t$ (avec $t_{max} \in [32,34]$).

Barboteur sur le dessus, donc pas de contamination possible, et pas d'entrée de matière, uniquement une sortie. Un tube ne peut que perdre de la matière.

Milieu liquide, composition cf. Bigey et al. (2021) <u>et Diego pour confirmation de la composition</u>.

A $t_0$ et à $t_{27}$, comptage au cytomètre de flux le nombre de cellules. L'effectif de la population à $t_{27}$ constitue un proxy de la fitness.

### Limites du protocole

- N'y a t il pas une production d'autres gazs (*e.g.* O<sub>2</sub>) qui entrainent une perte de masse mais qui n'est pas liée à la production de CO<sub>2</sub> ? Est-ce que l'hypothèse forte consistant à admettre que l'essentiel des gazs produits est du CO<sub>2</sub> est valable ?

- Milieu différent du levain réel, mais comme le montre Bigey et al (2021) *(et d'autres papiers ? Michel et al. 2022 ?)*, bien que synthétique, ce milieu peut mettre en évidence des phénotypes variables entre souches, notamment entre souches commerciales et sauvages.

- **Mesure de fitness** : c'est une approche prospective car on mesure des traits d'histoire de vie des souches isolées sans tenir compte d'une potentielle compétition directe entre souches. Une hypothèse forte permettant d'inférer une fitness à partir de ces mesures est que la compétition entre souches n'est qu'indirecte par une efficacité variable dans la consommation des ressources.

### Analyse

* Moyenne glissante et différence glissante

* Ajustement d'un modèle de croissance
  
  * Moindres carrés
  
  * Maximum likelihood
    
    * nlmer
  
  * Bayésien
    
    * Tuto : [Fitting differential equation models with Stan](https://shug3502.github.io/blog/DifferentialEqnsStan)
    * Sensitivity analysis : [adjustr](https://corymccartan.github.io/adjustr/), [rgiordan/StanSensitivity](https://github.com/rgiordan/StanSensitivity), 
    * [Hierarchical Gompertz model avec *brms*](https://discourse.mc-stan.org/t/hierarchical-gompertz-model/13724/13), [Stan coding question for Gompertz curve fitting](https://discourse.mc-stan.org/t/stan-coding-question-for-gompertz-curve-fitting/14618), [Estimating Non-Linear Models with brms](https://cran.r-project.org/web/packages/brms/vignettes/brms_nonlinear.html)
