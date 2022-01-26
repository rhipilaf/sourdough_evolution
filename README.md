# Evolution de populations de *Saccharomyces cerevisiae* au cours de la colonisation du levain

## Structure du répertoire

- Dans `rapports/` : petits rapports en .markdown pour nos réunions

- Dans `scripts/` : ensemble des scripts utilisés pour les analyses

- Dans `notes/` : j'y stocke ma prise de notes



| parameter  | fit_name                                 | logLik  | AIC    | deltaAIC |
| ----------:|:---------------------------------------- | -------:| ------:| --------:|
| cell_t27   | **~ sn + ft + wt + bs + bk + Rbm + Rc0** | -2460.2 | 5050.4 | 0.0      |
| cell_t27   | ~ ft + wt + bs + bk + Rbm + Rc0 + Rsn    | -3377.5 | 6779.1 | 1728.7   |
| cell_t27   | ~ ft + wt + bs + bk + Rbm + Rc0          | -3379.9 | 6781.7 | 1731.3   |
| cell_t27   | ~ sn + ft + wt + bs + bk                 | -3457.5 | 7041.0 | 1990.6   |
| co2max     | **~ ft + wt + bs + bk + Rbm + Rc0**      | -141.0  | 304.0  | 0.0      |
| co2max     | ~ ft + wt + bs + bk + Rbm + Rc0 + Rsn    | -141.0  | 306.0  | 2.0      |
| co2max     | ~ sn + ft + wt + bs + bk                 | -100.7  | 327.4  | 23.3     |
| co2max     | ~ sn + ft + wt + bs + bk + Rbm + Rc0     | -130.7  | 391.4  | 87.3     |
| death_prct | **~ sn + ft + wt + bs + bk + Rbm + Rc0** | -440.1  | 1010.3 | 0.0      |
| death_prct | ~ ft + wt + bs + bk + Rbm + Rc0          | -577.7  | 1177.3 | 167.0    |
| death_prct | ~ ft + wt + bs + bk + Rbm + Rc0 + Rsn    | -577.3  | 1178.5 | 168.2    |
| death_prct | ~ sn + ft + wt + bs + bk                 | -546.0  | 1217.9 | 207.6    |
| lag        | **~ ft + wt + bs + bk + Rbm + Rc0**      | -265.5  | 553.0  | 0.0      |
| lag        | ~ ft + wt + bs + bk + Rbm + Rc0 + Rsn    | -265.5  | 555.0  | 2.0      |
| lag        | ~ sn + ft + wt + bs + bk + Rbm + Rc0     | -223.4  | 576.8  | 23.8     |
| lag        | ~ sn + ft + wt + bs + bk                 | -240.1  | 606.3  | 53.3     |
| tvmax      | **~ sn + ft + wt + bs + bk + Rbm + Rc0** | -207.9  | 545.9  | 0.0      |
| tvmax      | ~ sn + ft + wt + bs + bk                 | -210.9  | 547.7  | 1.9      |
| tvmax      | ~ ft + wt + bs + bk + Rbm + Rc0 + Rsn    | -300.4  | 624.8  | 78.9     |
| tvmax      | ~ ft + wt + bs + bk + Rbm + Rc0          | -320.7  | 663.5  | 117.6    |
| vmax       | **~ sn + ft + wt + bs + bk**             | 59.1    | 7.8    | 0.0      |
| vmax       | ~ ft + wt + bs + bk + Rbm + Rc0 + Rsn    | -9.0    | 42.0   | 34.1     |
| vmax       | ~ ft + wt + bs + bk + Rbm + Rc0          | -13.4   | 48.8   | 40.9     |
| vmax       | ~ sn + ft + wt + bs + bk + Rbm + Rc0     | -13.2   | 156.4  | 148.5    |

> 
