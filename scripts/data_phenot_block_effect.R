
# BLOCK EFFECT

library(lme4)

traits <- c("cell_t27", "co2max", "death_prct", "lag", "tvmax", "vmax")


# Correction by the block effect

fits <- data_phenot_parms_clean %>%
  filter(parameter != "cell_t0") %>%
  mutate(flour_id = strains$flour_id[match(strain_name, strains$strain_name)],
         baker_id = strains$baker_id[match(strain_name, strains$strain_name)],
         flour_type = flours$mill_type[match(flour_id, flours$flour_id)],
         wheat_type = flours$wheat_type[match(flour_id, flours$flour_id)],
         backslopping = strains$backslopping[match(strain_name, strains$strain_name)],
         cell_t0 = data_cyto$cell_t0[match(robot_id, data_cyto$robot_id)]) %>%
  filter(!if_any(everything(), ~ is.na(.x))) %>%
  group_by(parameter) %>%
  nest() %>%
  mutate(fit_sn_ft_wt_bs_bk_Rbm_Rc0 = map(data, ~ lmer(value ~ strain_name + flour_type * wheat_type + backslopping + baker_id + (1|bloc_month) + (1|cell_t0), data = .)),
         fit_sn_Rbm_Rc0 = map(data, ~ lmer(value ~ strain_name + (1|bloc_month) + (1|cell_t0), data = .)),
         fit_sn = map(data, ~ lm(value ~ strain_name, data = .)),
         fit_Rsn_Rbm_Rc0 = map(data, ~ lmer(value ~ 1 + (1|strain_name) + (1|bloc_month) + (1|cell_t0), data = .)),
         fit_ft_wt_bs_bk_Rbm_Rc0 =  map(data, ~ lmer(value ~ flour_type * wheat_type + backslopping + baker_id + (1|bloc_month) + (1|cell_t0), data = .)),
         fit_sn_ft_wt_bs_bk = map(data, ~ lm(value ~ strain_name + flour_type * wheat_type + backslopping + baker_id, data = .)),
         fit_ft_wt_bs_bk_Rbm_Rc0_Rsn = map(data, ~ lmer(value ~ flour_type * wheat_type + backslopping + baker_id + (1|bloc_month) + (1|cell_t0) + (1|strain_name), data = .))) %>%
  pivot_longer(cols = starts_with("fit_"), names_to = "fit_name", values_to = "fit")


infos <- fits %>%
  mutate(summary = map(fit, glance)) %>%
  unnest(summary) %>%
  select(-data, -fit) %>%
  arrange(parameter, AIC) %>%
  group_by(parameter) %>%
  mutate(deltaAIC = AIC - min(AIC)) %>%
  select(parameter, fit_name, logLik, AIC, deltaAIC, sigma, REMLcrit, df.residual) %>%
  mutate_if(is.numeric, ~ round(.x, 1)) %>%
  mutate(fit_name = str_replace_all(fit_name, "^fit_", "~ "),
         fit_name = str_replace_all(fit_name, "_", " + "),
         fit_name = ifelse(deltaAIC < 2, str_replace_all(fit_name, "^", "**"), fit_name),
         fit_name = ifelse(deltaAIC < 2, str_replace_all(fit_name, "$", "**"), fit_name)) %>%
  knitr::kable(format = "pipe", align = "rlrrrrrr")

res_fitted_plots <- fits %>%
  mutate(summary = map(fit, augment),
         AIC = map_dbl(fit, ~ unlist(glance(.x))['AIC'])) %>%
  arrange(parameter, AIC) %>%
  unnest(summary) %>%
  select(-data, -fit) %>%
  group_by(parameter) %>%
  nest() %>%
  mutate(data = map(data, ~ .x %>% mutate(fit_name = fct_reorder(fit_name, AIC))), 
         plots = map2(data, parameter,
                    ~ ggplot(data = .x) +
                      aes(x = .fitted, y = .resid) +
                      geom_point() +
                      facet_grid(. ~ fit_name, scales = "free") +
                      theme_minimal() +
                      ggtitle(.y) +
                      theme(strip.text.y = element_text(angle = 0),
                            title = element_text(face = "bold")))) %>%
  select(plots, parameter) %>%
  as.list()

gridExtra::arrangeGrob(grobs = res_fitted_plots$plots, ncol = 1, nrow = length(unique(fits$parameter)), top = "") %>%
  myggsave(., filename = "output/res_fitted_plots", device = "pdf", width = 15, height = 20)


## 1. Estimation des effets blocs avec un modele par milieu
# On estime l'effet bloc uniquement avec le Sourdough.

LMERS0 = lmer(logtrait ~ Species + (1|Strain) +(1|Block), data=mydataS, REML=TRUE)


## 2. Modèle complet
# On travaille sur l'echelle log (résidus mieux pour certaines variables + 
# interprétation est la même pour tous : on estime un rapport). On estime 
# l'effet bloc pour le comparer à celui qu'on a estimé en 1.

LMER0 = lmer(logtrait ~Species * Habitat+(1|Strain)+(1|Block),data=mydataV,REML=TRUE)


## 3. Estimation de l'effet bloc en effet fixe
# Avec seulement les souches présentes dans les deux blocs
# On teste l'interaction, puis on estime l'effet dans le modèle additif si pas d'interaction
#on travaille sur la variable en log et on calcul des effets fixes Souche et block

LM0= lm(logtrait ~ Strain*Block,data=mydataV)
print(anova(LM0))

AnoInter<-rbind.data.frame(AnoInter,cbind.data.frame(trait=trait,pInter=anova(LM0)$P[3]))

LM0= lm(logtrait ~ Strain+Block,data=mydataV)


## 3'. Modèle complet sur echelle log
# on retire au préalable l'effet bloc
# On estime l'effet milieu pour chaque espèce.
# On stocke estimation et IC




strains_distrib_inblocs <- data_cyto %>%
  select(strain_name, bloc) %>%
  unique() %>%
  group_by(strain_name) %>%
  summarise(nb_blocs = n()) %>%
  arrange(nb_blocs)

data_phenot_parms_clean %>%
  mutate(parameter = add_units(parameter)) %>%
  ggplot(.) +
  aes(y = value, fill = bloc) +
  geom_boxplot(alpha = 0.3) +
  facet_wrap(~parameter, scales = "free") +
  theme_minimal()

bloc_effect_model <- data_phenot_parms_clean %>%
  group_by(parameter) %>%
  summarise(p.value = anova(lm(value ~ bloc))$'Pr(>F)'[1])





## Data structure : distribution of strains among blocks ####
links <- table(data_cyto$strain_name, data_cyto$bloc) %>%
  as.data.frame() %>%
  filter(Freq != 0)

nodes <- data.frame(name=c(as.character(links$Var1),
                           as.character(links$Var2)
) %>% unique()
)

links$IDsource <- match(links$Var1, nodes$name)-1 
links$IDtarget <- match(links$Var2, nodes$name)-1

p <- networkD3::sankeyNetwork(Links = links, Nodes = nodes,
                              Source = "IDsource", Target = "IDtarget",
                              Value = "Freq", NodeID = "name", 
                              sinksRight=FALSE, iterations = 100,
                              width = 500, height = 1500)
htmlwidgets::saveWidget(p, file=paste0(getwd(), "/output/sankeyBlocks.html"))
