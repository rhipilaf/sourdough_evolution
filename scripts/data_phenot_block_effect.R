
# BLOCK EFFECT

traits <- c("cell_t27", "co2max", "death_prct", "lag", "tvmax", "vmax")


# Classic analysis ####

## Fitting ####

fits_classic <- data_phenot_parms_clean %>%
  filter(parameter != "cell_t0") %>%
  mutate(flour_id = strains$flour_id[match(strain_name, strains$strain_name)],
         baker_id = strains$baker_id[match(strain_name, strains$strain_name)],
         flour_type = flours$mill_type[match(flour_id, flours$flour_id)],
         wheat_type = flours$wheat_type[match(flour_id, flours$flour_id)],
         backslopping = strains$backslopping[match(strain_name, 
                                                   strains$strain_name)],
         cell_t0 = data_cyto$cell_t0[match(robot_id, data_cyto$robot_id)]) %>%
  # filter(!if_any(everything(), ~ is.na(.x))) %>%
  group_by(parameter) %>%
  nest() %>%
  mutate(fit_sn_ft_wt_bs_bk_Rbm_Rc0 = map(data, ~ lmer(value ~ strain_name 
                                                       + flour_type * wheat_type 
                                                       + backslopping + baker_id 
                                                       + (1|bloc_month) 
                                                       + (1|cell_t0), 
                                                       data = .)),
         fit_sn_ft_wt_bs_bk = map(data, ~ lm(value ~ strain_name + flour_type 
                                             * wheat_type + backslopping + 
                                               baker_id, data = .)),
         fit_sn_Rbm_Rc0 = map(data, ~ lmer(value ~ strain_name + (1|bloc_month) 
                                           + (1|cell_t0), data = .)),
         fit_sn = map(data, ~ lm(value ~ strain_name, data = .)),
         fit_Rsn_Rbm_Rc0 = map(data, ~ lmer(value ~ 1 + (1|strain_name) 
                                            + (1|bloc_month) + (1|cell_t0), 
                                            data = .)),
         fit_ft_wt_bs_bk_Rbm_Rc0 =  map(data, ~ lmer(value ~ flour_type 
                                                     * wheat_type + backslopping 
                                                     + baker_id + (1|bloc_month) 
                                                     + (1|cell_t0), data = .)),
         fit_ft_wt_bk_Rbm_Rc0 =  map(data, ~ lmer(value ~ flour_type 
                                                  * wheat_type + baker_id 
                                                  + (1|bloc_month) 
                                                  + (1|cell_t0), data = .)),
         fit_ft_wt_bs_Rbm_Rc0 =  map(data, ~ lmer(value ~ flour_type 
                                                  * wheat_type + backslopping 
                                                  + (1|bloc_month) 
                                                  + (1|cell_t0), data = .)),
         fit_bk_Rbm_Rc0 =  map(data, ~ lmer(value ~ baker_id
                                            + (1|bloc_month) 
                                            + (1|cell_t0), data = .)),
         fit_bk_Rbm_Rc0_Rsn =  map(data, ~ lmer(value ~ baker_id
                                                + (1|bloc_month) 
                                                + (1|cell_t0)
                                                + (1|strain_name), data = .)),
         fit_ft_wt_bs_bk_Rbm_Rc0_Rsn = map(data, ~ lmer(value ~ flour_type 
                                                        * wheat_type 
                                                        + backslopping 
                                                        + baker_id 
                                                        + (1|bloc_month) 
                                                        + (1|cell_t0) 
                                                        + (1|strain_name), 
                                                        data = .))) %>%
  
  pivot_longer(cols = starts_with("fit_"), names_to = "fit_name", 
               values_to = "fit")

## Info ####

fits_classic_infos <- fits_classic %>%
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

## Res Fit plot ####

res_fit_plot <- fits_classic %>%
  mutate(summary = map(fit, augment),
         AIC = map_dbl(fit, ~ unlist(glance(.x))['AIC'])) %>%
  unnest(summary) %>%
  select(-data, -fit) %>%
  group_by(parameter) %>%
  nest() %>%
  mutate(data = map(data, ~ .x %>% mutate(fit_name = fct_reorder(fit_name, AIC))), 
         res_fitted_plot = map2(data, parameter,
                      ~ ggplot(data = .x) +
                        aes(x = .fitted, y = .resid) +
                        geom_point() +
                        facet_grid(. ~ fit_name, scales = "free") +
                        theme_minimal() +
                        ggtitle(.y) +
                        theme(strip.text.y = element_text(angle = 0),
                              title = element_text(face = "bold")))) %>%
  ungroup() %>%
  select(res_fitted_plot) %>%
  as.list() %>%
  gridExtra::arrangeGrob(grobs = .$res_fitted_plot, ncol = 1, nrow = length(unique(fits$parameter)), top = "") %>%
  myggsave(., filename = "output/data_phenot_classic_res_fitted_plots", device = "pdf", width = 15, height = 20)



# Correction by the block effect ####

## Block effects estimation ####

strains_in_both <- data_phenot_parms_clean %>%
  select(strain_name, bloc_month) %>%
  unique() %>%
  group_by(strain_name) %>%
  summarise(n = n()) %>%
  filter(n > 1) %>%
  select(strain_name) %>% 
  flatten_chr()

fits_blockeff_estim <- data_phenot_parms_clean %>%
  filter(parameter != "cell_t0",
         strain_name %in% strains_in_both) %>%
  mutate(flour_id = strains$flour_id[match(strain_name, strains$strain_name)],
         baker_id = strains$baker_id[match(strain_name, strains$strain_name)],
         flour_type = flours$mill_type[match(flour_id, flours$flour_id)],
         wheat_type = flours$wheat_type[match(flour_id, flours$flour_id)],
         backslopping = strains$backslopping[match(strain_name, strains$strain_name)],
         cell_t0 = data_cyto$cell_t0[match(robot_id, data_cyto$robot_id)]) %>%
  # filter(!if_any(everything(), ~ is.na(.x))) %>%
  group_by(parameter) %>%
  nest() %>%
  
  # Estimation of the block effect as a random effect and as a fixed effect. baker_id is removed of the model because there is only one baker in this data
  mutate(fit_block_ranef = map(data, ~ lmer(value ~ flour_type * wheat_type + backslopping + (1|bloc_month) + (1|cell_t0) + (1|strain_name), data = .)),
         fit_block_fixef = map(data, ~ lm(value ~ bloc_month * cell_t0 * strain_name, data = .))) %>%
  mutate(ranef_bloc_int = map_dbl(fit_block_ranef, ~ summary(.x)$coeff[1,1]),
         ranef_bloc_1 = map_dbl(fit_block_ranef, ~ ranef(.x)$bloc_month[1,1]),
         ranef_bloc_2 = map_dbl(fit_block_ranef, ~ ranef(.x)$bloc_month[2,1]),
         fixef_bloc_int = map_dbl(fit_block_fixef, 
                                  ~ .x$coefficients["bloc_month2"]), 
         fixef_bloc_1 = - fixef_bloc_int/2, 
         fixef_bloc_2 = fixef_bloc_int/2, 
         p = map_dbl(fit_block_fixef, ~ anova(.x)$P[1])) %>%
  select(-data, -fit_block_fixef) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p, method = "BH"))

fits_blockeff_estim %>%
  mutate(summary = map(fit_block_ranef, augment)) %>%
  unnest(summary) %>%
  select(parameter, flour_type, wheat_type, backslopping, .fitted, .resid) %>%
  ggplot() +
  aes(x = .fitted, y = .resid, color = paste(flour_type, wheat_type, backslopping)) +
  geom_point() +
  facet_wrap(~ parameter, scales = "free") +
  theme_minimal() +
  scale_color_discrete(name = "") +
  theme(strip.text.y = element_text(angle = 0),
        title = element_text(face = "bold"),
        legend.position = "bottom")
myggsave(filename = "output/data_phenot_blockeff_estim_res_fitted_plot", width = 7, height = 7)


## Model fitting after block effect removal ####

fits_blockeff_corr <- data_phenot_parms_clean %>%
  filter(parameter != "cell_t0") %>%
  mutate(flour_id = strains$flour_id[match(strain_name, strains$strain_name)],
         baker_id = strains$baker_id[match(strain_name, strains$strain_name)],
         flour_type = flours$mill_type[match(flour_id, flours$flour_id)],
         wheat_type = flours$wheat_type[match(flour_id, flours$flour_id)],
         backslopping = strains$backslopping[match(strain_name, 
                                                   strains$strain_name)],
         cell_t0 = data_cyto$cell_t0[match(robot_id, data_cyto$robot_id)]) %>%
  mutate(fixef_bloc_1 = fits_blockeff_estim$fixef_bloc_1[match(parameter, 
                                                               fits_blockeff_estim$parameter)],
         fixef_bloc_2 = fits_blockeff_estim$fixef_bloc_2[match(parameter, 
                                                               fits_blockeff_estim$parameter)],
         valueC = ifelse(bloc_month == 1, 
                         value - fixef_bloc_1, 
                         value - fixef_bloc_2)) 

ggplot(fits_blockeff_corr) +
  aes(x = bloc_month, y = value) +
  geom_boxplot() +
  facet_wrap(~ parameter, scales = "free")

  # filter(!if_any(everything(), ~ is.na(.x))) %>%
  group_by(parameter) %>%
  nest() %>%
  mutate(fit_sn_ft_wt_bs_bk_Rc0 = map(data, ~ lmer(valueC ~ strain_name 
                                                       + flour_type * wheat_type 
                                                       + backslopping + baker_id 
                                                       + (1|cell_t0), 
                                                       data = .)),
         fit_sn_ft_wt_bs_bk = map(data, ~ lm(valueC ~ strain_name + flour_type 
                                             * wheat_type + backslopping + 
                                               baker_id, data = .)),
         fit_sn_Rc0 = map(data, ~ lmer(valueC ~ strain_name + (1|cell_t0), 
                                       data = .)),
         fit_sn = map(data, ~ lm(valueC ~ strain_name, data = .)),
         fit_Rsn_Rc0 = map(data, ~ lmer(valueC ~ 1 + (1|strain_name) 
                                        + (1|cell_t0), 
                                            data = .)),
         fit_ft_wt_bs_bk_Rc0 =  map(data, ~ lmer(valueC ~ flour_type 
                                                     * wheat_type + backslopping 
                                                     + baker_id + (1|cell_t0), 
                                                     data = .)),
         fit_ft_wt_bk_Rc0 =  map(data, ~ lmer(valueC ~ flour_type 
                                                  * wheat_type + baker_id 
                                                  + (1|cell_t0), data = .)),
         fit_ft_wt_bs_Rc0 =  map(data, ~ lmer(valueC ~ flour_type 
                                                  * wheat_type + backslopping 
                                                  + (1|cell_t0), data = .)),
         fit_bk_Rc0 =  map(data, ~ lmer(valueC ~ baker_id
                                            + (1|cell_t0), data = .)),
         fit_bk_Rc0_Rsn =  map(data, ~ lmer(valueC ~ baker_id
                                                + (1|cell_t0)
                                                + (1|strain_name), data = .)),
         fit_ft_wt_bs_bk_Rc0_Rsn = map(data, ~ lmer(valueC ~ flour_type 
                                                        * wheat_type 
                                                        + backslopping 
                                                        + baker_id 
                                                        + (1|cell_t0) 
                                                        + (1|strain_name), 
                                                        data = .))) %>%
  
  pivot_longer(cols = starts_with("fit_"), names_to = "fit_name", 
               values_to = "fit")

## Info ####

fits_blockeff_corr_infos <- fits_blockeff_corr %>%
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




# Data structure : distribution of strains among blocks ####
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

# Bazar! ####

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