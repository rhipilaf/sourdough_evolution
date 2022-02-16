

# Weight loss trends ####
data_phenot %>%
  ggplot(.) +
  aes(x = time, y = weight_loss, group = robot_id) +
  geom_point(shape = ".") +
  geom_line(alpha = 0.3) +
  theme_minimal()
myggsave(filename = "output/data_phenot_overview_weight", width = 8, height = 6)


# CO2 cumul trends ####
data_phenot %>%
  mutate(excluded = ifelse(robot_id %in% excluded_replicates$robot_id, TRUE, FALSE)) %>%
  ggplot() +
  aes(x = time, y = co2_cumul, group = robot_id, color = excluded) +
  geom_point(shape = ".") +
  geom_line(aes(size = excluded), alpha = 0.1) +
  scale_color_manual(values = c("black","red")) +
  scale_size_manual(values = c(0.05,2)) +
  theme_minimal()
myggsave(filename = "output/data_phenot_overview_co2_cumul", width = 8, height = 6)


# CO2 flow rate trends ####
ggplot(data_phenot) +
  aes(x = time, y = co2_flowrate_3p, group = robot_id) +
  geom_point(shape = ".") +
  geom_line(alpha = 0.1) +
  theme_minimal()
myggsave(filename = "output/data_phenot_overview_co2_flowrate", width = 8, height = 6)


# Estimated parameters (raw data, including outliers) ####

## Global distribution of parameters ####
data_phenot_parms %>%
  mutate(parameter = add_units(parameter)) %>%
  ggplot(.) +
  aes(x = value) +
  geom_density(fill = alpha("red", 0.5)) +
  facet_wrap(~parameter, scales = "free") +
  theme_minimal()
myggsave(filename = "output/global_distrib_parms", width = 6, height = 5)

## Intra-strain variance ####
data_phenot_parms %>%
  mutate(parameter = add_units(parameter)) %>%
  ggplot(.) + 
  aes(y = value, group = strain_name, x = strain_name) +
  geom_point(alpha = 0.3) +
  geom_line(alpha = 0.3) +
  facet_wrap(~parameter, scales = "free") +
  theme_minimal() +
  theme(axis.text.x = element_blank())
myggsave(filename = "output/intra_strain_var", width = 6, height = 5)

## Intra strain variance distribution ####
data_phenot_parms %>%
  mutate(parameter = add_units(parameter)) %>%
  group_by(strain_name, parameter) %>%
  summarise(strain_var = var(value)) %>%
  ggplot(., aes(x = strain_var)) +
  geom_histogram(fill = alpha("green", 0.5)) +
  facet_wrap(~parameter, scales = "free") +
  theme_minimal()
myggsave(filename = "output/intra_strain_var_distrib", width = 6, height = 5)


# Estimated parameters (clean data, without outliers) ####

## Global distribution of parameters ####
data_phenot_parms_clean %>%
  mutate(parameter = add_units(parameter),
         strain_type = strains$strain_type[match(strain_name, strains$strain_name)]) %>%
  ggplot(.) +
  aes(x = value) +
  geom_density(fill = alpha("red", 0.5)) +
  geom_vline(aes(xintercept = value, color = strain_type)) +
  scale_color_manual(values = strain_types_cols) +
  facet_wrap(~parameter, scales = "free") +
  theme_minimal() +
  theme(legend.position = c(0.7,0.1), axis.text.x = element_text(angle = 90))
myggsave(filename = "output/global_distrib_parms_clean", width = 6, height = 5)

## Intra-strain variance ####
data_phenot_parms_clean %>%
  mutate(parameter = add_units(parameter)) %>%
  ggplot(.) + 
  aes(y = value, group = strain_name, x = strain_name) +
  geom_point(alpha = 0.3) +
  geom_line(alpha = 0.3) +
  facet_wrap(~parameter, scales = "free") +
  theme_minimal() +
  theme(axis.text.x = element_blank())
myggsave(filename = "output/intra_strain_var_clean", width = 6, height = 5)

## Intra strain variance distribution ####
data_phenot_parms_clean %>%
  mutate(parameter = add_units(parameter)) %>%
  group_by(strain_name, parameter) %>%
  summarise(strain_var = var(value)) %>%
  ggplot(., aes(x = strain_var)) +
  geom_histogram(fill = alpha("green", 0.5)) +
  facet_wrap(~parameter, scales = "free") +
  theme_minimal()
myggsave(filename = "output/intra_strain_var_distrib_clean", width = 6, height = 5)

# Effects ####
## Inoculation effect ####
data_phenot_parms_clean %>%
  filter(parameter %in% c("lag","tvmax","cell_t0")) %>%
  select(robot_id, parameter, value) %>% 
  unique() %>%
  pivot_wider(names_from = parameter, values_from = value) %>%
  pivot_longer(cols = c("lag","tvmax"), names_to = "parameter", values_to = "value") %>%
  ggplot(aes(y = value, x = cell_t0)) +
  facet_wrap(~ parameter, scales = "free") +
  geom_point() +
  theme_minimal()
myggsave(filename = "output/effects_inoculation", width = 6, height = 5)


## Flour type effect ####
data_phenot_parms_clean %>%
  filter(parameter != "cell_t0") %>%
  mutate(flour_id = strains$flour_id[match(strain_name, strains$strain_name)], 
         flour_type = flours$mill_type[match(flour_id, flours$flour_id)]) %>%
  filter(!is.na(flour_type)) %>%
  ggplot() +
  aes(x = flour_type, y = value) +
  geom_boxplot() +
  facet_wrap(~parameter, scales = "free") +
  theme_minimal()
myggsave(filename = "output/effects_flour", width = 7, height = 6)


## Wheat type effect ####
data_phenot_parms_clean %>%
  filter(parameter != "cell_t0") %>%
  mutate(flour_id = strains$flour_id[match(strain_name, strains$strain_name)], 
         wheat_type = flours$wheat_type[match(flour_id, flours$flour_id)]) %>%
  filter(!is.na(wheat_type)) %>%
  ggplot() +
  aes(x = wheat_type, y = value) +
  geom_boxplot() +
  facet_wrap(~parameter, scales = "free") +
  theme_minimal()
myggsave(filename = "output/effects_wheat", width = 7, height = 6)


## Baker effect ####
data_phenot_parms_clean %>%
  filter(parameter != "cell_t0") %>%
  mutate(baker_id = strains$baker_id[match(strain_name, strains$strain_name)],
         strain_type = strains$strain_type[match(strain_name, strains$strain_name)]) %>%
  filter(!is.na(baker_id)) %>%
  ggplot() +
  aes(x = baker_id, y = value) +
  geom_boxplot() +
  geom_point(aes(color = strain_type)) +
  scale_color_manual(values = strain_types_cols) +
  facet_wrap(~parameter, scales = "free") +
  theme_minimal()
myggsave(filename = "output/effects_baker", width = 7, height = 6)


## Backslopping effect ####
data_phenot_parms_clean %>%
  filter(parameter != "cell_t0") %>%
  mutate(backslopping = strains$backslopping[match(strain_name, strains$strain_name)]) %>%
  filter(!is.na(backslopping)) %>%
  ggplot() +
  aes(x = backslopping, y = value) +
  geom_boxplot() +
  facet_wrap(~parameter, scales = "free") +
  theme_minimal()
myggsave(filename = "output/effects_backslopping", width = 7, height = 6)

## Block effect ####
data_phenot_parms_clean %>%
  filter(parameter != "cell_t0") %>%
  ggplot() +
  aes(x = bloc_month , y = value) +
  geom_boxplot() +
  facet_wrap(~parameter, scales = "free") +
  theme_minimal()
myggsave(filename = "output/effects_block_month", width = 7, height = 6)


strains_in_both <- data_phenot_parms_clean %>%
  select(strain_name, bloc_month) %>%
  unique() %>%
  group_by(strain_name) %>%
  summarise(n = n()) %>%
  filter(n > 1) %>%
  select(strain_name) %>% 
  flatten_chr()

data_phenot_parms_clean %>%
  filter(parameter != "cell_t0",
         strain_name %in% strains_in_both) %>%
  ggplot() +
  aes(x = bloc_month , y = value) +
  geom_boxplot() +
  facet_wrap(~parameter, scales = "free") +
  theme_minimal()
myggsave(filename = "output/effects_block_month_onlyshared", width = 7, height = 6)
