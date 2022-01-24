

# Weight loss trends ####
data_phenot %>%
  ggplot(.) +
  aes(x = time, y = weight_loss, group = robot_id) +
  geom_point(shape = ".") +
  geom_line(alpha = 0.3) +
  theme_minimal()
myggsave(filename = "output/data_phenot_overview_weight", width = 8, height = 6)


# CO2 cumul trends ####
ggplot(data_phenot) +
  aes(x = time, y = co2_cumul, group = robot_id) +
  geom_point(shape = ".") +
  geom_line(alpha = 0.1) +
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

# Inoculation effect ####
data_phenot_parms_clean %>%
  filter(parameter %in% c("lag","cell_t0")) %>%
  select(robot_id, parameter, value) %>% 
  unique() %>%
  pivot_wider(names_from = parameter, values_from = value) %>%
  ggplot(aes(y = lag, x = cell_t0)) +
  geom_point() +
  theme_minimal()

# Flour effect
data_phenot_parms_clean %>%
  filter(parameter != "cell_t0") %>%
  mutate(flour_type = )
  ggplot()

# Baker effect

# Backslopping effect