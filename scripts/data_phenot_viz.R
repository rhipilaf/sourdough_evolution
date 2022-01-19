

# Weight loss trends ####
ggplot(data_phenot) +
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

## Formatting data ####
data_phenot_parms_tidy <- data_phenot_parms %>%
  mutate(bloc = data_cyto$bloc[match(robot_id, data_cyto$robot_id)]) %>%
  pivot_longer(cols = c("co2max", "vmax", "tvmax", "t1g"),
               names_to = "parameter", values_to = "value") %>%
  mutate(parameter = case_when(parameter == "co2max" ~ paste(parameter, "(g)"),
                               parameter == "vmax" ~ paste(parameter, "(g/h)"),
                               parameter == "tvmax" ~ paste(parameter, "(h)"),
                               parameter == "t1g" ~ paste(parameter, "(h)")))

## Global distribution of parameters ####
data_phenot_parms_tidy %>%
  ggplot(.) +
  aes(x = value) +
  geom_density(fill = alpha("red", 0.5)) +
  facet_wrap(~parameter, scales = "free") +
  theme_minimal()
myggsave(filename = "output/global_distrib_parms", width = 6, height = 5)

## Intra-strain variance ####
data_phenot_parms_tidy %>%
  ggplot(.) + 
  aes(y = value, group = strain_name, x = strain_name) +
  geom_point(alpha = 0.3) +
  geom_line(alpha = 0.3) +
  facet_wrap(~parameter, scales = "free") +
  theme_minimal() +
  theme(axis.text.x = element_blank())
myggsave(filename = "output/intra_strain_var", width = 6, height = 5)

## Intra strain variance distribution ####
data_phenot_parms_tidy %>%
  group_by(strain_name, parameter) %>%
  summarise(strain_var = var(value)) %>%
  ggplot(., aes(x = strain_var)) +
  geom_histogram(fill = alpha("green", 0.5)) +
  facet_wrap(~parameter, scales = "free") +
  theme_minimal()
myggsave(filename = "output/intra_strain_var_distrib", width = 6, height = 5)


# Estimated parameters (clean data, without outliers) ####

## Formatting data ####
data_phenot_parms_clean_tidy <- data_phenot_parms_clean %>%
  mutate(bloc = data_cyto$bloc[match(robot_id, data_cyto$robot_id)]) %>%
  pivot_longer(cols = c("co2max", "vmax", "tvmax", "t1g"),
               names_to = "parameter", values_to = "value") %>%
  mutate(parameter = case_when(parameter == "co2max" ~ paste(parameter, "(g)"),
                               parameter == "vmax" ~ paste(parameter, "(g/h)"),
                               parameter == "tvmax" ~ paste(parameter, "(h)"),
                               parameter == "t1g" ~ paste(parameter, "(h)")))

## Global distribution of parameters ####
data_phenot_parms_clean_tidy %>%
  ggplot(.) +
  aes(x = value) +
  geom_density(fill = alpha("red", 0.5)) +
  facet_wrap(~parameter, scales = "free") +
  theme_minimal()
myggsave(filename = "output/global_distrib_parms_clean", width = 6, height = 5)

## Intra-strain variance ####
data_phenot_parms_clean_tidy %>%
  ggplot(.) + 
  aes(y = value, group = strain_name, x = strain_name) +
  geom_point(alpha = 0.3) +
  geom_line(alpha = 0.3) +
  facet_wrap(~parameter, scales = "free") +
  theme_minimal() +
  theme(axis.text.x = element_blank())
myggsave(filename = "output/intra_strain_var_clean", width = 6, height = 5)

## Intra strain variance distribution ####
data_phenot_parms_clean_tidy %>%
  group_by(strain_name, parameter) %>%
  summarise(strain_var = var(value)) %>%
  ggplot(., aes(x = strain_var)) +
  geom_histogram(fill = alpha("green", 0.5)) +
  facet_wrap(~parameter, scales = "free") +
  theme_minimal()
myggsave(filename = "output/intra_strain_var_distrib_clean", width = 6, height = 5)

# Block effect ####
data_phenot_parms_clean_tidy %>%
  ggplot(.) +
  aes(y = value, fill = bloc) +
  geom_boxplot(alpha = 0.3) +
  facet_wrap(~parameter, scales = "free") +
  theme_minimal()

bloc_effect_model <- data_phenot_parms_clean_tidy %>% 
  group_by(parameter) %>%
  summarise(p.value = anova(lm(value ~ bloc))$'Pr(>F)'[1])


# Inoculation effect ####
plot(t1g ~ cell_t0, data_phenot_parms_clean)
