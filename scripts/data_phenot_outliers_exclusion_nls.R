
## Excluding outliers
excluded_replicates <- data_phenot_parms %>%
  filter(robot_id %in% robot_id[parameter == "co2max" & value < 20] |
           robot_id %in% robot_id[parameter == "tvmax" & value > 21.75]
         # robot_id %in% robot_id[parameter == "tvmax" & value > 19]
  ) %>%
  filter(parameter %in% c("co2max","tvmax")) %>%
  pivot_wider(names_from = parameter, values_from = value) %>%
  select(robot_id, co2max, tvmax) %>% 
  mutate_if(is.numeric, ~ round(., 1)) %>%
  mutate(strain_name = data_cyto$strain_name[match(robot_id, 
                                                   data_cyto$robot_id)],
         reason = data_cyto$notes[match(robot_id, data_cyto$robot_id)]) %>%
  unique() %>%
  arrange(strain_name, robot_id)

## Excluding outlier replicates
data_phenot_clean <- data_phenot %>% 
  filter(!(robot_id %in% excluded_replicates$robot_id))

data_phenot_parms_clean <- data_phenot_parms %>% 
  filter(!(robot_id %in% excluded_replicates$robot_id))

if (file.exists(u)) 