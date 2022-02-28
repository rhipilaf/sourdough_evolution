
# Computes of CO2 cumulation and CO2 flow rate at each time t

source("scripts/functions_co2_estimation.R")

data_phenot %<>%
  group_by(robot_id) %>%
  arrange(robot_id, time) %>%
  mutate(co2_cumul = weight_to_cumul(w = weight_loss, start = 2, t.ref = c(2), 
                                     vol = 15.2),
         co2_flowrate_3p = moving_diff(x = co2_cumul, ti = time, l = 3)) %>%
  ungroup()


if(file.exists("data/data_robot/data_phenot_parms_std.rds")) {
  
  data_phenot_parms <- readRDS("data/data_robot/data_phenot_parms_std.rds")
  
  
} else {
  
  # Extracts statistiques from the CO2 cumulation (g) and CO2 flow rate (g/h)
  data_phenot_parms <- data_phenot %>%
    group_by(robot_id) %>%
    summarise(co2max = max(co2_cumul), # Maximum cumulated CO2
              vmax = max(co2_flowrate_3p), # Maximum CO2 flow rate
              tvmax = time[which(co2_flowrate_3p == max(co2_flowrate_3p))][1], # Time at the maximum CO2 flow rate
              lag = time[min(which(co2_cumul > 1))]) %>% # Latency : time at which cumulated CO2 reaches 1 for the first time
    ungroup() %>%
    mutate(strain_name = data_cyto$strain_name[match(robot_id, data_cyto$robot_id)],
           cell_t0 = data_cyto$cell_t0[match(robot_id, data_cyto$robot_id)],
           pop_size = data_cyto$pop_size[match(robot_id, data_cyto$robot_id)],
           death_prct = data_cyto$death_prct[match(robot_id, data_cyto$robot_id)],
           bloc = data_cyto$bloc[match(robot_id, data_cyto$robot_id)],
           bloc_month = data_cyto$bloc_month[match(robot_id, data_cyto$robot_id)]) %>%
    pivot_longer(cols = c("co2max", "vmax", "tvmax", "lag", 
                          "cell_t0", "pop_size", "death_prct"),
                 names_to = "parameter", values_to = "value") %>%
    arrange(robot_id, parameter)
  
  
  saveRDS(data_phenot_parms, file = "data/data_robot/data_phenot_parms_std.rds")
  
}
