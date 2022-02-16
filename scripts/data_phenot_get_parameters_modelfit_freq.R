
# Computes of CO2 cumulation and CO2 flow rate at each time t

source("scripts/functions_co2_estimation.R")

data_phenot %<>%
  group_by(robot_id) %>%
  arrange(robot_id, time) %>%
  mutate(co2_cumul = weight_to_cumul(w = weight_loss, start = 2, t.ref = c(2), 
                                     vol = 15.2),
         co2_flowrate_3p = moving_diff(x = co2_cumul, ti = time, l = 3)) %>%
  ungroup()


if(file.exists("data/data_robot/data_phenot_parms_nls.rds")) {
  
  data_phenot_parms <- readRDS("data/data_robot/data_phenot_parms_nls.rds")
  data_phenot_parms_nls_info <- readRDS("data/data_robot/data_phenot_parms_nls_info.rds")
  data_phenot_parms_nls_preds <- readRDS("data/data_robot/data_phenot_parms_nls_preds.rds")
  
} else {
  
  # Gompertz growth model (classic model)
  # Gompertz <- function(x, y0, ymax, k, lag){
  #   result <- y0 + (ymax -y0)*exp(-exp(k*(lag-x)/(ymax-y0) + 1) )
  #   return(result)
  # }
  
  # partial derivative in x (=time) of the Gompertz curve
  #GompD <- mosaicCalc::D(y0 + (ymax -y0)*exp(-exp(k*(lag-x)/(ymax-y0) + 1)) ~ x)
  
  # Modified Gompertz growth model (Mohammadi et al. 2014)
  Gompertz <- function(t, pmax, vmax, lambda){
    
    result <- pmax*exp(-exp((vmax*exp(1))*(lambda-t)/(pmax) + 1))
    
    return(result)
    
  }

  # partial derivative in x (=time) of the Gompertz curve
  GompD <- mosaicCalc::D(pmax*exp(-exp((vmax*exp(1))*(lambda-t)/(pmax) + 1)) ~ t)
  #pmax*(exp(-exp((vmax*exp(1))*(lambda-t)/(pmax)+1))*(exp((vmax*exp(1))*(lambda-t)/(pmax)+1)*((vmax * exp(1))/(pmax))))
  
  # estimate growth parameters
  cat("Fitting Gompertz curves to the data...\n")
  data_phenot_parms_fits <- data_phenot %>%
    group_by(robot_id) %>%
    nest() %>%
    mutate(fit = purrr::map(data, ~ nls.multstart::nls_multstart(co2_cumul ~ Gompertz(time, y0, ymax, k, lag),
                                                                 data = .x,
                                                                 iter = 1000,
                                                                 start_lower = c(y0 = 0, ymax = 20, k = 2, lag = 4),
                                                                 start_upper = c(y0 = 0.5, ymax = 25, k = 4.5, lag = 12),
                                                                 supp_errors = 'Y',
                                                                 na.action = na.omit,
                                                                 lower = c(y0 = 0, ymax = 0, k = 0, lag = 0))))
  
  # get models goodness-of-fit indicators
  data_phenot_parms_nls_info <- data_phenot_parms_fits %>%
    mutate(summary = map(fit, glance)) %>%
    unnest(summary) %>%
    select(-data, -fit)
  
  # get models parameters
  data_phenot_parms <- data_phenot_parms_fits %>%
    mutate(p = map(fit, tidy)) %>%
    unnest(p) %>%
    select(-data, -fit)
  
  # get confidence intervals
  data_phenot_parms_CI <- data_phenot_parms_fits %>%
    mutate(cis = map(fit, nlstools::confint2),
           cis = map(cis, data.frame)) %>%
    unnest(cis) %>%
    rename(conf.low = X2.5.., conf.high = X97.5..) %>%
    group_by(robot_id) %>%
    mutate(term = c('y0', 'ymax', 'k', 'lag')) %>%
    ungroup() %>%
    select(-data, -fit)
  
  # get predictions
  data_phenot_parms_nls_preds <- data_phenot_parms_fits %>%
    mutate(., p = map(fit, augment)) %>%
    unnest(p) %>%
    select(-data, -fit)
  
  # merge parameters and CI estimates
  data_phenot_parms <- merge(data_phenot_parms, data_phenot_parms_CI, 
                             by = intersect(names(data_phenot_parms), 
                                            names(data_phenot_parms_CI)))
  
  # estimation of tvmax and data formatting
  data_phenot_parms %<>%
    select(robot_id, term, estimate) %>%
    group_by(robot_id) %>%
    summarise(tvmax = optimize(GompD, 
                               y0 = estimate[term == "y0"], 
                               ymax = estimate[term == "ymax"], 
                               lag = estimate[term == "lag"],  
                               k = estimate[term == "k"],  
                               maximum = T, interval = c(0, 40))$maximum) %>%
    pivot_longer(cols = "tvmax", names_to = "term", values_to = "estimate") %>%
    bind_rows(., data_phenot_parms) %>%
    mutate(term = case_when(term == "k" ~ "co2max",
                            term == "lag" ~ "lag",
                            term == "tvmax" ~ "tvmax",
                            term == "ymax" ~ "vmax",
                            TRUE ~ term)) %>%
    bind_rows(., data_cyto %>%
                select(robot_id, cell_t0, cell_t27, death_prct) %>%
                unique() %>%
                pivot_longer(cols = c("cell_t0", "cell_t27","death_prct"), 
                             names_to = "term", values_to = "estimate")) %>%
    arrange(robot_id, term) %>%
    mutate(strain_name = data_cyto$strain_name[match(robot_id, data_cyto$robot_id)],
           bloc = data_cyto$bloc[match(robot_id, data_cyto$robot_id)]) %>%
    rename(parameter = term, value = estimate) %>%
    relocate(robot_id, strain_name, bloc)
    


  # ggplot() +
  #   geom_function(fun = Gompertz, 
  #                 args = list(y0 = 0.736, ymax = 23.3, lag = 10.6, k = 8.99)) +
  #   geom_function(fun = GompD, 
  #                 args = list(y0 = 0.736, ymax = 23.3, lag = 10.6, k = 8.99)) +
  #   xlim(0, 30)
  
  saveRDS(data_phenot_parms, "data/data_robot/data_phenot_parms_nls.rds")
  saveRDS(data_phenot_parms_nls_info, "data/data_robot/data_phenot_parms_nls_info.rds")
  saveRDS(data_phenot_parms_nls_preds, "data/data_robot/data_phenot_parms_nls_preds.rds")
  
  
}