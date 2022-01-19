

if(file.exists("data/data_robot/data_phenot_parms_nls.rds")) {
  
  data_phenot_parms <- readRDS("data/data_robot/data_phenot_parms_nls.rds")
  data_phenot_parms_nls_info <- readRDS("data/data_robot/data_phenot_parms_nls_info.rds")
  data_phenot_parms_nls_preds <- readRDS("data/data_robot/data_phenot_preds.rds")
  
} else {
  
  source("scripts/functions_co2_estimation.R")
  
  # Computes of CO2 cumulation and CO2 flow rate at each time t
  data_phenot %<>%
    group_by(robot_id) %>%
    arrange(robot_id, time) %>%
    mutate(co2_cumul = weight_to_cumul(w = weight_loss, start = 2, t.ref = c(2), 
                                       vol = 15.2),
           co2_flowrate_3p = moving_diff(x = co2_cumul, ti = time, l = 3)) %>%
    ungroup()
  
  # Extracts parameters
  
  Gompertz <- function(x, y0, ymax, k, lag){
    result <- y0 + (ymax -y0)*exp(-exp(k*(lag-x)/(ymax-y0) + 1) )
    return(result)
  }
  
  cat("Fitting Gompertz curves to the data...\n")
  data_phenot_parms_fits <- data_phenot %>%
    group_by(robot_id) %>%
    nest() %>%
    mutate(fit = purrr::map(data, ~ nls.multstart::nls_multstart(co2_cumul ~ Gompertz(time, y0, ymax, k, lag),
                                                  data = .x,
                                                  iter = 1000,
                                                  start_lower = c(y0 = 0, ymax = 20, k = 2, lag = 4),
                                                  start_upper = c(y0 = 0.5, ymax = 25, k = 4.5, lag = 12),
                                                  # supp_errors = 'Y',
                                                  na.action = na.omit,
                                                  lower = c(y0 = 0, ymax = 0, k = 0, lag = 0))))
  
  data_phenot_parms_nls_info <- data_phenot_parms_fits %>%
    mutate(summary = map(fit, glance)) %>%
    unnest(summary) %>%
    select(-data, -fit)
  
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
  
  # merge parameters and CI estimates
  data_phenot_parms <- merge(data_phenot_parms, data_phenot_parms_CI, 
                             by = intersect(names(data_phenot_parms), 
                                            names(data_phenot_parms_CI)))
  
  # get predictions
  data_phenot_parms_nls_preds <- data_phenot_parms_fits %>%
    mutate(., p = map(fit, augment)) %>%
    unnest(p) %>%
    select(-data, -fit)
  
  
  # mutate(y0 = map_dbl(fit, ~ .x$m$getPars()['y0']),
  #        ymax = map_dbl(fit, ~ .x$m$getPars()['ymax']),
  #        k = map_dbl(fit, ~ .x$m$getPars()['k']),
  #        lag = map_dbl(fit, ~ .x$m$getPars()['lag']),
  #        deviance = map_dbl(fit, ~ .x$m$deviance()),
  #        )
  
  # ungroup() %>%
  #   mutate(strain_name = data_cyto$strain_name[match(robot_id, data_cyto$robot_id)],
  #          cell_t0 = data_cyto$cell_t0[match(robot_id, data_cyto$robot_id)],
  #          cell_t27 = data_cyto$cell_t27[match(robot_id, data_cyto$robot_id)],
  #          death_prct = data_cyto$death_prct[match(robot_id, data_cyto$robot_id)])
  
  saveRDS(data_phenot_parms, "data/data_robot/data_phenot_parms_nls.rds")
  saveRDS(data_phenot_parms_nls_info, "data/data_robot/data_phenot_parms_nls_info.rds")
  saveRDS(data_phenot_parms_nls_preds, "data/data_robot/data_phenot_preds.rds")
  
  
}