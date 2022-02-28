
source("scripts/data_phenot_get_parameters_nls.R")
dpp_nls <- data_phenot_parms

source("scripts/data_phenot_get_parameters_std.R")
dpp_std <- data_phenot_parms


dpp_comp <- bind_rows(bind_cols(method = "nls", dpp_nls), 
                      bind_cols(method = "std", dpp_std)) %>%
  filter(p.value < 0.01 | is.na(p.value),
         parameter %in% c("co2max", "lag", "tvmax", "vmax"))

ggplot(dpp_comp) +
  aes(x = robot_id, y = value, group = robot_id, color = method) +
  geom_line(color = "black", alpha = 0.5) +
  geom_point(alpha = 0.5, size = 0.5) +
  facet_wrap(~parameter ,scales = "free") +
  theme(axis.text.x = element_blank())
myggsave(filename = "output/compare_param_estim_methods",
        width = 8, height = 6)
