
# Overview of weight loss trends
data_phenot_overview_weight <- ggplot(data_phenot) +
  aes(x = time, y = weight_loss, group = robot_id) +
  geom_point(shape = ".") +
  geom_line(alpha = 0.3) +
  theme_minimal()

ggsave(plot = data_phenot_overview_weight, filename = "output/data_phenot_overview_weight.pdf", width = 8, height = 6)


# Overview of CO2 cumul trends
data_phenot_overview_co2_cumul <- ggplot(data_phenot) +
  aes(x = time, y = co2_cumul, group = robot_id, color = date_hour) +
  geom_point(shape = ".") +
  geom_line(alpha = 0.1) +
  scale_color_continuous(type = "viridis") +
  theme_minimal()

ggsave(plot = data_phenot_overview_co2_cumul, filename = "output/data_phenot_overview_co2_cumul.pdf", width = 8, height = 6)


# Overview of CO2 cumul trends
data_phenot_overview_co2_flowrate <- ggplot(data_phenot) +
  aes(x = time, y = co2_flowrate, group = robot_id, color = date_hour) +
  geom_point(shape = ".") +
  geom_line(alpha = 0.1) +
  scale_color_continuous(type = "viridis") +
  theme_minimal()

ggsave(plot = data_phenot_overview_co2_flowrate, filename = "output/data_phenot_overview_co2_flowrate.pdf", width = 8, height = 6)
