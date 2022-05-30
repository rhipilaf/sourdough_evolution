

intersp <- read.csv("data/data_lauriane/mydata_yeast_all.csv") %>%
  mutate(flour_id = ifelse(grepl("Temoin",flour_id), "Temoin", flour_id),
         bk_fl = paste0(baker_id)) %>%
  group_by(baker_id, flour_id, species) %>%
  summarise(n = n())

ggplot(intersp) +
  geom_bar(aes(x = n, y = bakers, fill = species), stat = "identity") +
  facet_wrap(~baker_id) +
  theme(legend.position = "none")
  