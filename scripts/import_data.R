
# Data

strains <- read.csv("data/strains.csv")
bakers <- read.csv("data/bakers.csv")
flours <- read.csv("data/flours.csv")


# Fitness data

data_cyto <- readxl::read_xlsx("data/data_robot/data_cyto.xlsx") %>%
  select(c("date","strain_name", "robot", "robot_id","robot_inoc_id",
           "position", "bloc","cell_t0", "cell_t27","death_prct","notes")) %>%
  mutate(bloc = as.factor(as.character(bloc)))


# Phenotypic data

# data_phenot <- read.csv("data/data_robot/data_phenot.csv") %>%
#   mutate(strain_name = data_cyto$strain_name[match(robot_id, data_cyto$robot_id)])
source("scripts/data_phenot_get_parameters.R")

cyto_what_is_this <- c("robot","modalities","robot_inoc_id","classement_id","position","MTA","MTF","robot_before")



