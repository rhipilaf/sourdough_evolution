
# METADATA

strains <- readxl::read_xlsx("data/strains.xlsx")
bakers <- readxl::read_xlsx("data/bakers.xlsx")
flours <- readxl::read_xlsx("data/flours.xlsx")


# PHENOTYPIC & FITNESS DATA

data_cyto <- readxl::read_xlsx("data/data_robot/data_cyto.xlsx") %>%
  select(c("date","strain_name", "robot", "robot_id","robot_inoc_id",
           "position", "bloc","cell_t0", "cell_t27","death_prct","notes")) %>%
  mutate(bloc = as.factor(as.character(bloc)),
         cell_t0 = as.numeric(cell_t0),
         cell_t27 = as.numeric(cell_t27),
         death_prct = as.numeric(death_prct),
         date = as.Date(date, format = "%Y_%m_%d"))

data_phenot <- read.csv("data/data_robot/data_phenot.csv") %>%
  mutate(date_hour = as.Date(date_hour, format = "%d/%m/%Y %H:%M:%S"))

# source("scripts/data_phenot_get_parameters_std.R")
source("scripts/data_phenot_get_parameters_modelfit_freq.R")
# source("scripts/data_phenot_get_parameters_modelfit_bayes.R")

## Excluding outliers
excluded_replicates <- data_phenot_parms %>% 
  filter(co2max < 20) %>% 
  select(robot_id) %>% 
  mutate(reason = data_cyto$notes[match(robot_id, data_cyto$robot_id)])

## Excluding outlier replicates
data_phenot_clean <- data_phenot %>% filter(!(robot_id %in% excluded_replicates$robot_id))
data_phenot_parms_clean <- data_phenot_parms %>% filter(!(robot_id %in% excluded_replicates$robot_id))

