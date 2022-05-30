

# METADATA ####

data_phenot_parms_clean_LBMB <- readxl::read_xlsx("data/eqtable_strains_labels.xlsx", col_types = "text")
eqt_chrom <- readxl::read_xlsx("data/eqtable_chrom.xlsx")

strains <- readxl::read_xlsx("data/strains.xlsx", col_types = 'text') %>%
  mutate(backslopping = factor(backslopping, levels = c("T1", "T5", "T32"))) %>%
  filter(baker_id %in% c("LB", "AM", "MB", "PDC"))
bakers <- readxl::read_xlsx("data/bakers.xlsx", col_types = 'text')
flours <- readxl::read_xlsx("data/flours.xlsx", col_types = 'text')
backsloppings <- readxl::read_xlsx("data/backsloppings.xlsx", 
                                   col_types = 'text') %>%
  mutate(bs_time = factor(bs_time, levels = c("Ti", "Tf")))

chrlen <- read.table("data/data_geno/01_refgenome/S288c_ABC.length", 
                     col.names = c("chr_long","length")) %>%
  mutate(cumlength = cumsum(length),
         chr = convnames(chr_long, eqt_chrom, from = "vcf", to = "std_short"))

strains_scer_all <- readxl::read_xlsx("data/strains_scer_all.xlsx", 
                                      sheet = 1) %>%
  mutate_at(vars(Scerevisiae, dup, subset, bread, included, evoexp), as.logical) %>%
  mutate(group = case_when(evoexp ~ "boul_evoexp",
                           bread & type %in% c("Dough", "Sourdough") ~ "boul_sdgh",
                           bread & type == "Commercial" ~ "boul_comm",
                           !bread ~ "other"),
         group = factor(group, levels = c("boul_evoexp", "boul_sdgh",
                                          "boul_comm", "other"),
                        labels = c("Boulangerie/Evolution expérimentale", 
                                   "Boulangerie/Levain", 
                                   "Boulangerie/Commerciale",
                                   "Autre origine")))

meta_evoexp <- list()
meta_evoexp$data_raw <- readxl::read_xlsx("data/data_lauriane/data_evoexp.xlsx",
                                          sheet = "data_raw") %>%
  separate(sdgh_id, into = c("flour_id","baker_id"), sep = "-") %>%
  mutate(backslopping = toupper(backslopping),
         date = as.character(date),
         hour = gsub("^(.+) ","", as.character(hour)),
         time = ymd_hms(paste(date, hour)))

meta_evoexp$data_rich <- readxl::read_xlsx("data/data_lauriane/data_evoexp.xlsx",
                                           sheet = "data_rich") %>%
  mutate(backslopping = paste0("T", backslopping),
         T_sourdough = as.numeric(T_sourdough),
         T_fournil = as.numeric(T_sourdough))

meta_evoexp <- full_join(meta_evoexp$data_raw, meta_evoexp$data_rich,
                 by = c("flour_id", "baker_id", "backslopping")) %>%
  relocate(flour_id, baker_id, backslopping, time)


# PHENOTYPIC & FITNESS DATA ####

data_cyto <- readxl::read_xlsx("data/data_robot/data_cyto.xlsx", 
                               col_types = 'text') %>%
  filter(project == "Lauriane") %>% # only strains from Lauriane
  select(c("date","strain_name", "robot", "robot_id","robot_inoc_id",
           "position", "bloc","cell_t0", "cell_t27","death_prct","notes")) %>%
  rename(pdead = death_prct) %>%
  mutate(bloc = as.factor(as.character(round(as.numeric(bloc)))),
         cell_t0 = as.numeric(cell_t0),
         cell_t27 = as.numeric(cell_t27),
         pdead = as.numeric(pdead),
         npop = cell_t27 * (1 - pdead/100),
         ntot = cell_t27,
         date = as.Date(date, format = "%Y_%m_%d"),
         bloc_month = as.factor(case_when(bloc %in% 3:5 ~ 1,
                                          bloc %in% 7:9 ~ 2))) %>%
  select(-cell_t27)

data_phenot <- read.csv("data/data_robot/data_phenot.csv") %>%
  filter(robot_id %in% data_cyto$robot_id) %>%
  mutate(strain_name = data_cyto$strain_name[match(robot_id, 
                                                   data_cyto$robot_id)],
         date_hour = as.Date(date_hour, format = "%d/%m/%Y %H:%M:%S")) %>%
  # relocate(robot_id, strain_name) %>% 
  as_tibble()

## Phenotypic parameters import
# source("scripts/data_phenot_get_parameters_std.R")
# source("scripts/data_phenot_get_parameters_modelfit_freq.R")
# source("scripts/data_phenot_get_parameters_modelfit_bayes.R")


# DOUGH (we do the hypothesis that the dough density is at a density of 0).

C_t0_g <- 1e7 #Number of cells initially inoculated per gram of dough
C_t0_mg <- C_t0_g / 1000 # /1000 = To convert to mg(µL)

data_dough_meta <- xlsx::read.xlsx("data/data_dough/data_dough.xlsx", 
                                   sheetName = "metadata") %>%
  pivot_longer(cols = starts_with("time"), 
               names_to = "time", values_to = "hour") %>%
  mutate(time = as.numeric(gsub("[a-z]|(_)","",time))) %>%
  group_by(strain_name, replicate) %>%
  mutate(time_real = as.numeric(hour - hour[time == 0], "hours"))

data_dough_man <- xlsx::read.xlsx("data/data_dough/data_dough.xlsx", 
                                  sheetName = "count_manual") %>%
  mutate_at(vars(names(.)[grep("(^X)|(^hmax)", names(.))]), ~ as.character(.)) %>%
  pivot_longer(cols = names(.)[grep("^X", names(.))],
               values_to = "cfu") %>%
  separate(name, into = c("time","dilution","petriplate"), sep = "_") %>% # 
  mutate(dilution = 10^as.numeric(dilution), # converts as power of ten
         time = as.numeric(gsub("X|h","",time)),
         comment = ifelse(grepl("[a-z]", cfu), cfu, NA),
         cfu = as.numeric(ifelse(grepl("[a-zA-Z]", cfu), 0, cfu)),
         method = "manual",) %>% # removes characters
  complete(strain_name, replicate, time) %>%
  
  # For the real strains, the theoritical cells concentration (10^-7/g) is set for t0
  mutate(cfu = ifelse(grepl("^LB", strain_name) & time == 0 & dilution == 1, 
                      C_t0_mg, cfu),
         cfu = ifelse(grepl("TEMOIN", strain_name) & replicate == "a" & 
                        time == 0 & dilution == 1, 
                      0, cfu)) %>%
  filter(!is.na(cfu)) %>%
  select(strain_name, replicate, time, petriplate, cfu, dilution, method)

data_dough_auto <- xlsx::read.xlsx("data/data_dough/data_dough.xlsx", 
                                   sheetName = "count_auto") %>%
  select(strain_name, sample, cfu, dilution) %>%
  .[-nrow(.),] %>%
  separate(sample, c("flour_id","backslopping","replicate",
                     "time","petriplate")) %>%
  mutate(baker_id = "LB",
         time = as.numeric(gsub("H|h", "", time)),
         dilution = 1/dilution,
         replicate = tolower(replicate),
         method = "auto") %>%
  select(strain_name, replicate, time, petriplate, cfu, dilution, method)

data_dough <- rbind(data_dough_auto, data_dough_man) %>%
  mutate(dilution = 1/dilution) %>%
  separate(col = strain_name, c("baker_id","flour_id","backslopping","tmp"),
           remove = F, sep = "-") %>%
  mutate(baker_id = ifelse(grepl("(TEMOIN)|(CONTROLE)", baker_id), 
                           NA, baker_id)) %>%
  select(-tmp) %>%
  left_join(data_dough_meta, by = c("strain_name","replicate","time")) %>%
  mutate(C_mg = ifelse(dilution < 1, cfu/(dilution*50), cfu),
         C_g = C_mg * 1000,
         # Weighing observations according to their level of reliability depending on the level of dilution
         w = log10(dilution) + max(abs(log10(dilution))) + 1) %>% 
  arrange(baker_id, flour_id, backslopping, time, dilution)


rm(data_dough_auto, data_dough_man, data_dough_meta)





# GENOMIC DATA ####


