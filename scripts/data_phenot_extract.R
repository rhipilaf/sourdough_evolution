
path = "data/data_robot/raw_data_lauriane/"
files = list.files(path)
output_file_name = "data/data_robot/data_phenot.csv"


# Importing all files in a list
data_robot_phenot_list <- lapply(files, function(x) read.table(paste0(path,x), 
                                                               header=FALSE, 
                                                               sep = "\t", 
                                                               dec = "."))
names(data_robot_phenot_list) <- sub(files, pattern = "_tab.txt", replacement = "")

# Adding a new column giving the name of the file it's coming from to each data.frame
data_robot_phenot_list <- lapply(seq_along(data_robot_phenot_list), 
                                 function(i) cbind(robot_id = rep(names(data_robot_phenot_list)[i], 
                                                              nrow(data_robot_phenot_list[[i]])),
                                                   data_robot_phenot_list[[i]]))

# Stacking them all together in a single data.frame
data_robot_phenot_df <- do.call(rbind.data.frame, data_robot_phenot_list)[,1:4]
names(data_robot_phenot_df) <- c("robot_id","date_hour","time","weight_loss")

# Extracting the position and tube format information from the file name.
data_robot_phenot_df %<>%
  mutate(position = str_extract(robot_id, pattern = "[0-9]{3}$"),
         tube_format = str_extract(robot_id, pattern = "^.+(?=-[0-9]{8})"))

# Writing the raw data as a csv file
write.csv(x = data_robot_phenot_df, file = output_file_name, row.names = F)

