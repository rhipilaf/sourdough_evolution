#!/usr/bin/Rscript

args <- commandArgs()

input = args[6]

path_output = args[7]
output = sub(".bam","_compressed.bam", sub(pattern = "^.+/", replacement = "", x = input, perl = T))

if (file.exists(paste0(path_output,"/",output))) stop(paste0(path_output,"/",output," already exists."))

depth <- read.table(input, col.names = c("chr","location","depth"))

step = 500
window = 1000

locations <- seq(step, nrow(depth)-step, step)

depth_comp <- data.frame(location = 0, depth = 0)

cat("Processing", input)

for(i in 1:length(locations)) {
  
  tmp <- mean(depth$depth[(locations[i]-500):(locations[i]+500)])
  depth_comp[i,] <- c(locations[i], tmp)
  
}

write.table(x = depth_comp, file = paste0(path_output,"/",output), row.names = F)

cat(" -DONE\n")

