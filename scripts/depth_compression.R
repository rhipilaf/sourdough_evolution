##!/usr/bin/env -S Rscript --vanilla

args <- commandArgs()

input = args[1]
output = args[2]
depth <- read.table(input, col.names = c("chr","location","depth"))

step = 500
window = 1000

locations <- seq(step, nrow(depth)-step, step)

depth_comp <- data.frame(location = 0, depth = 0)


cat("Processing", file)

for(i in 1:length(locations)) {
  
  tmp <- mean(depth$depth[(locations[i]-500):(locations[i]+500)])
  depth_comp[i,] <- c(locations[i], tmp)
  
}

write.table(x = depth_comp, file = output, row.names = F)

cat(" -DONE\n")

