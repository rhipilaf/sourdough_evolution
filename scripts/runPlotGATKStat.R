#!/usr/bin/env Rscript

# Retrieve args
args <- commandArgs(TRUE)

dth <-list(
  QD = 25,
  FS = 10,
  SOR = 2,
  MQ = 50,
  MQRS = c(-2.5, 2.5),
  RPRS = c(-3,3)
)

# Should contain 2 values
if( length(args) < 2) {
  stop("Missing arguments.")
}
if( length(args) > 3) {
  stop("Too many arguments.")
}
# 1) Input path
# 2) Output path (pdf)
# 3) Edit threshold
input <- NULL
if( file.exists(args[1]) ) {
  input <- args[1]
} else {
  if( dir.exists(args[1]) ) {
    input <- dir(path = args[1], full.names = TRUE)
  } else {
    dn <- dirname(args[1])
    fn <- basename(args[1])
    input <- dir(path = dn, pattern = fn, full.names = TRUE)
  }
}
if( length(input) == 0 ) {
  stop("Failed to find input stat file(s).")
}
output <- args[2]

if( length(args) == 3 ) {
  # Parse new threshold values
  tmp <- strsplit(args[3], ";")[[1]]
  for(kv in tmp) {
    th <- strsplit(kv, "=")[[1]]
    dth[[th[1]]] <- as.numeric(strsplit(th[2], ",")[[1]])
  }
}
  

# Required plot function
source("./plotGATKStat.R")

pdf(output, width=8, height=9)
for(i in input) {
  exp_id <- basename(i)
  m <- read.table(i, stringsAsFactors = F, header = T, fill = T)
  plotGATKStat(m, exp_id=exp_id, QD=dth[["QD"]], FS=dth[["FS"]], SOR=dth[["SOR"]],
               MQ=dth[["MQ"]], MQRS=dth[["MQRS"]], RPRS=dth[["RPRS"]])
}

