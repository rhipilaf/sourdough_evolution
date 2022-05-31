#!/usr/bin/env Rscript

# Retrieve args
args <- commandArgs(TRUE)

# Expected args
# 1) Input GATK file
# 2) Input sequence length
# 3) Ouput file path (write.table)
if( length(args) != 3 ) {
  stop("The number of input arguments must be 3.")
}
input <- args[1]
in.len <- args[2]
output <- args[3]
if( !file.exists(input) ) {
  stop("The provided input file does not exists (arg #1).")
}
if( !file.exists(in.len) ) {
  stop("The provided input length file does not exists (arg #2).")
}

# Load data
m <- read.table(input, header = TRUE, stringsAsFactors = FALSE)
len <- read.table(in.len, header = FALSE, stringsAsFactors = FALSE)

# Sub-functions
bial_alt_freq <- function(x) {
  x <- as.numeric(x)
  sx <- sum(x)
  if(sx == 0) {
    return(-1)
  } else {
    return(x[2]/sx)
  }
}

# Get sample Col index
col.spl <- grep("\\.AD$", names(m))

# Get allelic freq
m.alt.freq <- apply(m[col.spl], 2, function(co) {
  return(sapply(strsplit(co, ","), bial_alt_freq))
})

# Compute the cumulative location
len.fac <- c(0, cumsum(len[-nrow(len),2]))
names(len.fac) <- len[,1]
len.delta <- len.fac[m[,1]]

# Prepare names
cn <- sub(".AD$", "", colnames(m.alt.freq))
cn <- sub("^X(\\d)", "\\1", cn)
cn <- paste(rep("af.", length(cn)), cn, sep="")

m.tmp <- rbind(
  m.alt.freq,
  matrix(NA, nrow = nrow(len), ncol = ncol(m.alt.freq))
)

colnames(m.tmp) <- cn

mf <- data.frame(
  seq.id = c(m[,1], len[,1]),
  seq.loc = c(m[,2], len[,2]),
  cum.loc = c(m[,2]+len.delta, cumsum(len[,2])),
  all.ref = c(m[,3], rep(NA, nrow(len))),
  all.alt = c(m[,4], rep(NA, nrow(len)))
)
mf <- cbind.data.frame(mf, m.tmp)

write.table(mf, file=output, row.names = FALSE, quote = FALSE)
