library(data.table)

# Identify samples with an F coefficient within 3SD of the cohort mean

# Input and output files
args <- commandArgs(trailingOnly=True)
input_dir <- args[1]

# Load data
dat <- fread(paste(input_dir, "HETEROZYGOSITY.het", sep='/'))

# Identify samples
valid <- dat[F<=mean(F)+3*sd(F) & F>=mean(F)-3*sd(F)] 
invalid <- dat[F>mean(F)+3*sd(F) | F<mean(F)-3*sd(F)]

# print FID and IID for valid samples
fwrite(valid[,c("FID","IID")], paste(input_dir, "HETEROZYGOSITY.valid.sample", sep='/'), sep="\t") 
fwrite(invalid[,c("FID","IID","F")], paste(input_dir, "HETEROZYGOSITY.invalid.sample", sep='/'), sep="\t")
