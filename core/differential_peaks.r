library(DiffBind)

args = commandArgs(trailingOnly=TRUE)
ssheet = read.table(args[1], header=TRUE, stringsAsFactors=TRUE)
results_dir = args[2]

# 