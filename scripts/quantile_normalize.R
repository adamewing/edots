#!/usr/bin/Rscript

args   <- commandArgs(trailingOnly = TRUE)
if (length(args) == 1) {
    incsv  <- args[1]
    outcsv <- gsub('.csv', '.quantilenorm.csv', incsv)

    library(limma)

    d <- read.table(incsv,header=TRUE, sep=',', row.names=1)
    d.norm <- as.data.frame(normalizeBetweenArrays(as.matrix(d), method='quantile'))
    write.csv(d.norm,file=outcsv)
}
