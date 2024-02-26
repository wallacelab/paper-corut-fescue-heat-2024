# Purpose: Plot dispersion estimates from DESeq2
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library("tidyverse")
library("DESeq2")

# load deseq2 data
dds <- readRDS(snakemake@input[[1]])

png(snakemake@output[[1]], width=12, height=10, units="in", res=350)
plotDispEsts(dds)
dev.off()