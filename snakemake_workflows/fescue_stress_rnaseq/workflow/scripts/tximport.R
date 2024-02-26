# Purpose: R script to import Salmon quant files using tximport for use with DESeq2
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# Load libraries
library(tidyverse)
library(tximport)

# Load metadata
sample_metadata <- read.csv(snakemake@input[["metadata"]], header=TRUE, row.names=1, sep="\t")
print(sample_metadata)

files = file.path(sample_metadata$path, "quant.sf")
names(files) <- rownames(sample_metadata)

# Import Salmon quant files
txi <- tximport(files, type='salmon', countsFromAbundance = 'no', txOut=TRUE)

# Save results
base::saveRDS(                       
  object = txi,  # The txi object
  file = snakemake@output[["txi"]])