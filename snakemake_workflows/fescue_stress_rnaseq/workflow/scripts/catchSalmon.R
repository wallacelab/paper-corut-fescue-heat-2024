# Redirect output and messages to log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(edgeR)
library(tidyverse)

# Assign the list of salmon quantification files
quant_files <- snakemake@input[["salmon_quant"]]

print(quant_files)

# Check if quant files were found, if not, raise an error
if (length(quant_files) == 0) {
    stop("No quant.sf files were provided!")
}

catch <- catchSalmon(paths = quant_files)

sample_metadata <- read.table(
  snakemake@input[["sample_metadata"]],
  header = TRUE,
  row.names = "sample",
  check.names = FALSE
)

# Check sample consistency between catch and sample_metadata
if (!all(colnames(catch) == rownames(sample_metadata))) {
  message("Warning: samples in catch are not the same as in sample_metadata")
}

# Save the scaled count matrix as RDS file
saveRDS(catch, snakemake@output[[1]])

# Close the log sink
sink(type = "message")
sink()
