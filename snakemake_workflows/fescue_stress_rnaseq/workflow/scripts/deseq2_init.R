# Purpose: Initialize DESeq2 object for differential expression analysis
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type="message")

# load libraries
library("tidyverse")
library("tximport")
library("DESeq2")

# set up parallelization
parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}

print("Reading in tximport object...")
# read in the tximport object
txi <- readRDS(snakemake@input[["txi"]])

print("Reading in sample metadata...")
# read in the sample metadata
col_data <- read.table(
  snakemake@input[["sample_metadata"]],
  header = TRUE,
  row.names = "sample",
  check.names = FALSE
)

# order by rownames
col_data <- col_data[order(row.names(col_data)), , drop = FALSE] 

print("Relevelling variables of interest...")
# properly set the base level to the configuration in config.yaml, avoiding
# the default behaviour of choosing the alphabetical minimum level
for (vof in names(snakemake@config[["params"]][["deseq2"]][["variables_of_interest"]])) {
  # Convert the column to factor
  col_data[[vof]] <- as.factor(col_data[[vof]])
  
  # Print diagnostics to check the column levels and the desired base level
  print(paste("Checking variable:", vof))
  print(paste("Levels in data:", toString(levels(col_data[[vof]]))))
  base_level <- snakemake@config[["params"]][["deseq2"]][["variables_of_interest"]][[vof]][["base_level"]]
  print(paste("Base level from config:", base_level))
  
  # Relevel
  if (!base_level %in% levels(col_data[[vof]])) {
    warning(paste("Base level", base_level, "not found in", vof))
    next
  }
  cat("Relevelling", vof, "with base level", base_level, "\n")
  col_data[[vof]] <- relevel(col_data[[vof]], base_level)
}

# Ensure correct ordering of rownames
if (!identical(rownames(col_data), colnames(txi$counts))) {
  stop("Order of rownames in col_data doesn't match colnames of txi$counts")
}
# set the rownames of the col_data to the colnames of the txi object
rownames(col_data) <- colnames(txi$counts)

print(col_data)

# read in the design formula
design_formula <- snakemake@config[["params"]][["deseq2"]][["model"]]
cat("\n", design_formula, "\n")

print("Creating DESeq2 object...")
dds <- DESeqDataSetFromMatrix(round(txi$counts), 
    col_data, 
    as.formula(design_formula)
)

# Perform minimal filtering if activated
if(snakemake@config[["params"]][["deseq2"]][["filter"]][["minimal_filtering"]][["activate"]]) {
    
    print("Performing minimal filtering...")
    
    # extract min_counts from config file
    min_counts = snakemake@config[["params"]][["deseq2"]][["filter"]][["minimal_filtering"]][["min_counts"]]
    
    # remove uninformative columns based on min_counts
    dds <- dds[rowSums(counts(dds)) > min_counts, ]
    
    print(paste0("Number of genes after minimal filtering: ", nrow(dds)))
}

# Perform additional filtering if activated
if(snakemake@config[["params"]][["deseq2"]][["filter"]][["additional_filtering"]][["activate"]]) {
    
    print("Performing additional filtering...")
    
    # extract min_samples and min_counts_per_sample from config file
    min_samples = snakemake@config[["params"]][["deseq2"]][["filter"]][["additional_filtering"]][["min_samples"]]
    min_counts_per_sample = snakemake@config[["params"]][["deseq2"]][["filter"]][["additional_filtering"]][["min_counts_per_sample"]]
    
    # Filter out genes where there are less than 4 samples with counts greater than or equal to 10.
    keep <- rowSums(counts(dds) >= min_counts_per_sample) >= min_samples
    dds <- dds[keep,]
    
    print(paste0("Number of genes after additional filtering: ", nrow(dds)))
}

# Ensure non-empty dataset
if (nrow(dds) == 0 || ncol(dds) == 0) {
  stop("Data is empty after filtering")
}

# filter out outlier samples if specified
if(!is.null(snakemake@params[["exclude"]])) {
    samples_to_exclude <- snakemake@params[["exclude"]]
    print("Filtering out samples...")
    # Get the samples to keep
    samples_to_keep <- setdiff(colnames(dds), samples_to_exclude)
    # Subset the DESeqDataSet to include only the samples to keep
    dds <- dds[, samples_to_keep] 
}

print(colData(dds))

print("Filtering is complete...")

print("")

# subset deseq2 object if specified
# Check if subsetting is activated
if(snakemake@config[["params"]][["deseq2"]][["subset"]][["activate"]]){
    print("Subsetting deseq2 object...")
    # extract subset variables from config file
    subset_vars <- snakemake@config[["params"]][["deseq2"]][["subset"]][["filter"]]

    # Iterate over the subset conditions
    for(var_name in names(subset_vars)) {
      # Subset the DESeqDataSet
      dds <- dds[, dds[[var_name]] == subset_vars[[var_name]]]

      # Subset the sample_metadata
      colData(dds) <- colData(dds)[colData(dds)[[var_name]] == subset_vars[[var_name]], ]
    }
}

print("")
print(colData(dds))
print("Subsetting is complete...")

print("Grouping variables of interest...")
# Group columns if specified
if (snakemake@config[["params"]][["deseq2"]][["grouping"]][["activate"]]) {
    
    # extract column names from the config file
    grouping_vars <- snakemake@config[["params"]][["deseq2"]][["grouping"]][["vars"]]
    print(grouping_vars)
    # Use paste() function with do.call() to group columns
    colData(dds)$group <- do.call(paste, c(lapply(grouping_vars, function(x) colData(dds)[,x]), sep = ""))
    # Convert group column to factor
    colData(dds)$group <- factor(colData(dds)$group)
    print(colData(dds))

    print("Changing design formula to ~group ...")
    design(dds) <- as.formula(snakemake@config[["params"]][["deseq2"]][["grouping"]][["model"]])
    # dds <- DESeq(dds)
}

print(colData(dds))
print("Grouping is complete...")

print("Performing DESeq2 normalization and preprocessing...")
print(design(dds))

# normalization and preprocessing
dds <- DESeq(dds, parallel = parallel)

print("Saving DESeq2 object...")
# Write dds object as RDS
saveRDS(dds, file = snakemake@output[[1]])

# Close logging
sink(type="message")
sink()