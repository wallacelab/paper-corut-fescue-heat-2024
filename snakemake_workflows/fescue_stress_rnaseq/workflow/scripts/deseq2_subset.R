# Purpose: Subset the DESeq2 object based on the defined variables and run DESeq2
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type="message")

# load libraries
library("tidyverse")
library("DESeq2")

# set up parallelization
parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}

# Load the dds object
dds <- readRDS(snakemake@input[['dds']])

# Subset based on the defined variables
analysis <- snakemake@params$analysis

# Print levels of time_point before subsetting
print(paste("Original levels of time_point:", paste(levels(dds$time_point), collapse=", ")))
# If you're subsetting based on other criteria, print their levels as well:
print(paste("Original levels of condition:", paste(levels(dds$condition), collapse=", ")))
print(paste("Original levels of genotype:", paste(levels(dds$genotype), collapse=", ")))
print(paste("Original levels of treatment:", paste(levels(dds$treatment), collapse=", ")))

subset_conditions <- list()
for (variable in names(analysis$subset)) {
    subset_conditions <- c(subset_conditions, paste0("(dds$", variable, '=="', analysis$subset[[variable]], '")'))
}
subset_string <- paste(subset_conditions, collapse=" & ")
dds_subset <- dds[, eval(parse(text=subset_string))]
print(paste("Subsetted data based on", subset_string))

# Exclude specified time points
if (!is.null(analysis$exclude_time_points)) {
    exclude_conditions <- paste0("dds_subset$time_point %in% c(", paste(shQuote(analysis$exclude_time_points), collapse=","), ")")
    dds_subset <- dds_subset[, !eval(parse(text=exclude_conditions))]
    print(paste("Excluded time points:", paste(analysis$exclude_time_points, collapse=", ")))
}

# Drop unused levels based on subsetting and exclusion
for (variable in names(analysis$subset)) {
    if (!is.null(dds_subset[[variable]])) {
        dds_subset[[variable]] <- droplevels(dds_subset[[variable]])
        print(paste("Dropped unused levels for", variable))
    }
}
if (!is.null(analysis$exclude_time_points)) {
    dds_subset$time_point <- droplevels(dds_subset$time_point)
    print("Dropped unused levels for time_point")
}


if (snakemake@config[["params"]][["deseq2"]][["grouping"]][["activate"]]) {
    dds_subset$group <- droplevels(dds_subset$group)
}

print(paste("Levels of time_point after subsetting:", paste(levels(dds_subset$time_point), collapse=", ")))
print(paste("Levels of condition after subsetting:", paste(levels(dds_subset$condition), collapse=", ")))
print(paste("Levels of genotype after subsetting:", paste(levels(dds_subset$genotype), collapse=", ")))
print(paste("Levels of treatment after subsetting:", paste(levels(dds_subset$treatment), collapse=", ")))

  
# Set the design formula
design_formula <- as.formula(analysis$formula)
design(dds_subset) <- design_formula
print(design(dds_subset))

variables <- all.vars(design_formula)
for (variable in variables) {
    if (is.factor(dds_subset[[variable]])) {
        print(paste("Levels of", variable, ":", paste(levels(dds_subset[[variable]]), collapse=", ")))
    }
}

# Filter genes if filtering is enabled
if (analysis$filter$enable) {
  keep <- rowSums(counts(dds_subset) >= analysis$filter$min_counts) >= analysis$filter$min_samples
  dds_subset <- dds_subset[keep,]
}

# Run DESeq2
dds_subset <- DESeq(dds_subset, test="LRT", 
                     reduced = as.formula(analysis$reduced_formula), 
                     parallel = parallel)

# Variance Stabilizing Transformation
vsd <- vst(dds_subset, blind=FALSE)
  
# Save the vsd object
saveRDS(vsd, file= snakemake@output[["vsd_subset"]])

# Save the dds object after processing
saveRDS(dds_subset, file= snakemake@output[["dds_subset"]])

# Print results names for each analysis
cat(paste0("Results for analysis: ", analysis$name, "\n"))
print(resultsNames(dds_subset))


# Close logging
sink(type="message")
sink()