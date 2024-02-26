# Purpose: Run DESeq2 differential expression analysis
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type="message")

library(DESeq2)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(pheatmap)
library(gridExtra)
library(ggpubr)
library(ggsci)
library(RColorBrewer)
library(reshape2)
library(EnhancedVolcano)
library(ashr)
library(DOSE)
library(ComplexHeatmap)
library(IHW)

# Source the functions script
source(snakemake@params[['funcs_script']])

# Load the dds object
dds <- readRDS(snakemake@input[['dds_subset']])
vsd <- readRDS(snakemake@input[['vsd_subset']])

# Extract analysis name
analysis_name <- gsub("dds_(.*).rds", "\\1", basename(snakemake@input[['dds_subset']]))

# Set thresholds
padj_threshold <- snakemake@params[["padj_threshold"]]
lfc_threshold <- snakemake@params[["lfc_threshold"]]
p_cutoff <- snakemake@params[["p_cutoff"]]
q_cutoff <- snakemake@params[["q_cutoff"]]
minGSSize <- snakemake@params[["minGSSize"]]

# Analysis
result_names <- resultsNames(dds)[!resultsNames(dds) == "Intercept"]

for (res_name in result_names) {
    # Extract the result for this specific comparison
    res <- results(dds, name=res_name, test="Wald", parallel=TRUE)

    # Create the directory if doesn't exist
    if (!dir.exists(snakemake@output[["tables_dir"]])) {
        dir.create(snakemake@output[["tables_dir"]], recursive = TRUE, showWarnings = FALSE)
    }
    if (!dir.exists(snakemake@output[["rds_dir"]])) {
        dir.create(snakemake@output[["rds_dir"]], recursive = TRUE, showWarnings = FALSE)
    }
    # Save the results
    write.csv(res, 
          file = paste0(snakemake@output[["tables_dir"]], "/res.", analysis_name, ".", res_name, ".csv"),
          row.names = TRUE)
    saveRDS(res, paste0(snakemake@output[["rds_dir"]], "/res.", analysis_name, ".", res_name, ".rds"))

    # Results with independent hypothesis weighting
    res.ihw <- results(dds, name=res_name, test="Wald", filterFun=ihw, parallel = TRUE)

    # Save the results
    write.csv(res.ihw,
          file = paste0(snakemake@output[["tables_dir"]], "/res.", analysis_name, ".", res_name, ".ihw.csv"),
          row.names = TRUE)
    saveRDS(res.ihw, paste0(snakemake@output[["rds_dir"]], "/res.", analysis_name, ".", res_name, ".ihw.rds"))

    # Perform log2 fold change shrinkage
    res.shrink <- lfcShrink(dds, coef=res_name, res=res, type="ashr", parallel = TRUE)

    # Save the results
    write.csv(res.shrink,
          file = paste0(snakemake@output[["tables_dir"]], "/res.", analysis_name, ".", res_name, ".shrink.csv"),
          row.names = TRUE)
    saveRDS(res.shrink, paste0(snakemake@output[["rds_dir"]], "/res.", analysis_name, ".", res_name, ".shrink.rds"))

}

# Close the log file
sink(type="message")
sink()
