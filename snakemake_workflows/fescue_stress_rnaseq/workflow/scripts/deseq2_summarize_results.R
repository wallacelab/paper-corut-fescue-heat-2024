# Purpose: Summarize DESeq2 results for each analysis and save the results as a csv file
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type="message")

library(tidyverse)
library(DESeq2)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(ggsci)

# Load the dds object
dds <- readRDS(snakemake@input[['dds_subset']])

# Load transcript hits information
tcp_hits <- read_csv(snakemake@input[['tcp_hits']])

# Remove space from column names and replace it with underscore
colnames(tcp_hits) <- gsub(" ", "_", colnames(tcp_hits))

# Extract analysis name
analysis_name <- snakemake@params$analysis$name

# Set thresholds
padj_threshold <- snakemake@params[["padj_threshold"]]
lfc_threshold <- snakemake@params[["lfc_threshold"]]

# Define directories
plot_dir <- snakemake@output[["tables_dir"]]
if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
}

# For each result_name, load the corresponding res and res.shrink
result_names <- resultsNames(dds)[!resultsNames(dds) == "Intercept"]

for (res_name in result_names) {
    # Load DESeq2 results
    res_file_path <- paste0(snakemake@input[["rds_dir"]], "/res.", analysis_name, ".", res_name, ".rds")
    print(paste("Loading", res_file_path))
    res <- readRDS(res_file_path)

    # Load lfc shrinkage results
    res_shrink_file_path <- paste0(snakemake@input[["rds_dir"]], "/res.", analysis_name, ".", res_name, ".shrink.rds")
    print(paste("Loading", res_shrink_file_path))  # Corrected to res_shrink_file_path for accurate message
    res.shrink <- readRDS(res_shrink_file_path)

    # Load IHW results
    res_ihw_file_path <- paste0(snakemake@input[["rds_dir"]], "/res.", analysis_name, ".", res_name, ".ihw.rds")
    print(paste("Loading", res_ihw_file_path))
    res.ihw <- readRDS(res_ihw_file_path)

    # Convert DESeq2 results to data frames
    df_res <- as.data.frame(res)
    df_res_shrink <- as.data.frame(res.shrink)
    df_res_ihw <- as.data.frame(res.ihw )

    # Convert row names to a column
    df_res$id <- rownames(df_res)
    df_res_shrink$id <- rownames(df_res_shrink)
    df_res_ihw$id <- rownames(df_res_ihw)

    # Add suffixes to columns (except 'baseMean' if it's the same across all data frames)
    df_res_shrink <- df_res_shrink %>% rename_with(~paste0(., ".shrink"), c(-baseMean, -id, -pvalue))
    df_res_ihw <- df_res_ihw %>% rename_with(~paste0(., ".ihw"), c(-baseMean, -id, -pvalue))

    # Perform a full join to merge all data frames by '"id", "baseMean", and "pvalue" column (which contains row names)
    merged_df <- full_join(df_res, df_res_shrink, by = c("id", "baseMean", "pvalue"))
    merged_df <- full_join(merged_df, df_res_ihw, by = c("id", "baseMean", "pvalue"))

    # Restore the row names from the 'id' column
    rownames(merged_df) <- merged_df$id
    merged_df$id <- NULL

    # Drop  padj.shrink, log2FoldChange.ihw, lfcSE.ihw, stat.ihw (they are same as padj, log2FoldChange, lfcSE, stat)
    merged_df <- merged_df %>% select(-padj.shrink, -log2FoldChange.ihw, -lfcSE.ihw, -stat.ihw)

    # Add transcript hits information to merged results table
    merged_df <- merge(merged_df, tcp_hits, by.x = "row.names", by.y = "TFQuery", all.x = TRUE)
    rownames(merged_df) <- merged_df$Row.names
    merged_df$Row.names <- NULL

    # Filter out rows based on padj_threshold and lfc_threshold
    merged_df <- merged_df %>% filter(padj <= padj_threshold, abs(log2FoldChange) >= lfc_threshold)

    # Order rows by padj
    merged_df <- merged_df[order(merged_df$padj), ]

    # Save the merged_df as a csv file
    write.csv(merged_df, 
            paste0(snakemake@output[["tables_dir"]], "/res_summary.signif_degs.", analysis_name, ".", res_name, ".csv"), 
            row.names = TRUE)
}

sink(type="message")
sink()