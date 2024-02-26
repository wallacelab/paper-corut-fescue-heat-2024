# Purpose: Run go enrichment on the unique DEGs between two comparisons
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type="message")

library(tidyverse)
library(DESeq2)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(ggsci)
library(clusterProfiler)
library(enrichplot)
library(DOSE)

# Source the functions script
source(snakemake@params[['funcs_script']])

# Define directories
lists_dir <- snakemake@input[["lists_dir"]]

# Load the dds object
dds <- readRDS(snakemake@input[['dds']])

comparisons_of_interest <- snakemake@params[["comparisons_of_interest"]]

# Print the first element of the comparisons_of_interest list
first_comparison <- comparisons_of_interest[[1]]
print(comparisons_of_interest[[1]])

second_comparison <- comparisons_of_interest[[2]]
print(comparisons_of_interest[[2]])

# Load the term_to_gene and term_to_name objects
term_to_gene <- readRDS(snakemake@input[["term_to_gene"]])
term_to_name <- readRDS(snakemake@input[["term_to_name"]])

# Set thresholds
padj_threshold <- snakemake@params[["padj_threshold"]]
lfc_threshold <- snakemake@params[["lfc_threshold"]]
p_cutoff <- snakemake@params[["p_cutoff"]]
q_cutoff <- snakemake@params[["q_cutoff"]]
minGSSize <- snakemake@params[["minGSSize"]]
n_category <- snakemake@params[["n_category"]]

# Define directories
plot_dir <- snakemake@output[["plots_dir"]]
print(plot_dir)
if (!dir.exists(plot_dir)) {
    print("Creating plot directory")
    dir.create(plot_dir, recursive = TRUE, showWarnings = TRUE)
}

run_signif_de_genes_list_combine <- function(dir_list, dds, term_to_gene, term_to_name) {

    signif_de_genes_list_shrink_combined <- list()

    for (dir in dir_list) {

        print(dir)

        # Get the list of rds files in the directory
        signif_de_genes_list_shrink <- list.files(dir, pattern = "signif_de_genes_list.shrink.*", full.names = TRUE)

        # Read rds file
        signif_de_genes_list_shrink <- readRDS(signif_de_genes_list_shrink)

        # Get genotype name from directory name and add to the list name
        # List name should be like H2.320, H12.320, D1.320, D3.320, W1.320
        # Genotype name is the string after last "_" in the directory name
        genotype <- strsplit(dir, "_")[[1]][length(strsplit(dir, "_")[[1]])]

        # Add genotype name to the list name
        names(signif_de_genes_list_shrink) <- paste0(genotype, ".", names(signif_de_genes_list_shrink))

        print(signif_de_genes_list_shrink)

        # Combine the list of dataframes into one dataframe
        signif_de_genes_list_shrink_combined <- c(signif_de_genes_list_shrink_combined, signif_de_genes_list_shrink)
    }

    return (signif_de_genes_list_shrink_combined)
}

# Create a function to compare the heat and normal lists for one time point and genotype combination
get_unique_DEGs <- function(list_1, list_2, time_point, genotype) {
  # Get the heat and normal lists for the time point and genotype combination
  list_of_interest <- list_1[[paste(genotype, time_point, sep=".")]]

  list_to_compare <- list_2[[paste(genotype, time_point, sep=".")]]

  # Get the genes that are in heat but not normal
  diff_degs <- setdiff(list_of_interest, list_to_compare)
  
  # Check if heat.not.normal is empty
  if (length(diff_degs) == 0) {
    print(paste("No genes in", time_point, genotype, sep=" "))
    return(data.frame(gene=NA, time_point=time_point, genotype=genotype, stringsAsFactors=FALSE))
  } else {
    # Create a dataframe with the genes that are in heat but not normal
    return(data.frame(gene=diff_degs, time_point=time_point, genotype=genotype, stringsAsFactors=FALSE))
  }
}

# Get the lists for the first and second comparisons
first_comparison_list <- lists_dir[grep(first_comparison, lists_dir)]
second_comparison_list <- lists_dir[grep(second_comparison, lists_dir)]

# Run the function to combine the lists
first_comparison_list_combine <- run_signif_de_genes_list_combine(first_comparison_list, dds, term_to_gene, term_to_name)
print(first_comparison_list_combine)

# Run the function to combine the lists
second_comparison_list_combine <- run_signif_de_genes_list_combine(second_comparison_list, dds, term_to_gene, term_to_name)
print(second_comparison_list_combine)

# Get the time points and genotypes
time.points <- c("HR2", "HR12", "D1", "D3", "W1")
genotypes <- c("303", "315", "316", "320")

# Use lapply to apply the function to all time point and genotype combinations
first_comparison_unique_degs <- lapply(time.points, function(time_point) {
  lapply(genotypes, function(genotype) {
    get_unique_DEGs(first_comparison_list_combine, second_comparison_list_combine, time_point, genotype)
  })
})

# Flatten the list of lists into a single list of dataframes
first_comparison_unique_degs <- unlist(first_comparison_unique_degs, recursive = FALSE)


# Combine the list of dataframes into one dataframe
first_comparison_unique_degs_df <- do.call(rbind, first_comparison_unique_degs)


# Convert to dataframe if it's not already
first_comparison_unique_degs_df  <- as.data.frame(first_comparison_unique_degs_df )

print(str(first_comparison_unique_degs_df))

# Remove rows where gene is NA
# first_comparison_unique_degs_df <- subset(first_comparison_unique_degs_df , !is.na(first_comparison_unique_degs_df $gene))

# Get unique genotype counts
unique_genotype_counts <- first_comparison_unique_degs_df %>%
  group_by(gene) %>%
  summarize(unique_genotype_count = n_distinct(genotype))

# Merge this back with the original dataframe if needed
first_comparison_unique_degs_df <- merge(first_comparison_unique_degs_df, unique_genotype_counts, by = "gene", all.x = TRUE)
print(first_comparison_unique_degs_df)

# filter out genes that are in 2 or fewer genotypes
first_comparison_unique_degs_df_filtered <- first_comparison_unique_degs_df[first_comparison_unique_degs_df$unique_genotype_count > 2, ]

# order by time point
first_comparison_unique_degs_df_filtered <- first_comparison_unique_degs_df_filtered[order(first_comparison_unique_degs_df_filtered$time_point), ]

print(first_comparison_unique_degs_df_filtered[order(first_comparison_unique_degs_df_filtered$gene), ])
# get gene names and run go enrichment
gene.names <- first_comparison_unique_degs_df_filtered$gene
print(gene.names)

# Run go enrichment using enricher function from clusterprofiler
first_comparison_unique_degs_df_enrich <- enricher(gene.names, 
                     TERM2GENE=term_to_gene,
                     TERM2NAME=term_to_name,
                     universe=rownames(counts(dds)),
                     minGSSize=minGSSize)
print(first_comparison_unique_degs_df_enrich)
# Plot
first_comparison_unique_degs_df_enrich_lollipop <- generate_GO_enrich_lollipop(first_comparison_unique_degs_df_enrich, 
                                                                             n_category=n_category, font_size=14, 
                                                                             label_format=45, title="")
print(first_comparison_unique_degs_df_enrich_lollipop)
# Save the plot
ggsave(paste0(plot_dir, "/DE_in.", first_comparison, ".notDE_in.", second_comparison, ".enrich.lolipop.png"), 
    first_comparison_unique_degs_df_enrich_lollipop, width=12, height=10, dpi=350)

# Close the log file
sink(type="message")
sink()