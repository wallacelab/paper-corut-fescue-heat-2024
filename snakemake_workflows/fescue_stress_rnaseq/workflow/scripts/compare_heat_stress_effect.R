# Purpose: Compare the effect of heat stress on gene expression in infected and uninfected plants at different time points.
# It finds the genes that are differentially expressed in infected and uninfected plants at each time point and then compares the genes that are 
# differentially expressed in infected and uninfected plants at each time point to find the genes that are unique to infected plants at each time point.
# It uses DETs that are shared in at least 3 comparisons in E+ plants.
# It then performs GO enrichment analysis on the unique genes to find the biological processes that are enriched in the unique genes.

# Open log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type="message")

library(tidyverse)
library(DESeq2)
library(VennDetail)
library(clusterProfiler)
library(enrichplot)
library(DOSE)
library(ggsci)

# Source the functions script
source(snakemake@params[['funcs_script']])

dds <- readRDS(snakemake@input[['dds']])
term_to_gene <- readRDS(snakemake@input[["term_to_gene"]])
term_to_name <- readRDS(snakemake@input[["term_to_name"]])

# Function to create file pattern
create_pattern <- function(base_pattern, time_point) {
  return(paste0(base_pattern, time_point, ".shrink.rds"))
}

# Function to generate names for list elements based on the pattern and file names
generate_names_for_list <- function(file_names, base_pattern, time_point) {
  modified_pattern <- gsub("\\.\\*", "|", base_pattern)
  names_for_list <- gsub(paste0(modified_pattern, time_point, "|.shrink.rds"), "", basename(file_names))
  names_for_list <- paste0(names_for_list, time_point)
  print(names_for_list)
  return(names_for_list)
}

# Function to plot Venn diagram
get_venn_diagram <- function(sig_genes, comparison, time_point) {
  ven <- venndetail(sig_genes)
  return(ven)
}

# Function to process files for a given time point
process_time_point <- function(time_point, base_pattern, dir, padj_threshold, lfc_threshold, comparison, min_shared) {
  # Get the files for the time point
  res_shrink_relative <- list.files(dir, 
                            pattern = create_pattern(base_pattern, time_point),
                            recursive = TRUE)
  
  # Get the full path of the files
  res_shrink <- file.path(dir, res_shrink_relative)
  # print(res_shrink)
  # Read the files
  res_shrink_list <- lapply(res_shrink, readRDS)

  # Get the names for the list elements
  names_for_list <- generate_names_for_list(res_shrink, base_pattern, time_point)
  
  # Ensure that the length of names matches the length of the list
  if(length(res_shrink_list) == length(names_for_list)) {
    # Set names for the list elements
    names(res_shrink_list) <- names_for_list
  } else {
    stop("Error: The number of elements and the number of names do not match.")
  }

  # Print the list
  # print(res_shrink_list)

  # Create a list to store significant genes
  sig_genes <- list()

  # for each element in res_shrink_list
  for (i in 1:length(res_shrink_list)) {
      # Get the name of the element
      name <- names(res_shrink_list[i])
      print(name)

      # Get the deseq2 results for that element
      res <- res_shrink_list[[i]]
      print(res)

      # Get significant results
      res_shrink_signif <- res[res$padj < padj_threshold & !is.na(res$padj) & abs(res$log2FoldChange) >= lfc_threshold, ]
      
      # Get the significant genes
      res_shrink_signif_genes <- rownames(res_shrink_signif)
      
      # Add the significant genes to the sig_genes list
      sig_genes[[name]] <- res_shrink_signif_genes
  }

  # Plot Venn diagram
  venn <- get_venn_diagram(sig_genes, comparison, time_point)

  # Find the genes shared in at least 3 comparisons
  print(head(venn@result))
  shared_genes <- getSet(venn, min=min_shared)$Detail
  print(comparison)
  print(time_point)
  print(paste0("Number of shared genes in at least ", min_shared, " genotypes: "))
  print(length(shared_genes))
  return(shared_genes)
}

# Main function
run_analysis <- function(dds_dir_list, results_dir_list, time_point, base_pattern, padj_threshold, lfc_threshold, comparison, min_shared) {
  dir <- dirname(results_dir_list[1])
  print(dir)

  # Process files for that time point
  shared_genes <- process_time_point(time_point, base_pattern, dir, padj_threshold, lfc_threshold, comparison, min_shared)
  return(shared_genes)
}

# Set parameters (could be read from a config file or snakemake parameters)
time_points <- c("HR2", "HR12", "D1", "D3", "W1")
padj_threshold <- snakemake@params[["padj_threshold"]]
lfc_threshold <- snakemake@params[["lfc_threshold"]]
p_cutoff <- snakemake@params[["p_cutoff"]]
q_cutoff <- snakemake@params[["q_cutoff"]]
n_category <- snakemake@params[["n_category"]]

res_dir <- snakemake@output[["res_dir"]]
if (!dir.exists(res_dir)) {
    print("Creating results directory")
    dir.create(res_dir, recursive = TRUE, showWarnings = TRUE)
}

plot_dir <- snakemake@output[["plots_dir"]]
if (!dir.exists(plot_dir)) {
    print("Creating plot directory")
    dir.create(plot_dir, recursive = TRUE, showWarnings = TRUE)
}

# Run analysis

# res.heat_uninf_vs_inf
base_pattern_infected_heat_vs_normal <- "res.infected_heat_vs_normal_.*conditionHS.time_point"
infected_heat_vs_normal.HR2.shared <- run_analysis(snakemake@input$dds_subset, snakemake@input$rds_dir, 
            time_point="HR2", base_pattern_infected_heat_vs_normal, 
            padj_threshold, lfc_threshold, 
            comparison= "infected_heat_vs_normal",
            min_shared=1)

infected_heat_vs_normal.HR12.shared <- run_analysis(snakemake@input$dds_subset, snakemake@input$rds_dir,
            time_point="HR12", base_pattern_infected_heat_vs_normal, 
            padj_threshold, lfc_threshold, 
            comparison= "infected_heat_vs_normal",
            min_shared=1)

infected_heat_vs_normal.D1.shared <- run_analysis(snakemake@input$dds_subset, snakemake@input$rds_dir,
            time_point="D1", base_pattern_infected_heat_vs_normal, 
            padj_threshold, lfc_threshold, 
            comparison= "infected_heat_vs_normal",
            min_shared=1)  

infected_heat_vs_normal.D3.shared <- run_analysis(snakemake@input$dds_subset, snakemake@input$rds_dir,
            time_point="D3", base_pattern_infected_heat_vs_normal, 
            padj_threshold, lfc_threshold, 
            comparison= "infected_heat_vs_normal",
            min_shared=1)

infected_heat_vs_normal.W1.shared <- run_analysis(snakemake@input$dds_subset, snakemake@input$rds_dir,
            time_point="W1", base_pattern_infected_heat_vs_normal, 
            padj_threshold, lfc_threshold, 
            comparison= "infected_heat_vs_normal",
            min_shared=1)
        

# res.uninfected_heat_vs_normal
base_pattern_uninfected_heat_vs_normal <- "res.uninfected_heat_vs_normal_.*conditionHS.time_point"
uninfected_heat_vs_normal.HR2.shared <- run_analysis(snakemake@input$dds_subset, snakemake@input$rds_dir, 
            time_point="HR2", base_pattern_uninfected_heat_vs_normal, 
            padj_threshold, lfc_threshold, 
            comparison= "uninfected_heat_vs_normal",
            min_shared=0)

uninfected_heat_vs_normal.HR12.shared <- run_analysis(snakemake@input$dds_subset, snakemake@input$rds_dir,
            time_point="HR12", base_pattern_uninfected_heat_vs_normal, 
            padj_threshold, lfc_threshold, 
            comparison= "uninfected_heat_vs_normal",
            min_shared=0)

uninfected_heat_vs_normal.D1.shared <- run_analysis(snakemake@input$dds_subset, snakemake@input$rds_dir,
            time_point="D1", base_pattern_uninfected_heat_vs_normal, 
            padj_threshold, lfc_threshold, 
            comparison= "uninfected_heat_vs_normal",
            min_shared=0)

uninfected_heat_vs_normal.D3.shared <- run_analysis(snakemake@input$dds_subset, snakemake@input$rds_dir,
            time_point="D3", base_pattern_uninfected_heat_vs_normal, 
            padj_threshold, lfc_threshold, 
            comparison= "uninfected_heat_vs_normal",
            min_shared=0)

uninfected_heat_vs_normal.W1.shared <- run_analysis(snakemake@input$dds_subset, snakemake@input$rds_dir,
            time_point="W1", base_pattern_uninfected_heat_vs_normal, 
            padj_threshold, lfc_threshold, 
            comparison= "uninfected_heat_vs_normal",
            min_shared=0)

print("")

# HR2
## Compare infected_heat_vs_normal.HR2.shared and uninfected_heat_vs_normal.HR2.shared and get unique genes in infected_heat_vs_normal.HR2.shared
infected_heat_vs_normal.HR2.unique <- setdiff(infected_heat_vs_normal.HR2.shared, uninfected_heat_vs_normal.HR2.shared)
print("HR2")
print("Number of E+ specific DETs: ")
print(length(infected_heat_vs_normal.HR2.unique))

## Run GO enrichment analysis
infected_heat_vs_normal.HR2.unique.GO <- clusterProfiler::enricher(infected_heat_vs_normal.HR2.unique, 
                     TERM2GENE=term_to_gene,
                     TERM2NAME=term_to_name,
                     universe=rownames(counts(dds)),
                     pvalueCutoff = p_cutoff,
                     pAdjustMethod = "BH",
                     qvalueCutoff = q_cutoff)
print(infected_heat_vs_normal.HR2.unique.GO)
## Save results
saveRDS(infected_heat_vs_normal.HR2.unique.GO, file.path(res_dir, "infected_heat_vs_normal.HR2.unique.GO.rds"))

## Generate dotplot
infected_heat_vs_normal.HR2.unique.GO.dotplot <- generate_GO_enrich_dotplot(infected_heat_vs_normal.HR2.unique.GO, 
                                                                             n_category=n_category, font_size=14, 
                                                                             label_format=45, title="")                                                                       
## Save dotplot
ggplot2::ggsave(paste0(plot_dir, "/infected_heat_vs_normal.HR2.unique.GO.dotplot.png"), 
                infected_heat_vs_normal.HR2.unique.GO.dotplot, width=12, height=14, dpi=350)

# HR12
## Compare infected_heat_vs_normal.HR12.shared and uninfected_heat_vs_normal.HR12.shared and get unique genes in infected_heat_vs_normal.HR12.shared
infected_heat_vs_normal.HR12.unique <- setdiff(infected_heat_vs_normal.HR12.shared, uninfected_heat_vs_normal.HR12.shared)
print("HR12")
print("Number of E+ specific DETs: ")
print(length(infected_heat_vs_normal.HR12.unique))

## Run GO enrichment analysis
infected_heat_vs_normal.HR12.unique.GO <- clusterProfiler::enricher(infected_heat_vs_normal.HR12.unique, 
                     TERM2GENE=term_to_gene,
                     TERM2NAME=term_to_name,
                     universe=rownames(counts(dds)),
                     pvalueCutoff = p_cutoff,
                     pAdjustMethod = "BH",
                     qvalueCutoff = q_cutoff)
print(infected_heat_vs_normal.HR12.unique.GO)
## Save results
saveRDS(infected_heat_vs_normal.HR12.unique.GO, file.path(res_dir, "infected_heat_vs_normal.HR12.unique.GO.rds"))

## Generate dotplot
infected_heat_vs_normal.HR12.unique.GO.dotplot <- generate_GO_enrich_dotplot(infected_heat_vs_normal.HR12.unique.GO, 
                                                                             n_category=n_category, font_size=14, 
                                                                             label_format=30, title="")
## Save dotplot
ggplot2::ggsave(paste0(plot_dir, "/infected_heat_vs_normal.HR12.unique.GO.dotplot.png"), 
                infected_heat_vs_normal.HR12.unique.GO.dotplot, width=10.5, height=11, dpi=350)

# D1
## Compare infected_heat_vs_normal.D1.shared and uninfected_heat_vs_normal.D1.shared and get unique genes in infected_heat_vs_normal.D1.shared
infected_heat_vs_normal.D1.unique <- setdiff(infected_heat_vs_normal.D1.shared, uninfected_heat_vs_normal.D1.shared)
print("D1")
print("Number of E+ specific DETs: ")
print(length(infected_heat_vs_normal.D1.unique))

## Run GO enrichment analysis
infected_heat_vs_normal.D1.unique.GO <- clusterProfiler::enricher(infected_heat_vs_normal.D1.unique, 
                     TERM2GENE=term_to_gene,
                     TERM2NAME=term_to_name,
                     universe=rownames(counts(dds)),
                     pvalueCutoff = p_cutoff,
                     pAdjustMethod = "BH",
                     qvalueCutoff = q_cutoff)
print(infected_heat_vs_normal.D1.unique.GO)
## Save results
saveRDS(infected_heat_vs_normal.D1.unique.GO, file.path(res_dir, "infected_heat_vs_normal.D1.unique.GO.rds"))
## Generate dotplot
infected_heat_vs_normal.D1.unique.GO.dotplot <- generate_GO_enrich_dotplot(infected_heat_vs_normal.D1.unique.GO, 
                                                                             n_category=n_category, font_size=14, 
                                                                             label_format=45, title="")

## Save dotplot
ggplot2::ggsave(paste0(plot_dir, "/infected_heat_vs_normal.D1.unique.GO.dotplot.png"), 
                infected_heat_vs_normal.D1.unique.GO.dotplot, width=12, height=14, dpi=350)
                
# D3
## Compare infected_heat_vs_normal.D3.shared and uninfected_heat_vs_normal.D3.shared and get unique genes in infected_heat_vs_normal.D3.shared
infected_heat_vs_normal.D3.unique <- setdiff(infected_heat_vs_normal.D3.shared, uninfected_heat_vs_normal.D3.shared)

print("D3")
print("Number of E+ specific DETs: ")
print(length(infected_heat_vs_normal.D3.unique))

## Run GO enrichment analysis
infected_heat_vs_normal.D3.unique.GO <- clusterProfiler::enricher(infected_heat_vs_normal.D3.unique, 
                     TERM2GENE=term_to_gene,
                     TERM2NAME=term_to_name,
                     universe=rownames(counts(dds)),
                     pvalueCutoff = p_cutoff,
                     pAdjustMethod = "BH",
                     qvalueCutoff = q_cutoff)
print(infected_heat_vs_normal.D3.unique.GO)
## Save results
saveRDS(infected_heat_vs_normal.D3.unique.GO, file.path(res_dir, "infected_heat_vs_normal.D3.unique.GO.rds"))


## Generate dotplot
infected_heat_vs_normal.D3.unique.GO.dotplot <- generate_GO_enrich_dotplot(infected_heat_vs_normal.D3.unique.GO, 
                                                                             n_category=n_category, font_size=14, 
                                                                             label_format=28, title="")

## Save dotplot
ggplot2::ggsave(paste0(plot_dir, "/infected_heat_vs_normal.D3.unique.GO.dotplot.png"), 
                infected_heat_vs_normal.D3.unique.GO.dotplot,  width=10.5, height=11, dpi=350)

# W1
## Compare infected_heat_vs_normal.W1.shared and uninfected_heat_vs_normal.W1.shared and get unique genes in infected_heat_vs_normal.W1.shared
infected_heat_vs_normal.W1.unique <- setdiff(infected_heat_vs_normal.W1.shared, uninfected_heat_vs_normal.W1.shared)

print("W1")
print("Number of E+ specific DETs: ")
print(length(infected_heat_vs_normal.W1.unique))

## Run GO enrichment analysis
infected_heat_vs_normal.W1.unique.GO <- clusterProfiler::enricher(infected_heat_vs_normal.W1.unique, 
                     TERM2GENE=term_to_gene,
                     TERM2NAME=term_to_name,
                     universe=rownames(counts(dds)),
                     pvalueCutoff = p_cutoff,
                     pAdjustMethod = "BH",
                     qvalueCutoff = q_cutoff)
print(infected_heat_vs_normal.W1.unique.GO)
## Save results
saveRDS(infected_heat_vs_normal.W1.unique.GO, file.path(res_dir, "infected_heat_vs_normal.W1.unique.GO.rds"))

## Generate dotplot
infected_heat_vs_normal.W1.unique.GO.dotplot <- generate_GO_enrich_dotplot(infected_heat_vs_normal.W1.unique.GO, 
                                                                             n_category=n_category, font_size=14, 
                                                                             label_format=35, title="")
## Save dotplot
ggplot2::ggsave(paste0(plot_dir, "/infected_heat_vs_normal.W1.unique.GO.dotplot.png"), 
                infected_heat_vs_normal.W1.unique.GO.dotplot, width=10.5, height=11, dpi=350)
                
# Close log file
sink()
sink(type = "message")

