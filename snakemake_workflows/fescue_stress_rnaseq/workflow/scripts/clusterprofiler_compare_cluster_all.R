# Purpose: Script to run compareCluster for all the comparisons and plot the results
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

# Define directories
lists_dir <- snakemake@input[["lists_dir"]]

# Load the dds object
dds <- readRDS(snakemake@input[['dds']])

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
if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
}

# Function to read rds file and add genotype name to the list name and 
# combine the list of dataframes into one dataframe and run compareCluster
run_compareCluster_combine <- function(dir_list, dds, term_to_gene, term_to_name) {

    signif_de_genes_list_shrink_combined <- list()

    for (dir in dir_list) {

        print(dir)

        # Get the list of rds files in the directory
        signif_de_genes_list_shrink <- list.files(dir, pattern = "signif_de_genes_list.shrink.*", full.names = TRUE)

        # Read rds file
        signif_de_genes_list_shrink <- readRDS(signif_de_genes_list_shrink)

        print(signif_de_genes_list_shrink)

        # Get genotype name from directory name and add to the list name
        # List name should be like HR2.320, HR12.320, D1.320, D3.320, W1.320
        # Genotype name is the string after last "_" in the directory name
        genotype <- strsplit(dir, "_")[[1]][length(strsplit(dir, "_")[[1]])]

        # Add genotype name to the list name
        names(signif_de_genes_list_shrink) <- paste0(genotype, ".", names(signif_de_genes_list_shrink))

        print(signif_de_genes_list_shrink)

        # Combine the list of dataframes into one dataframe
        signif_de_genes_list_shrink_combined <- c(signif_de_genes_list_shrink_combined, signif_de_genes_list_shrink)
    }

    # Create a dataframe from the named list
    signif_de_genes_list_shrink_combined_df <- do.call(rbind, lapply(names(signif_de_genes_list_shrink_combined), function(name) {
      components <- strsplit(name, split = ".", fixed = TRUE)[[1]]
      genotype <- components[1]
      time_point <- components[2]
    
      data.frame(
        gene = signif_de_genes_list_shrink_combined[[name]],
        time_point = time_point,
        genotype = genotype,
        stringsAsFactors = FALSE
      )
    }))

    compareCluster.combined <- compareCluster(gene~time_point+genotype ,data=signif_de_genes_list_shrink_combined_df, 
                     fun = enricher, 
                     TERM2GENE=term_to_gene,
                     TERM2NAME=term_to_name,
                     universe=rownames(counts(dds)),
                     minGSSize=minGSSize)
    compareCluster.combined.nona <- dplyr::filter(compareCluster.combined, Description != "NA")

    return(compareCluster.combined.nona)
}

# Run the function for comparison_heat_inf_vs_uninf
comparison_heat_inf_vs_uninf_dir_list <- lists_dir[grep("comparison_heat_inf_vs_uninf", lists_dir)]
print(comparison_heat_inf_vs_uninf_dir_list)

compareCluster.comparison_heat_inf_vs_uninf.combined.nona <- run_compareCluster_combine(comparison_heat_inf_vs_uninf_dir_list, dds, term_to_gene, term_to_name)

# Plot the results - Dot plot
compareCluster.comparison_heat_inf_vs_uninf.combined.nona.dotplot <- dotplot(compareCluster.comparison_heat_inf_vs_uninf.combined.nona, x="time_point", 
                                                    showCategory=n_category, label_format=60) + 
                                                    facet_grid(~genotype) + 
                                                    aes(x=fct_relevel(time_point, c('HR2', 'HR12', 'D1', 'D3', "W1"))) +
                                                    xlab(NULL) +
                                                    scale_color_gradientn(colours=c("#371ea3", "#46bac2", "#b3eebe"),
                                                    guide=guide_colorbar(reverse=TRUE, order=1)) +
                                                    guides(size = guide_legend(override.aes=list(shape=1))) +
                                                    theme(panel.grid.major.y = element_line(linetype='dotted', color='#808080'),
                                                    panel.grid.major.x = element_blank(),
                                                    legend.title=element_text(size=18), 
                                                    legend.text=element_text(size=16),
                                                    axis.text.y = element_text(size = 18),
                                                    axis.text.x = element_text(size = 18),
                                                    axis.title.x = element_text(size = 18),
                                                    strip.text.x = element_text(size = 20))

# Plot the results - Cnet plot
compareCluster.comparison_heat_inf_vs_uninf.combined.nona.cnetplot <- clusterProfiler::cnetplot(compareCluster.comparison_heat_inf_vs_uninf.combined.nona,
                            showCategory=15,
    						circular = FALSE,
                            cex.params = list(gene_label = 0.6, category_label = 2.5, gene_node=1, category_node=2),
                            node_label="category",
                            categorySize="pvalue",
                            layout = "nicely",
    						pie = "count",
    						legend_n = 2) +
                            labs(title = "",subtitle = "") +
                            theme(plot.title = element_text(size = 16, face = "bold"),
                                    plot.subtitle = element_text(size = 14),
                                    plot.caption = element_text(size = 12, face = "italic"),
			  	                    legend.title=element_text(size=18), 
                                    legend.text=element_text(size=16))

# Save the plot for comparison_heat_inf_vs_uninf
if (!dir.exists(paste0(plot_dir, "/comparison_heat_inf_vs_uninf/"))) {
    dir.create(paste0(plot_dir, "/comparison_heat_inf_vs_uninf/"), recursive = TRUE, showWarnings = FALSE)
}

ggsave(paste0(plot_dir, "/comparison_heat_inf_vs_uninf/", "compareCluster_all.comparison_heat_inf_vs_uninf", ".dotplot.png"), 
                compareCluster.comparison_heat_inf_vs_uninf.combined.nona.dotplot, 
                width = 22, height = 18, dpi=350)
ggsave(paste0(plot_dir, "/comparison_heat_inf_vs_uninf/", "compareCluster_all.comparison_heat_inf_vs_uninf", ".cnetplot.png"), 
                compareCluster.comparison_heat_inf_vs_uninf.combined.nona.cnetplot, 
                width = 22, height = 18, dpi=350)


# Run the function for comparison_infected_heat_vs_normal
comparison_infected_heat_vs_normal_dir_list <- lists_dir[grep("comparison_infected_heat_vs_normal", lists_dir)]
print(comparison_infected_heat_vs_normal_dir_list)

compareCluster.comparison_infected_heat_vs_normal.combined.nona <- run_compareCluster_combine(comparison_infected_heat_vs_normal_dir_list, dds, term_to_gene, term_to_name)

# Plot the results - Dot plot
compareCluster.comparison_infected_heat_vs_normal.combined.nona.dotplot <- dotplot(compareCluster.comparison_infected_heat_vs_normal.combined.nona, x="time_point", 
                                                    showCategory=n_category, label_format=60) + 
                                                    facet_grid(~genotype) + 
                                                    aes(x=fct_relevel(time_point, c('HR2', 'HR12', 'D1', 'D3', "W1"))) +
                                                    xlab(NULL) +
                                                    scale_color_gradientn(colours=c("#371ea3", "#46bac2", "#b3eebe"),
                                                    guide=guide_colorbar(reverse=TRUE, order=1)) +
                                                    guides(size = guide_legend(override.aes=list(shape=1))) +
                                                    theme(panel.grid.major.y = element_line(linetype='dotted', color='#808080'),
                                                    panel.grid.major.x = element_blank(),
                                                    legend.title=element_text(size=18), 
                                                    legend.text=element_text(size=16),
                                                    axis.text.y = element_text(size = 18),
                                                    axis.text.x = element_text(size = 18),
                                                    axis.title.x = element_text(size = 18),
                                                    strip.text.x = element_text(size = 20))

# Plot the results - Cnet plot
compareCluster.comparison_infected_heat_vs_normal.combined.nona.cnetplot <- clusterProfiler::cnetplot(compareCluster.comparison_infected_heat_vs_normal.combined.nona,
                            showCategory=15,
    						circular = FALSE,
                            cex.params = list(gene_label = 0.6, category_label = 2.5, gene_node=1, category_node=2),
                            node_label="category",
                            categorySize="pvalue",
                            layout = "nicely",
    						pie = "count",
    						legend_n = 2) +
                            labs(title = "",subtitle = "") +
                            theme(plot.title = element_text(size = 16, face = "bold"),
                                    plot.subtitle = element_text(size = 14),
                                    plot.caption = element_text(size = 12, face = "italic"),
			  	                    legend.title=element_text(size=18), 
                                    legend.text=element_text(size=16))

# Save the plot for comparison_infected_heat_vs_normal
if (!dir.exists(paste0(plot_dir, "/comparison_infected_heat_vs_normal/"))) {
    dir.create(paste0(plot_dir, "/comparison_infected_heat_vs_normal/"), recursive = TRUE, showWarnings = FALSE)
}

ggsave(paste0(plot_dir, "/comparison_infected_heat_vs_normal/", "compareCluster_all.comparison_infected_heat_vs_normal", ".dotplot.png"), 
                compareCluster.comparison_infected_heat_vs_normal.combined.nona.dotplot, 
                width = 24, height = 32, dpi=350)
ggsave(paste0(plot_dir, "/comparison_infected_heat_vs_normal/", "compareCluster_all.comparison_infected_heat_vs_normal", ".cnetplot.png"), 
                compareCluster.comparison_infected_heat_vs_normal.combined.nona.cnetplot, 
                width = 24, height = 26, dpi=350)

# Run the function for comparison_normal_inf_vs_uninf
comparison_normal_inf_vs_uninf_dir_list <- lists_dir[grep("comparison_normal_inf_vs_uninf", lists_dir)]
print(comparison_normal_inf_vs_uninf_dir_list)

compareCluster.comparison_normal_inf_vs_uninf.combined.nona <- run_compareCluster_combine(comparison_normal_inf_vs_uninf_dir_list, dds, term_to_gene, term_to_name)

# Plot the results
compareCluster.comparison_normal_inf_vs_uninf.combined.nona.dotplot <- dotplot(compareCluster.comparison_normal_inf_vs_uninf.combined.nona, x="time_point", 
                                                    showCategory=n_category, label_format=60) + 
                                                    facet_grid(~genotype) + 
                                                    aes(x=fct_relevel(time_point, c('HR2', 'HR12', 'D1', 'D3', "W1"))) +
                                                    xlab(NULL) +
                                                    scale_color_gradientn(colours=c("#371ea3", "#46bac2", "#b3eebe"),
                                                    guide=guide_colorbar(reverse=TRUE, order=1)) +
                                                    guides(size = guide_legend(override.aes=list(shape=1))) +
                                                    theme(panel.grid.major.y = element_line(linetype='dotted', color='#808080'),
                                                    panel.grid.major.x = element_blank(),
                                                    legend.title=element_text(size=18), 
                                                    legend.text=element_text(size=16),
                                                    axis.text.y = element_text(size = 18),
                                                    axis.text.x = element_text(size = 18),
                                                    axis.title.x = element_text(size = 18),
                                                    strip.text.x = element_text(size = 20))

# Plot the results - Cnet plot
compareCluster.comparison_normal_inf_vs_uninf.combined.nona.cnetplot <- clusterProfiler::cnetplot(compareCluster.comparison_normal_inf_vs_uninf.combined.nona,
                            showCategory=15,
    						circular = FALSE,
                            cex.params = list(gene_label = 0.6, category_label = 2.5, gene_node=1, category_node=2),
                            node_label="category",
                            categorySize="pvalue",
                            layout = "nicely",
    						pie = "count",
    						legend_n = 2) +
                            labs(title = "",subtitle = "") +
                            theme(plot.title = element_text(size = 16, face = "bold"),
                                    plot.subtitle = element_text(size = 14),
                                    plot.caption = element_text(size = 12, face = "italic"),
			  	                    legend.title=element_text(size=18), 
                                    legend.text=element_text(size=16))


# Save the plot for comparison_normal_inf_vs_uninf
if (!dir.exists(paste0(plot_dir, "/comparison_normal_inf_vs_uninf/"))) {
    dir.create(paste0(plot_dir, "/comparison_normal_inf_vs_uninf/"), recursive = TRUE, showWarnings = FALSE)
}

ggsave(paste0(plot_dir, "/comparison_normal_inf_vs_uninf/", "compareCluster_all.comparison_normal_inf_vs_uninf", ".dotplot.png"), 
                compareCluster.comparison_normal_inf_vs_uninf.combined.nona.dotplot, 
                width = 22, height = 20, dpi=350)
ggsave(paste0(plot_dir, "/comparison_normal_inf_vs_uninf/", "compareCluster_all.comparison_normal_inf_vs_uninf", ".cnetplot.png"), 
                compareCluster.comparison_normal_inf_vs_uninf.combined.nona.cnetplot, 
                width = 22, height = 20, dpi=350)

# Run the function for comparison_uninfected_heat_vs_normal
comparison_uninfected_heat_vs_normal_dir_list <- lists_dir[grep("comparison_uninfected_heat_vs_normal", lists_dir)]
print(comparison_uninfected_heat_vs_normal_dir_list)

compareCluster.comparison_uninfected_heat_vs_normal.combined.nona <- run_compareCluster_combine(comparison_uninfected_heat_vs_normal_dir_list, dds, term_to_gene, term_to_name)

# Plot the results
compareCluster.comparison_uninfected_heat_vs_normal.combined.nona.dotplot <- 
dotplot(compareCluster.comparison_uninfected_heat_vs_normal.combined.nona, x="time_point", 
                                                    showCategory=n_category, label_format=60) + 
                                                    facet_grid(~genotype) + 
                                                    aes(x=fct_relevel(time_point, c('HR2', 'HR12', 'D1', 'D3', "W1"))) +
                                                    xlab(NULL) +
                                                    scale_color_gradientn(colours=c("#371ea3", "#46bac2", "#b3eebe"),
                                                    guide=guide_colorbar(reverse=TRUE, order=1)) +
                                                    guides(size = guide_legend(override.aes=list(shape=1))) +
                                                    theme(panel.grid.major.y = element_line(linetype='dotted', color='#808080'),
                                                    panel.grid.major.x = element_blank(),
                                                    legend.title=element_text(size=18), 
                                                    legend.text=element_text(size=16),
                                                    axis.text.y = element_text(size = 18),
                                                    axis.text.x = element_text(size = 18),
                                                    axis.title.x = element_text(size = 18),
                                                    strip.text.x = element_text(size = 20))

# Plot the results - Cnet plot
compareCluster.comparison_uninfected_heat_vs_normal.combined.nona.cnetplot <- clusterProfiler::cnetplot(compareCluster.comparison_uninfected_heat_vs_normal.combined.nona,
                            showCategory=15,
    						circular = FALSE,
                            cex.params = list(gene_label = 0.6, category_label = 2.5, gene_node=1, category_node=2),
                            node_label="category",
                            categorySize="pvalue",
                            layout = "nicely",
    						pie = "count",
    						legend_n = 2) +
                            labs(title = "",subtitle = "") +
                            theme(plot.title = element_text(size = 16, face = "bold"),
                                    plot.subtitle = element_text(size = 14),
                                    plot.caption = element_text(size = 12, face = "italic"),
			  	                    legend.title=element_text(size=18), 
                                    legend.text=element_text(size=16))

# Save the plot for comparison_uninfected_heat_vs_normal
if (!dir.exists(paste0(plot_dir, "/comparison_uninfected_heat_vs_normal/"))) {
    dir.create(paste0(plot_dir, "/comparison_uninfected_heat_vs_normal/"), recursive = TRUE, showWarnings = FALSE)
}

ggsave(paste0(plot_dir, "/comparison_uninfected_heat_vs_normal/", "compareCluster_all.comparison_uninfected_heat_vs_normal", ".dotplot.png"), 
                compareCluster.comparison_uninfected_heat_vs_normal.combined.nona.dotplot, 
                width = 24, height = 32, dpi=350)
ggsave(paste0(plot_dir, "/comparison_uninfected_heat_vs_normal/", "compareCluster_all.comparison_uninfected_heat_vs_normal", ".cnetplot.png"), 
                compareCluster.comparison_uninfected_heat_vs_normal.combined.nona.cnetplot,
                width = 24, height = 26, dpi=350)

# Close the log file
sink(type="message")
sink()