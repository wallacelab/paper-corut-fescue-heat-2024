# Purpose: This script is used to perform the compareCluster function from clusterProfiler. 
# It takes as input the results of the DESeq2 analysis and the term_to_gene and term_to_name objects. 
# It generates the compareCluster plot and saves it in the plots directory. 
# It also saves the list of significant genes in the lists directory. 
# The list of significant genes is used by the clusterprofiler_enrichment.R script to generate the enrichment plots.
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

# Load the dds object
print(paste("Loading", snakemake@input[['dds_subset']]))
dds <- readRDS(snakemake@input[['dds_subset']])

# Extract analysis name
analysis_name <- snakemake@params$analysis$name
coef <- snakemake@params$analysis$coef_of_interest

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

lists_dir <- snakemake@output[["lists_dir"]]
if (!dir.exists(lists_dir)) {
    dir.create(lists_dir, recursive = TRUE, showWarnings = FALSE)
}

# For each result_name, load the corresponding res and res.shrink
result_names <- resultsNames(dds)[!resultsNames(dds) == "Intercept"]

result_names <- result_names[grep(coef, result_names)]

signif_de_genes_list <- list()
signif_de_genes_list_shrink <- list()

for (res_name in result_names) {
    res_file_path <- paste0(snakemake@input[["rds_dir"]], "/res.", analysis_name, ".", res_name, ".rds")
    print(paste("Loading", res_file_path))
    res <- readRDS(res_file_path)

    # Get significant results
    res.signif <- res[res$padj < padj_threshold & !is.na(res$padj) & abs(res$log2FoldChange) >= lfc_threshold, ]
    
    # Get signif gene names
    res.signif.genes <- rownames(res.signif)

    res_shrink_file_path <- paste0(snakemake@input[["rds_dir"]], "/res.", analysis_name, ".", res_name, ".shrink.rds")
    print(paste("Loading", res_shrink_file_path))  # Corrected to res_shrink_file_path for accurate message
    res.shrink <- readRDS(res_shrink_file_path)

    # Get significant results - shrink
    res.shrink.signif <- res.shrink[res.shrink$padj < padj_threshold & !is.na(res.shrink$padj) & abs(res.shrink$log2FoldChange) >= lfc_threshold, ]

    # Get signif gene names - shrink
    res.shrink.signif.genes <- rownames(res.shrink.signif)

    # Determine the position for insertion based on the specified conditions
    # The initial position is after the last element currently in the list
    insert_pos <- length(signif_de_genes_list) + 1

    # Adjust the position based on the specified ordering rules
    if (grepl("HR2", res_name) && !grepl("HR12", res_name)) {
        insert_pos <- 1
    } else if (grepl("HR12", res_name)) {
        # If there is already an "HR2" element, "HR12" should follow it
        if (any(grepl("HR2", names(signif_de_genes_list)))) {
            insert_pos <- which(grepl("HR2", names(signif_de_genes_list))) + 1
        } else {
            insert_pos <- 1
        }
    } else if (grepl("D1", res_name)) {
        # If there is already an "HR12" element, "D1" should follow it
        if (any(grepl("HR12", names(signif_de_genes_list)))) {
            insert_pos <- which(grepl("HR12", names(signif_de_genes_list))) + 1
        } else if (any(grepl("HR2", names(signif_de_genes_list)))) {
            insert_pos <- which(grepl("HR2", names(signif_de_genes_list))) + 1
        } else {
            insert_pos <- 1
        }
    } else if (grepl("D3", res_name)) {
        # Similar logic for "D3", it should follow "D1" if it exists
        if (any(grepl("D1", names(signif_de_genes_list)))) {
            insert_pos <- which(grepl("D1", names(signif_de_genes_list))) + 1
        } else if (any(grepl("HR12", names(signif_de_genes_list)))) {
            insert_pos <- which(grepl("HR12", names(signif_de_genes_list))) + 1
        } else if (any(grepl("HR2", names(signif_de_genes_list)))) {
            insert_pos <- which(grepl("HR2", names(signif_de_genes_list))) + 1
        } else {
            insert_pos <- 1
        }
    } else if (grepl("W1", res_name)) {
        # "W1" should always be last, so its position is correct by default
    }

    # Now insert the elements at the determined position
    signif_de_genes_list <- append(signif_de_genes_list, setNames(list(res.signif.genes), res_name), insert_pos - 1)
    signif_de_genes_list_shrink <- append(signif_de_genes_list_shrink, setNames(list(res.shrink.signif.genes), res_name), insert_pos - 1)
    
    print(names(signif_de_genes_list))
    print(names(signif_de_genes_list_shrink))
    
    # Change the names of the results based on the specified conditions. For example, if name contains HR2, change it to HR2
    # If name contains HR12, change it to HR12
    # If name contains D1, change it to D1
    # If name contains D3, change it to D3
    # If name contains W1, change it to W1
    if (grepl("HR2", res_name)) {
        names(signif_de_genes_list)[names(signif_de_genes_list) == res_name] <- "HR2"
        names(signif_de_genes_list_shrink)[names(signif_de_genes_list_shrink) == res_name] <- "HR2"
    } else if (grepl("HR12", res_name)) {
        names(signif_de_genes_list)[names(signif_de_genes_list) == res_name] <- "HR12"
        names(signif_de_genes_list_shrink)[names(signif_de_genes_list_shrink) == res_name] <- "HR12"
    } else if (grepl("D1", res_name)) {
        names(signif_de_genes_list)[names(signif_de_genes_list) == res_name] <- "D1"
        names(signif_de_genes_list_shrink)[names(signif_de_genes_list_shrink) == res_name] <- "D1"
    } else if (grepl("D3", res_name)) {
        names(signif_de_genes_list)[names(signif_de_genes_list) == res_name] <- "D3"
        names(signif_de_genes_list_shrink)[names(signif_de_genes_list_shrink) == res_name] <- "D3"
    } else if (grepl("W1", res_name)) {
        names(signif_de_genes_list)[names(signif_de_genes_list) == res_name] <- "W1"
        names(signif_de_genes_list_shrink)[names(signif_de_genes_list_shrink) == res_name] <- "W1"
    }

    print(names(signif_de_genes_list))
    print(names(signif_de_genes_list_shrink))
}

print(signif_de_genes_list_shrink)

# Save the list of significant genes
saveRDS(signif_de_genes_list_shrink, paste0(lists_dir, "/", "signif_de_genes_list.shrink.", analysis_name, ".rds"))
saveRDS(signif_de_genes_list, paste0(lists_dir, "/", "signif_de_genes_list.", analysis_name, ".rds"))

# Run compareCluster
compareCluster.res <- compareCluster(geneCluster=signif_de_genes_list_shrink, 
                                              fun = enricher,
                                              universe = rownames(counts(dds)), 
                                              TERM2GENE=term_to_gene,
                                              TERM2NAME=term_to_name)

# Check if compareCluster.res is empty. If empty then print a message and exit
if (length(compareCluster.res) == 0) {
    print("compareCluster.res is empty")
    q()
}

# Remove rows with NA in Description column
compareCluster.res.nona <- dplyr::filter(compareCluster.res, Description != "NA")

# Save the compareCluster results as RDS
saveRDS(compareCluster.res.nona, paste0(lists_dir, "/", "compareCluster.res.nona.", analysis_name, ".GO.rds"))

print("Generating the dotplot...")

if (analysis_name=="uninfected_heat_vs_normal_303"){
    # Plot the results
    compareCluster.res.nona.dotplot <- clusterProfiler::dotplot(compareCluster.res.nona, showCategory=15, label_format=45) + 
                                xlab(NULL) +
                                # scale_color_gradientn(colours=c("#371ea3", "#46bac2", "#b3eebe"),
                                # guide=guide_colorbar(reverse=TRUE, order=1)) +
                                guides(size = guide_legend(override.aes=list(shape=1))) +
                                theme(panel.grid.major.y = element_line(linetype='dotted', color='#808080'),
                                panel.grid.major.x = element_blank(),
                                plot.title = element_text(size = 16, face = "bold"),
                                plot.caption = element_text(size = 12, face = "italic"),
                                legend.title=element_text(size=22), 
                                legend.text=element_text(size=20),
                                axis.text.y = element_text(size = 23),
                                axis.text.x = element_text(size = 23),
                                axis.title.x = element_text(size = 22)) +
                                scale_colour_gradientn(colours = pal_npg()(4))

    print("Saving the dotplot...")
    # Save the plot
    ggsave(paste0(plot_dir, "/", "compareCluster.", analysis_name, ".dotplot.png"), compareCluster.res.nona.dotplot, width=14.5, height=18, dpi=350)
} else {
   # Plot the results
    compareCluster.res.nona.dotplot <- clusterProfiler::dotplot(compareCluster.res.nona, showCategory=15, label_format=40) + 
                                xlab(NULL) +
                                # scale_color_gradientn(colours=c("#371ea3", "#46bac2", "#b3eebe"),
                                # guide=guide_colorbar(reverse=TRUE, order=1)) +
                                guides(size = guide_legend(override.aes=list(shape=1))) +
                                theme(panel.grid.major.y = element_line(linetype='dotted', color='#808080'),
                                panel.grid.major.x = element_blank(),
                                plot.title = element_text(size = 16, face = "bold"),
                                plot.caption = element_text(size = 12, face = "italic"),
                                legend.title=element_text(size=22), 
                                legend.text=element_text(size=20),
                                axis.text.y = element_text(size = 23),
                                axis.text.x = element_text(size = 23),
                                axis.title.x = element_text(size = 22)) +
                                scale_colour_gradientn(colours = pal_npg()(4))

    print("Saving the dotplot...")
    # Save the plot
    ggsave(paste0(plot_dir, "/", "compareCluster.", analysis_name, ".dotplot.png"), compareCluster.res.nona.dotplot, width=14.5, height=18, dpi=350)
}

print("Generating the cnetplot...")
png(paste0(plot_dir, "/", "compareCluster.", analysis_name, ".cnetplot.png"), width=20, height=14, units = "in", res=350)
compareCluster.res.nona.cnetplot <- clusterProfiler::cnetplot(compareCluster.res.nona,
                            showCategory=15,
    						circular = FALSE,
                            cex.params = list(gene_label = 0.8, category_label = 1, gene_node=1, category_node=1),
                            node_label="category",
                            categorySize="pvalue",
                            layout = "nicely",
    						pie = "count") +
                            labs(title = "",subtitle = "") +
                            theme(plot.title = element_text(size = 16, face = "bold"),
                                    plot.subtitle = element_text(size = 14),
                                    plot.caption = element_text(size = 12, face = "italic"),
			  	                    legend.title=element_text(size=18), 
                                    legend.text=element_text(size=16))

print("Saving the cnetplot...")
# Save the plot
# ggsave(paste0(plot_dir, "/", "compareCluster.", analysis_name, ".cnetplot.png"), compareCluster.res.nona.cnetplot, width=18, height=14, dpi=350)
dev.off()

print("Done!")

# Close the log file
sink(type="message")
sink()