# Purpose: Run GO enrichment analysis using clusterProfiler and generate plots
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type="message")

library(tidyverse)
library(DESeq2)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(ggsci)
library(RColorBrewer)
library(reshape2)
library(DOSE)
library(ComplexHeatmap)

# Source the functions script
source(snakemake@params[['funcs_script']])

# Load the dds object
dds <- readRDS(snakemake@input[['dds_subset']])

# Extract analysis name
analysis_name <- gsub("dds_(.*).rds", "\\1", basename(snakemake@input[['dds_subset']]))

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
if (!dir.exists(paste0(plot_dir, "/dotplot"))) {
    dir.create(paste0(plot_dir, "/dotplot"), recursive = TRUE, showWarnings = FALSE)
}
if (!dir.exists(paste0(plot_dir, "/barplot"))) {
    dir.create(paste0(plot_dir, "/barplot"), recursive = TRUE, showWarnings = FALSE)
}
if (!dir.exists(paste0(plot_dir, "/lollipop"))) {
    dir.create(paste0(plot_dir, "/lollipop"), recursive = TRUE, showWarnings = FALSE)
}
if (!dir.exists(paste0(plot_dir, "/cnet"))) {
    dir.create(paste0(plot_dir, "/cnet"), recursive = TRUE, showWarnings = FALSE)
}
if (!dir.exists(paste0(plot_dir, "/heat"))) {
    dir.create(paste0(plot_dir, "/heat"), recursive = TRUE, showWarnings = FALSE)
}
if (!dir.exists(paste0(plot_dir, "/emap"))) {
    dir.create(paste0(plot_dir, "/emap"), recursive = TRUE, showWarnings = FALSE)
}

if (!dir.exists(paste0(plot_dir, "/compareCluster"))) {
    dir.create(paste0(plot_dir, "/compareCluster"), recursive = TRUE, showWarnings = FALSE)
}

# For each result_name, load the corresponding res and res.shrink
result_names <- resultsNames(dds)[!resultsNames(dds) == "Intercept"]

for (res_name in result_names) {
    # Load the results
    res_file_path <- paste0(snakemake@input[["rds_dir"]], "/res.", analysis_name, ".", res_name, ".rds")
    res <- readRDS(res_file_path)
    
    # Get the significant results
    res.signif <- res[res$padj < padj_threshold & !is.na(res$padj) & abs(res$log2FoldChange) >= lfc_threshold, ]
    
    print(res)
    
    # Load the shrinked results
    res_shrink_file_path <- paste0(snakemake@input[["rds_dir"]], "/res.", analysis_name, ".", res_name, ".shrink.rds")
    res.shrink <- readRDS(res_shrink_file_path)

    # Get the significant results
    res.shrink.signif <- res.shrink[res.shrink$padj < padj_threshold & !is.na(res.shrink$padj) & abs(res.shrink$log2FoldChange) >= lfc_threshold, ]
    
    print(res.shrink)

    # Create the directory if doesn't exist
    if (!dir.exists(snakemake@output[["tables_dir"]])) {
        dir.create(snakemake@output[["tables_dir"]], recursive = TRUE, showWarnings = FALSE)
    }
    
    ##### GO Enrichment Analysis #####
    # Run GO enrichment analysis
    GO.enrich.res <- run_GO_enricher(res, 
                                     GO.list=term_to_gene, 
                                     GO.terms=term_to_name,
                                     universe = rownames(counts(dds)),
                                     padj_threshold = padj_threshold, 
                                     lfc_threshold = lfc_threshold, 
                                     p_cutoff = p_cutoff, 
                                     q_cutoff = q_cutoff,
                                     minGSSize = minGSSize)
    # Write the results to a file
    write.csv(GO.enrich.res, 
            paste0(snakemake@output[["tables_dir"]], "/GO_enrichment.", analysis_name, ".", res_name, ".csv"), 
            row.names = TRUE)
    
    # Run GO enrichment analysis on shrinked results
    GO.enrich.res.shrink <- run_GO_enricher(res.shrink, 
                                            GO.list=term_to_gene, 
                                            GO.terms=term_to_name,
                                            universe = rownames(counts(dds)),
                                            padj_threshold = padj_threshold, 
                                            lfc_threshold = lfc_threshold, 
                                            p_cutoff = p_cutoff, 
                                            q_cutoff = q_cutoff,
                                            minGSSize = minGSSize)
    # Write the results to a file
    write.csv(GO.enrich.res.shrink, 
            paste0(snakemake@output[["tables_dir"]], "/GO_enrichment.", analysis_name, ".", res_name, ".shrink.csv"), 
            row.names = TRUE)

    # Generate plots

    ## Dotplot
    print("Generating dotplot...")
    dotplot <- generate_GO_enrich_dotplot(go_results=GO.enrich.res, n_category=n_category, font_size=14, label_format=45, title="") + 
                                        theme(legend.title=element_text(size=18), 
                                              legend.text=element_text(size=16),
                                              axis.text.y = element_text(size = 18),
                                              axis.text.x = element_text(size = 18),
                                              axis.title.x = element_text(size = 18))
    ggsave(paste0(plot_dir, "/dotplot/", "dotplot.", res_name, ".png"), dotplot, width=12, height=14, dpi=350)

    dotplot_shrink <- generate_GO_enrich_dotplot(go_results=GO.enrich.res.shrink, n_category=n_category, font_size=14, label_format=45, title="") + 
                                        theme(legend.title=element_text(size=18), 
                                              legend.text=element_text(size=16),
                                              axis.text.y = element_text(size = 18),
                                              axis.text.x = element_text(size = 18),
                                              axis.title.x = element_text(size = 18))
    ggsave(paste0(plot_dir, "/dotplot/", "dotplot.shrink.", res_name, ".png"), dotplot_shrink, width=12, height=14, dpi=350)

    ## Barplot
    print("Generating barplot...")
    barplot <- generate_GO_enrich_barplot(go_results=GO.enrich.res, n_category=n_category, font_size=14, width=45, title="")  + 
                                        theme(legend.title=element_text(size=18), 
                                              legend.text=element_text(size=16),
                                              axis.text.y = element_text(size = 18),
                                              axis.text.x = element_text(size = 18),
                                              axis.title.x = element_text(size = 18))
    ggsave(paste0(plot_dir, "/barplot/", "barplot.", res_name, ".png"), barplot, width=12, height=14, dpi=350)

    barplot_shrink <- generate_GO_enrich_barplot(go_results=GO.enrich.res.shrink, n_category=n_category, font_size=14, width=45, title="") + 
                                        theme(legend.title=element_text(size=18), 
                                              legend.text=element_text(size=16),
                                              axis.text.y = element_text(size = 18),
                                              axis.text.x = element_text(size = 18),
                                              axis.title.x = element_text(size = 18))
    ggsave(paste0(plot_dir, "/barplot/", "barplot.shrink.", res_name, ".png"), barplot_shrink, width=12, height=14, dpi=350)

    ## Lollipop Plot
    print("Generating lollipop plot...")
    lollipop <- generate_GO_enrich_lollipop(go_results=GO.enrich.res, n_category=n_category, font_size=14, label_format=45, title="") + 
                                        theme(legend.title=element_text(size=18), 
                                              legend.text=element_text(size=16),
                                              axis.text.y = element_text(size = 18),
                                              axis.text.x = element_text(size = 18),
                                              axis.title.x = element_text(size = 18))
    ggsave(paste0(plot_dir, "/lollipop/", "lollipop.", res_name, ".png"), lollipop, width=12, height=14, dpi=350)

    lollipop_shrink <- generate_GO_enrich_lollipop(go_results=GO.enrich.res.shrink, n_category=n_category, font_size=14, label_format=45, title="") + 
                                        theme(legend.title=element_text(size=18), 
                                              legend.text=element_text(size=16),
                                              axis.text.y = element_text(size = 18),
                                              axis.text.x = element_text(size = 18),
                                              axis.title.x = element_text(size = 18))
    ggsave(paste0(plot_dir, "/lollipop/", "lollipop.shrink.", res_name, ".png"), lollipop_shrink, width=12, height=14, dpi=350)

    ## Cnet Plot
    print("Generating cnet plot...")
    cnet_plot <- generate_GO_enrich_cnet(go_results=GO.enrich.res, deseq2_results=res, n_category=n_category, font_size=14, title="")
    ggsave(paste0(plot_dir, "/cnet/", "cnet_plot.", res_name, ".png"), cnet_plot, width=12, height=14, dpi=350)

    cnet_plot_shrink <- generate_GO_enrich_cnet(go_results=GO.enrich.res.shrink, deseq2_results=res.shrink, n_category=n_category, font_size=14, title="")
    ggsave(paste0(plot_dir, "/cnet/", "cnet_plot.shrink.", res_name, ".png"), cnet_plot_shrink, width=12, height=14, dpi=350)

    ## Enrich Heat Plot
    print("Generating heat plot...")
    heat_plot <- generate_GO_enrich_heat(go_results=GO.enrich.res, deseq2_results=res, n_category=n_category, font_size=14, width=45, title="")
    ggsave(paste0(plot_dir, "/heat/", "heat_plot.", res_name, ".png"), heat_plot, width=12, height=14, dpi=350)
    
    heat_plot_shrink <- generate_GO_enrich_heat(go_results=GO.enrich.res.shrink, deseq2_results=res.shrink, n_category=n_category, font_size=14, width=45, title="")
    ggsave(paste0(plot_dir, "/heat/", "heat_plot.shrink.", res_name, ".png"), heat_plot_shrink, width=12, height=14, dpi=350)
    
    ## Emap Plot
    print("Generating emap plot...")
    emap <- generate_GO_enrich_emap(go_results=GO.enrich.res, deseq2_results=res, n_category=n_category, font_size=14, title="")
    ggsave(paste0(plot_dir, "/emap/", "emap.", res_name, ".png"), emap, width=12, height=14, dpi=350)

    emap_shrink <- generate_GO_enrich_emap(go_results=GO.enrich.res.shrink, deseq2_results=res.shrink, n_category=n_category, font_size=14, title="")
    ggsave(paste0(plot_dir, "/emap/", "emap.shrink.", res_name, ".png"), emap_shrink, width=12, height=14, dpi=350)   

    ## compareCluster_upDOWN
    print("Running compareCluster GO enrichment for up/down regulated genes...")
    compareCluster_upDOWN <- run_compareCluster_upDOWN(res=res.signif, dds=dds, 
                                                        term_to_gene=term_to_gene, term_to_name=term_to_name,
                                                        p_cutoff = p_cutoff, q_cutoff = q_cutoff,
                                                        minGSSize = minGSSize)
    
    # if not compareCluster_upDOWN is NULL
    if (!is.null(compareCluster_upDOWN)) {
        print(compareCluster_upDOWN)
        print("Generating compareCluster dotplot for up/down regulated genes...")
        compareCluster_upDOWN.dotplot <- generate_GO_enrich_dotplot(compareCluster_upDOWN, 
                                                                                n_category=n_category, 
                                                                                font_size=14, label_format=45, 
                                                                                title="") + 
                                                                                xlab(NULL) + 
                                                                                theme(legend.title=element_text(size=18), 
                                                                                      legend.text=element_text(size=16),
                                                                                      axis.text.y = element_text(size = 20),
                                                                                      axis.text.x = element_text(size = 20),
                                                                                      axis.title.x = element_text(size = 20),
                                                                                      panel.grid.major.y = element_line(linetype='dotted', color='#808080'),
                                                                                      panel.grid.major.x = element_blank()) + labs(caption = "")
        ggsave(paste0(plot_dir, "/compareCluster/", "compareCluster_upDOWN.dotplot.", res_name, ".png"), compareCluster_upDOWN.dotplot, width=14, height=18, dpi=350)
        # print("Generating compareCluster cnetplot for up/down regulated genes...")
        # compareCluster_upDOWN.cnetplot <- cnetplot(compareCluster_upDOWN, 
        #                                            showCategory=n_category,
        #                                            circular = FALSE,
        #                                            cex.params = list(gene_label = 1.5, category_label = 2.5, gene_node=1.5, category_node=2.5),
        #                                            node_label="category",
        #                                            categorySize="pvalue",
        #                                            layout = "nicely") +
        #                                            labs(title = "",subtitle = "") +
        #                                            theme(plot.title = element_text(size = 16, face = "bold"),
        #                                                  plot.subtitle = element_text(size = 14),
        #                                                  plot.caption = element_text(size = 12, face = "italic"),
        #                                                  legend.title=element_blank(),
        #                                                  legend.position="top",
        #                                                  legend.text=element_text(size=18)) + 
        #                                                              guides(size = guide_legend(override.aes=list(shape=1)))
        # print("Saving compareCluster cnetplot for up/down regulated genes...")
        # ggsave(paste0(plot_dir, "/compareCluster/", "compareCluster_upDOWN.cnetplot.", res_name, ".png"), compareCluster_upDOWN.cnetplot, width=16, height=14, dpi=350)
    
    }  else {
        print("compareCluster_upDOWN is NULL. Skipping...")
    }
    
    print("Running compareCluster GO enrichment for up/down regulated genes... (shrink)")
    compareCluster_upDOWN_shrink <- run_compareCluster_upDOWN(res=res.shrink.signif, dds=dds, 
                                                        term_to_gene=term_to_gene, term_to_name=term_to_name,
                                                        p_cutoff = p_cutoff, q_cutoff = q_cutoff,
                                                        minGSSize = minGSSize) 
    

    # if not compareCluster_upDOWN_shrink is NULL
    if (!is.null(compareCluster_upDOWN_shrink)) {
        print(compareCluster_upDOWN_shrink)
        print("Generating compareCluster dotplot for up/down regulated genes... (shrink)")
        compareCluster_upDOWN_shrink.dotplot <- generate_GO_enrich_dotplot(compareCluster_upDOWN_shrink, 
                                                                                    n_category=n_category, 
                                                                                    font_size=14, label_format=35, 
                                                                                    title="") + 
                                                                                    xlab(NULL) + 
                                                                                    theme(legend.title=element_text(size=20), 
                                                                                          legend.text=element_text(size=18),
                                                                                          axis.text.y = element_text(size = 22),
                                                                                          axis.text.x = element_text(size = 22),
                                                                                          axis.title.x = element_text(size = 22),
                                                                                          panel.grid.major.y = element_line(linetype='dotted', color='#808080'),
                                                                                          panel.grid.major.x = element_blank()) + labs(caption = "")
        print("Saving compareCluster dotplot for up/down regulated genes... (shrink)")
        ggsave(paste0(plot_dir, "/compareCluster/", "compareCluster_upDOWN.dotplot.shrink.", res_name, ".png"), compareCluster_upDOWN_shrink.dotplot, width=14, height=16, dpi=350)
        
        # print("Generating compareCluster cnetplot for up/down regulated genes... (shrink)")
        # compareCluster_upDOWN_shrink.cnetplot <- cnetplot(compareCluster_upDOWN_shrink, 
        #                                               showCategory=n_category,
        #                                               circular = FALSE,
        #                                               cex.params = list(gene_label = 1.5, category_label = 2.5, gene_node=1.5, category_node=2.5),
        #                                               node_label="category",
        #                                               categorySize="pvalue",
        #                                               layout = "nicely") +
        #                                               labs(title = "",subtitle = "") +
        #                                               theme(plot.title = element_text(size = 16, face = "bold"),
        #                                                     plot.subtitle = element_text(size = 14),
        #                                                     plot.caption = element_text(size = 12, face = "italic"),
        #                                                     legend.title=element_blank(),
        #                                                     legend.position="top",
        #                                                     legend.text=element_text(size=18)) + 
        #                                                                 guides(size = guide_legend(override.aes=list(shape=1)))
        # print("Saving compareCluster cnetplot for up/down regulated genes... (shrink)")
        # ggsave(paste0(plot_dir, "/compareCluster/", "compareCluster_upDOWN.cnetplot.shrink.", res_name, ".png"), compareCluster_upDOWN_shrink.cnetplot, width=16, height=14, dpi=350)
    
    }  else {
        print("compareCluster_upDOWN_shrink is NULL. Skipping...")
    }   
    print(paste0("Done: ", res_name))
}

# Close logging
sink(type="message")
sink()
