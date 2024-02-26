# Date: 2023-11-02
# Author: Adnan Kivanc Corut

# This script contains common functions for functional enrichment analysis and clusterProfiler package

# ====== GO Enrichment Analysis Functions ====== #

# GO enrichment analysis function
#"""
#This function runs GO enrichment analysis using enricher function from clusterProfiler package
#Parameters:
#    sig.deseq2.res: DESeq2 results
#    GO.list: GO list (TERM2GENE)
#    GO.terms: GO terms (TERM2NAME)
#    p_cutoff: p-value cutoff
#    q_cutoff: q-value cutoff
#"""
run_GO_enricher <- function(deseq2.res, GO.list, GO.terms, universe, padj_threshold = 0.05, lfc_threshold, p_cutoff=0.05, q_cutoff=0.05, minGSSize = 5, maxGSSize = 5000){

  print(nrow(deseq2.res))
  # Select significant DE genes
  deseq2.signif.res <- deseq2.res[deseq2.res$padj < padj_threshold & !is.na(deseq2.res$padj), ]
  print(nrow(deseq2.signif.res))
  # If lfc_threshold is specified, select DE genes with log2 fold change >= lfc_threshold
  if (!missing(lfc_threshold)) {
    deseq2.signif.res <- deseq2.signif.res[abs(deseq2.signif.res$log2FoldChange) >= lfc_threshold, ]
  }
  print(nrow(deseq2.signif.res))
  print(min(abs(deseq2.signif.res$log2FoldChange)))
  # Get the list of significant DE genes
  signif.de.genes <- rownames(deseq2.signif.res)

  # if universe is missing run this
  if (missing(universe)) {
    # Run GO enrichment analysis
    GO.enrich.res <- enricher(signif.de.genes, 
                              TERM2GENE=GO.list, 
                              TERM2NAME = GO.terms, 
                              pvalueCutoff = p_cutoff,
                              pAdjustMethod = "BH",
                              minGSSize = minGSSize,
                              maxGSSize = maxGSSize,
                              qvalueCutoff = q_cutoff)
    return(GO.enrich.res)
  } else {
    # Run GO enrichment analysis
    GO.enrich.res <- enricher(signif.de.genes, 
                              TERM2GENE=GO.list, 
                              TERM2NAME = GO.terms, 
                              pvalueCutoff = p_cutoff,
                              pAdjustMethod = "BH",
                              minGSSize = minGSSize,
                              maxGSSize = maxGSSize,
                              qvalueCutoff = q_cutoff,
                              universe = universe)
    return(GO.enrich.res)
  }

}

# GO enrichment dotplot function
#"""
#This function generates a dotplot of enriched GO terms
#
#Parameters:
#    go_results: GO enrichment results
#    n_category: number of categories to show
#    font_size: font size
#    label_format: label format
#    title: title of the plot
#"""
generate_GO_enrich_dotplot <- function(go_results, n_category=20, font_size=16, label_format=45, title=""){
  
  # Check if go_results is empty
  if (nrow(go_results) == 0 || is.null(go_results)) {
    print("No results to plot")
    return()
  }

  # Remove the rows with NA in Description
	go_results_nona <- dplyr::filter(go_results, Description != "NA")

  # Generate dotplot
	dot_plot <- clusterProfiler::dotplot(go_results_nona, showCategory=n_category, 
										font.size=font_size, title=title, label_format = label_format) + 
		labs(caption = bquote("* p-value cutoff= 0.05; q-value cutoff= 0.05")) +
		theme(plot.title = element_text(size = 16, face = "bold"),
          plot.caption = element_text(size = 12, face = "italic"),
			  	axis.text.y = element_text(size = 16),
			  	axis.text.x = element_text(size = 16),
				  axis.title.x = element_text(size = 18),
				  panel.grid.major = element_blank(),
        	panel.grid.minor = element_blank()) +
		scale_colour_gradientn(colours = pal_npg()(4))

  	return(dot_plot)
}

# GO enrichment barplot function
#"""
#This function generates a barplot of enriched GO terms
#
#Parameters:
#    go_results: GO enrichment results
#    n_category: number of categories to show
#    font_size: font size
#    width: width of the plot
#    title: title of the plot
#"""
generate_GO_enrich_barplot <- function(go_results, n_category=20, font_size=12, width=45, title="", subtitle=""){
  	
  # Check if go.res is empty
  if (nrow(go_results) == 0 || is.null(go_results)) {
    print("No results to plot")
    return()
  }    

  # Remove the rows with NA in Description
	go_results_nona <- dplyr::filter(go_results, Description != "NA")

	# Generate bar plot
  bar_plot <- barplot(go_results_nona, 
                      drop=TRUE, font.size = font_size,
                      showCategory=n_category, x = "Count") + 
  labs(title = title, 
		subtitle = subtitle, 
		y ="Count",
       	caption = bquote("* p-value cutoff= 0.05; q-value cutoff= 0.05")) +
  theme(plot.title = element_text(size = 16, face = "bold"),
        	plot.subtitle = element_text(size = 14),
        	plot.caption = element_text(size = 12, face = "italic"),
		axis.text.x = element_text(size = 12),
		axis.title.x = element_text(size = 16), 
		axis.text.y = element_text(size = 14),
		axis.title.y = element_text(size = 16)) +
	scale_fill_gradientn(colours = pal_npg()(4))

	return(bar_plot)
}

# GO enrichment lollipop function
#"""
#This function generates a lollipop plot of enriched GO terms
#
#Parameters:
#    go_results: GO enrichment results
#    n_category: number of categories to show
#    font_size: font size
#    label_format: label format
#    title: title of the plot
#"""
generate_GO_enrich_lollipop <- function(go_results, n_category=10, font_size=13, label_format=45, title=""){
	
  # Check if go.res is empty
  if (nrow(go_results) == 0 || is.null(go_results)) {
    print("No results to plot")
    return()
  }   

	# Remove the rows with NA in Description
	go_results_nona <- dplyr::filter(go_results, Description != "NA")
	
	# Add rich factor column
	go_enrich_result_with_rich_factor <- mutate(go_results_nona@result, richFactor = Count / as.numeric (sub("/\\d+", "", BgRatio)))	
	
	go_enrich_result_with_rich_factor <- go_enrich_result_with_rich_factor %>% 
		arrange(pvalue) %>%
  		slice_head(n = n_category)
	
	# Generate lollipop chart
	lollipop_chart <- ggplot(go_enrich_result_with_rich_factor, showCategory = n_category, 
	        aes(richFactor,
	        fct_reorder(Description, richFactor))) + 
	        geom_segment(aes(xend=0, yend = Description)) + 
	        geom_point(aes(color=p.adjust, size = Count)) + 
	        scale_color_gradientn(colours=c("#7e62a3", "#46bac2", "#f7ca64"),
	                    trans = "log10", 
	                    guide=guide_colorbar(reverse=TRUE, order=1)) + 
	        scale_size_continuous(range=c(2, 10)) +
			scale_y_discrete(labels = function(go_enrich_result_with_rich_factor) str_wrap(go_enrich_result_with_rich_factor, width = label_format)) +
	        xlab("Rich Factor") +
	        ylab(NULL) +
	        labs(title=title,
				caption = bquote("* p-value cutoff= 0.05; q-value cutoff= 0.05")) +
	        theme(plot.title = element_text(size = 16, face = "bold"),
	              plot.subtitle = element_text(size = 14),
	              plot.caption = element_text(size = 12, face = "italic"),
				  axis.text.y = element_text(size = 10)) +
	        theme_dose(12)

	return(lollipop_chart)
}

# GO enrichment cnetplot function
#"""
#This function generates a cnetplot of enriched GO terms
#
#Parameters:
#    go_results: GO enrichment results
#    deseq2_results: DESeq2 results
#    n_category: number of categories to show
#    font_size: font size
#    title: title of the plot
#"""
generate_GO_enrich_cnet <- function(go_results, deseq2_results, n_category=10, font_size=13, title="" , subtitle = ""){
	
  # Check if go.res is empty
  if (nrow(go_results) == 0 || is.null(go_results)) {
    print("No results to plot")
    return()
  }   

	# Remove the rows with NA in Description
	go_results_nona <- dplyr::filter(go_results, Description != "NA")

	# Get the log2FoldChange of DE genes
	logchange_col <- deseq2_results[["log2FoldChange"]]
	genes <- rownames(deseq2_results)
	foldChange <- setNames(logchange_col, genes)

	# Generate cnetplot
  cnet_plot <- cnetplot(go_results_nona,
                          color.params = list(foldChange = foldChange, edge=FALSE),
                          showCategory=n_category,
  						circular = FALSE,
                          cex.params = list(gene_label = 0.8, category_label = 1, gene_node=0.8, category_node=1),
                          layout = "nicely") +
                          labs(title = title, 
                             subtitle = subtitle) +
                          theme(plot.title = element_text(size = 16, face = "bold"),
                                  plot.subtitle = element_text(size = 14),
                                  plot.caption = element_text(size = 12, face = "italic"))
	
	return(cnet_plot) 
}

generate_GO_enrich_heat <- function(go_results, deseq2_results, n_category=10, font_size=13, width=45, title="" , subtitle = ""){

  # Check if go.res is empty
  if (nrow(go_results) == 0 || is.null(go_results)) {
    print("No results to plot")
    return()
  }   

  # Remove rows with NA values
  go_results_nona <- dplyr::filter(go_results, Description != "NA")

	# Get log2FoldChange values
	logchange_col <- deseq2_results[["log2FoldChange"]]
	genes <- rownames(deseq2_results)
	foldChange <- setNames(logchange_col, genes) 

  # Generate heat plot
  heat_plot <- heatplot(go_results_nona, 
                            foldChange=foldChange, 
                            showCategory=n_category) +
							scale_y_discrete(labels = function(go_results_nona) str_wrap(go_results_nona, width = width)) +
                            labs(title = title, 
                               subtitle = subtitle, 
                               caption = bquote("* p-value cutoff= 0.05; q-value cutoff= 0.05")) +
                            theme(plot.title = element_text(size = 16, face = "bold"),
                                    plot.subtitle = element_text(size = 14),
                                    plot.caption = element_text(size = 12, face = "italic"),
                                    panel.border = element_rect(color = "black", fill = NA),
									axis.text.x = element_text(size = 14),
									axis.text.y = element_text(size = 14)) + 
                  xlab("")

	return(heat_plot)
}

generate_GO_enrich_emap <- function(go_results, deseq2_results, n_category=10, font_size=13, title="" , subtitle = ""){
	
  # Check if go.res is empty
  if (nrow(go_results) == 0 || is.null(go_results)) {
    print("No results to plot")
    return()
  }   

	# Remove rows with NA values
  go_results_nona <- dplyr::filter(go_results, Description != "NA")
	
	# Calculate pairwise term similarity
  go_results_nona <- pairwise_termsim(go_results_nona)
    
	# Generate emap plot
  emap_plot <- emapplot(go_results_nona, 
                            showCategory = n_category, 
                            cex.params = list(category_node = 1, line=0.8)) +
                            labs(title = title, 
                               subtitle = subtitle) +
                            theme(plot.title = element_text(size = 16, face = "bold"),
                                    plot.subtitle = element_text(size = 14),
                                    plot.caption = element_text(size = 12, face = "italic"))

	return(emap_plot)
}

# ====== KEGG Enrichment Analysis Functions ====== #

# KEGG enrichment analysis function
#"""
#This function runs KEGG enrichment analysis using enrichKEGG function from clusterProfiler package
#Parameters:
#    tcp_to_gene_df: TCP to gene dataframe
#    res: deseq2 results
#    dds: DESeqDataSet
#    enrichment_type: enrichment type (pathway or module)
#    org: organism
#    p_cutoff: p-value cutoff
#    q_cutoff: q-value cutoff
#"""
extract_and_kegg_enrich <- function(tcp_to_gene_df, res, dds, enrichment_type = "pathway", org='', p_cutoff=0.1, q_cutoff=0.2, padj_threshold = 0.05, lfc_threshold = 0.0) {
  
  # Get significant DE transcripts
  signif_res <- res[res$padj < padj_threshold & !is.na(res$padj) & abs(res$log2FoldChange) >= lfc_threshold, ]
  print(signif_res)
  signif_DE_transcripts <- rownames(signif_res)

  print("Number of DE transcripts:")
  print(length(signif_DE_transcripts))

  # Extract the subset of the dataframe where target_id is in the vector list
  subset_df <- tcp_to_gene_df[tcp_to_gene_df$target_id %in% signif_DE_transcripts & tcp_to_gene_df$TAIR_hit != "", ]
  
  print(head(tcp_to_gene_df))

  cat("Subset dataframe:\n")
  print(head(subset_df))
  
  # Remove the ".X" pattern from the TAIR_hit column
  subset_df$TAIR_hit <- sub("\\.[0-9]$", "", subset_df$TAIR_hit)
  
  # Extract the unique TAIR_hit values from this subset
  tair_hits_list <- unique(subset_df$TAIR_hit)
  print(tair_hits_list)
  print("Number of Gene hits:")
  print(length(tair_hits_list))
  
  universe.dds <- rownames(counts(dds))
  universe.df <- tcp_to_gene_df[tcp_to_gene_df$target_id %in% universe.dds & tcp_to_gene_df$TAIR_hit != "", ]
  universe.df$TAIR_hit <- sub("\\.[0-9]$", "", universe.df$TAIR_hit)
  universe.list <- universe.df$TAIR_hit
  
  
  # Decide on enrichment type
  if (enrichment_type == "pathway") {
    enrichment_result <- enrichKEGG(gene         = tair_hits_list,
                                    organism     = org,
                                    universe     = universe.list,
                                    pvalueCutoff = p_cutoff)
  } else if (enrichment_type == "module") {
    enrichment_result <- enrichMKEGG(gene         = tair_hits_list,
                                     organism     = org,
                                     universe     = universe.list,
                                     pvalueCutoff = p_cutoff,
                                     qvalueCutoff = q_cutoff)
  } else {
    stop("Invalid enrichment type. Choose either 'pathway' or 'module'.")
  }
  
  return(enrichment_result)
}

# KEGG enrichment bubble plot function
#"""
#This function generates a bubble plot of enriched KEGG pathways
#
#Parameters:
#    kegg_results: KEGG enrichment results
#    n_category: number of categories to show
#"""
generate_kegg_enrich_bubble_plot <- function(kegg_results, n_category = 20) {

    # Check if kegg_results is NULL
    if (is.null(kegg_results)) {
      print("No results to plot")
      return()
    }

	  # Remove the rows with NA in Description
	  kegg_results_nona <- dplyr::filter(kegg_results, Description != "NA")
	  
	  # Remove - Arabidopsis thaliana (thale cress) from Description
	  kegg_results_nona@result$Description <- sub(" - Arabidopsis thaliana \\(thale cress\\)", "", kegg_results_nona@result$Description)
	  
	  # Add rich factor column
	  kegg_results_with_rich_factor <- mutate(kegg_results_nona@result, richFactor = Count / as.numeric (sub("/\\d+", "", BgRatio)))	
	  
	  kegg_results_with_rich_factor <- kegg_results_with_rich_factor %>% 
	  	arrange(pvalue) %>%
    		slice_head(n = n_category)
      
    # Create the bubble plot
    p <- ggplot(kegg_results_with_rich_factor, aes(x = richFactor, y = Description, size = Count * 4, color = p.adjust)) +
        geom_point(alpha=0.7) +
        scale_size_continuous(range = c(3, 10)) +
        scale_color_gradientn(colours=c("#7e62a3", "#46bac2", "#f7ca64"),
	          trans = "log10", 
	          guide=guide_colorbar(reverse=TRUE, order=1)) + 
        theme_bw() +
        labs(title = "",
             x = "Rich Factor",
             y = "Pathway",
             size = "Gene Count",
             color = "Adjusted p-value") +
        theme(legend.position="right",
              legend.title=element_text(size=18), 
              legend.text=element_text(size=16),
              axis.text.y = element_text(size = 16),
              axis.text.x = element_text(size = 16),
              axis.title.x = element_text(size = 18), 
		          axis.title.y = element_text(size = 18),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
    
    return(p)
}


# A compareCluster function for up/down regulated genes
#"""
#This function runs compareCluster function from clusterProfiler package for up/down regulated genes
#
#Parameters:
#    res: DESeq2 results
#    dds: DESeqDataSet
#    term_to_gene: term to gene list
#    term_to_name: term to name list
#    p_cutoff: p-value cutoff
#    q_cutoff: q-value cutoff
#    minGSSize: minimum gene set size
#    maxGSSize: maximum gene set size
#"""
run_compareCluster_upDOWN <- function(res, dds, term_to_gene, term_to_name, p_cutoff=0.05, q_cutoff=0.1, minGSSize = 5, maxGSSize = 5000){
    df <- data.frame(Gene=rownames(res), FC=res$log2FoldChange)
    df$group <- "UP-regulated"
    df$group[df$FC < 0] <- "Down-regulated"
    formula_res <- compareCluster(Gene~group, data=df, fun="enricher",
                                  TERM2GENE=term_to_gene,
                                  TERM2NAME=term_to_name,
                                  universe=rownames(counts(dds)),
                                  pvalueCutoff = p_cutoff,
                                  pAdjustMethod = "BH",
                                  minGSSize = minGSSize,
                                  maxGSSize = maxGSSize,
                                  qvalueCutoff = q_cutoff)

    # if no results, return empty dataframe
    if (is.null(formula_res)) {
      return(NULL)
    }
    
    formula_res <- dplyr::filter(formula_res, Description != "NA")
    return(formula_res)
}