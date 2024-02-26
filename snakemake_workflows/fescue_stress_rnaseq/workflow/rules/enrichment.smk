# Define the path for common enrichemt functions script
ENRICH_FUNCS_SCRIPT= os.path.abspath("./scripts/clusterprofiler_common_functions.R")

# ==============================================================================================
#     Importing GO Terms
# ==============================================================================================

rule import_go:
    input:
        go_info= "resources/ref/annotation/fescue.go_list.txt",
    output:
        term_to_gene= WORKDIR + "/results/go_enrichment/mappings/term_to_gene.rds",
        term_to_name= WORKDIR + "/results/go_enrichment/mappings/term_to_name.rds"
    conda:
        "../envs/enrich.yaml"
    message: "Importing GO terms..."
    log:
        WORKDIR + "/logs/go_enrichment/import_go.log"
    script:
        "../scripts/import_go.R"

# ==============================================================================================
#     Running GO Enrichment Analysis
# ==============================================================================================

rule run_go_enrichment:
    input:
        dds_subset = WORKDIR + "/results/deseq2/subset/dds_{analysis_name}.rds",
        rds_dir = WORKDIR + "/results/deseq2/diffexp/comparison_{analysis_name}",
        term_to_gene= WORKDIR + "/results/go_enrichment/mappings/term_to_gene.rds",
        term_to_name= WORKDIR + "/results/go_enrichment/mappings/term_to_name.rds",
        tcp_to_gene = "resources/ref/annotation/fescue.transcriptome_info.txt",
    output:
        tables_dir = report(
                        directory(WORKDIR + "/results/tables/enrichment/comparison_{analysis_name}/go"),
                        patterns=["{name}.csv"],
                        caption="../report/clusterProfiler-go_enrichment.rst",
                        category="Functional Enrichment-GO",
                        subcategory="{analysis_name}"),
        plots_dir = report(
                        directory(WORKDIR + "/results/plots/enrichment/comparison_{analysis_name}/go"),
                        patterns=["{name}.png"],
                        caption="../report/clusterProfiler-go_enrichment.rst",
                        category="Functional Enrichment-GO",
                        subcategory="{analysis_name}"),
    params:
        funcs_script = ENRICH_FUNCS_SCRIPT,
        padj_threshold = config["params"]["deseq2"]["padj_threshold"],
        lfc_threshold = config["params"]["deseq2"]["lfc_threshold"],
        p_cutoff = config["params"]["clusterProfiler"]["p_cutoff"],
        q_cutoff = config["params"]["clusterProfiler"]["q_cutoff"],
        minGSSize = config["params"]["clusterProfiler"]["minGSSize"],
        p_adj_method = config["params"]["clusterProfiler"]["p_adj_method"],
        n_category = config["params"]["clusterProfiler"]["n_category"],
    conda:
        "../envs/deseq2_diffexp.yaml"
    log:
        WORKDIR + "/logs/enrichment/go/run_go_enrichment_{analysis_name}.log"
    threads: 
        config["params"]["clusterProfiler"]["threads"]
    script:
        "../scripts/clusterprofiler_go_enrich.R"

# ==============================================================================================
#     Running KEGG Enrichment Analysis
# ==============================================================================================

rule run_kegg_enrichment:
    input:
        dds_subset = WORKDIR + "/results/deseq2/subset/dds_{analysis_name}.rds",
        rds_dir = WORKDIR + "/results/deseq2/diffexp/comparison_{analysis_name}",
        tcp_to_gene = "resources/ref/annotation/fescue.transcriptome_info.txt",
    output:
        tables_dir = report(
                        directory(WORKDIR + "/results/tables/enrichment/comparison_{analysis_name}/kegg"),
                        patterns=["{name}.csv"],
                        caption="../report/clusterProfiler-kegg_enrichment.rst",
                        category="Functional Enrichment-KEGG",
                        subcategory="{analysis_name}"),
        plots_dir = report(
                        directory(WORKDIR + "/results/plots/enrichment/comparison_{analysis_name}/kegg"),
                        patterns=["{name}.png"],
                        caption="../report/clusterProfiler-kegg_enrichment.rst",
                        category="Functional Enrichment-KEGG",
                        subcategory="{analysis_name}"),
    params:
        funcs_script = ENRICH_FUNCS_SCRIPT,
        padj_threshold = config["params"]["deseq2"]["lfc_threshold"],
        lfc_threshold = config["params"]["deseq2"]["lfc_threshold"],
        p_cutoff = config["params"]["clusterProfiler"]["p_cutoff"],
        q_cutoff = config["params"]["clusterProfiler"]["q_cutoff"],
        kegg_organism = config["params"]["clusterProfiler"]["kegg_organism"],
    conda:
        "../envs/kegg_enrich.yaml"
    log:
        WORKDIR + "/logs/enrichment/kegg/run_kegg_enrichment_{analysis_name}.log"
    threads: 
        config["params"]["clusterProfiler"]["threads"]
    script:
        "../scripts/clusterprofiler_kegg_enrich.R"

# ==============================================================================================
#     Running WikiPathways Enrichment Analysis
# ==============================================================================================

rule run_wikipathways_enrichment:
    input:
        dds_subset = WORKDIR + "/results/deseq2/subset/dds_{analysis_name}.rds",
        rds_dir = WORKDIR + "/results/deseq2/diffexp/comparison_{analysis_name}",
        tcp_to_gene = "resources/ref/annotation/fescue.transcriptome_info.txt",
    output:
        tables_dir = report(
                        directory(WORKDIR + "/results/tables/enrichment/comparison_{analysis_name}/wikipathways"),
                        patterns=["{name}.csv"],
                        caption="../report/clusterProfiler-wikipathways_enrichment.rst",
                        category="Functional Enrichment-WikiPathways",
                        subcategory="{analysis_name}"),
        plots_dir = report(
                        directory(WORKDIR + "/results/plots/enrichment/comparison_{analysis_name}/wikipathways"),
                        patterns=["{name}.png"],
                        caption="../report/clusterProfiler-wikipathways_enrichment.rst",
                        category="Functional Enrichment-WikiPathways",
                        subcategory="{analysis_name}"),
    params:
        funcs_script = ENRICH_FUNCS_SCRIPT,
        padj_threshold = config["params"]["deseq2"]["padj_threshold"],
        lfc_threshold = config["params"]["deseq2"]["lfc_threshold"],
        p_cutoff = config["params"]["clusterProfiler"]["p_cutoff"],
        q_cutoff = config["params"]["clusterProfiler"]["q_cutoff"],
        n_category = config["params"]["clusterProfiler"]["n_category"],
    conda:
        "../envs/wikipathways.yaml"
    log:
        WORKDIR + "/logs/enrichment/wikipathways/run_wikipathways_enrichment_{analysis_name}.log"
    threads: 
        config["params"]["clusterProfiler"]["threads"]
    script:
        "../scripts/clusterprofiler_wikipathways.R"

# ==============================================================================================
#     Running compareCluster
# ==============================================================================================

rule run_compareCluster:
    input:
        dds_subset = WORKDIR + "/results/deseq2/subset/dds_{analysis_name}.rds",
        rds_dir = WORKDIR + "/results/deseq2/diffexp/comparison_{analysis_name}",
        term_to_gene= WORKDIR + "/results/go_enrichment/mappings/term_to_gene.rds",
        term_to_name= WORKDIR + "/results/go_enrichment/mappings/term_to_name.rds",
    output:
        plots_dir = report(
                        directory(WORKDIR + "/results/plots/enrichment/comparison_{analysis_name}/compareCluster"),
                        patterns=["{name}.png"],
                        caption="../report/clusterProfiler-compareCluster.rst",
                        category="Functional Enrichment-compareCluster",
                        subcategory="{analysis_name}"),
        lists_dir = directory(WORKDIR + "/results/go_enrichment/compareCluster/comparison_{analysis_name}"),
    params:
        analysis=lambda wildcards: next(a for a in config["params"]["deseq2"]["analyses"] if a["name"] == wildcards.analysis_name),
        padj_threshold = config["params"]["deseq2"]["padj_threshold"],
        lfc_threshold = config["params"]["deseq2"]["lfc_threshold"],
        p_cutoff = config["params"]["clusterProfiler"]["p_cutoff"],
        q_cutoff = config["params"]["clusterProfiler"]["q_cutoff"],
        minGSSize = config["params"]["clusterProfiler"]["minGSSize"],
        p_adj_method = config["params"]["clusterProfiler"]["p_adj_method"],
        n_category = config["params"]["clusterProfiler"]["n_category"],
    conda:
        "../envs/deseq2_diffexp.yaml"
    log:
        WORKDIR + "/logs/enrichment/compareCluster/run_compareCluster_{analysis_name}.log"
    threads:
        config["params"]["clusterProfiler"]["threads"]
    script:
        "../scripts/clusterprofiler_compare_cluster.R"

# =========== # =========== # =========== # =========== # =========== # =========== # ========== #

rule run_compareCluster_all:
    input:
        dds = WORKDIR + "/results/deseq2/dds.rds",
        lists_dir = expand(
                        WORKDIR + "/results/go_enrichment/compareCluster/comparison_{analysis_name}",
                        analysis_name = [a["name"] for a in config["params"]["deseq2"]["analyses"]]
                    ),
        term_to_gene= WORKDIR + "/results/go_enrichment/mappings/term_to_gene.rds",
        term_to_name= WORKDIR + "/results/go_enrichment/mappings/term_to_name.rds",
    output:
        plots_dir = directory(WORKDIR + "/results/plots/enrichment/compareCluster_all"),
        done = touch(WORKDIR + "/results/plots/enrichment/compareCluster_all/compareCluster_all.done")
    params:
        padj_threshold = config["params"]["deseq2"]["padj_threshold"],
        lfc_threshold = config["params"]["deseq2"]["lfc_threshold"],
        p_cutoff = config["params"]["clusterProfiler"]["p_cutoff"],
        q_cutoff = config["params"]["clusterProfiler"]["q_cutoff"],
        minGSSize = config["params"]["clusterProfiler"]["minGSSize"],
        p_adj_method = config["params"]["clusterProfiler"]["p_adj_method"],
        n_category = config["params"]["clusterProfiler"]["n_category"],
        use_shrink = True
    conda:
        "../envs/deseq2_diffexp.yaml"
    log:
        WORKDIR + "/logs/enrichment/compareCluster_all/run_compareCluster_all.log"
    threads:
        config["params"]["clusterProfiler"]["threads"]
    script:
        "../scripts/clusterprofiler_compare_cluster_all.R"
    
# ==============================================================================================
#     Running Go Enrichment Analysis on Unique DEGs
# ==============================================================================================

rule run_go_enrichment_unique_degs:
    input:
        dds = WORKDIR + "/results/deseq2/dds.rds",
        lists_dir = expand(
                        WORKDIR + "/results/go_enrichment/compareCluster/comparison_{analysis_name}",
                        analysis_name = [a["name"] for a in config["params"]["deseq2"]["analyses"]]
                    ),
        term_to_gene= WORKDIR + "/results/go_enrichment/mappings/term_to_gene.rds",
        term_to_name= WORKDIR + "/results/go_enrichment/mappings/term_to_name.rds",
    output:
        plots_dir = directory(WORKDIR + "/results/plots/enrichment/go_enrichment_unique_degs"),
        done = touch(WORKDIR + "/results/plots/enrichment/go_enrichment_unique_degs/go_enrichment_unique_degs.done")
    params:
        # comparisons_of_interest = ("comparison_heat_inf_vs_uninf", "comparison_normal_inf_vs_uninf"),
        comparisons_of_interest = ("comparison_infected_heat_vs_normal", "comparison_uninfected_heat_vs_normal"),
        padj_threshold = config["params"]["deseq2"]["padj_threshold"],
        lfc_threshold = config["params"]["deseq2"]["lfc_threshold"],
        p_cutoff = config["params"]["clusterProfiler"]["p_cutoff"],
        q_cutoff = config["params"]["clusterProfiler"]["q_cutoff"],
        minGSSize = config["params"]["clusterProfiler"]["minGSSize"],
        p_adj_method = config["params"]["clusterProfiler"]["p_adj_method"],
        n_category = config["params"]["clusterProfiler"]["n_category"],
        funcs_script = ENRICH_FUNCS_SCRIPT,
    conda:
        "../envs/deseq2_diffexp.yaml"
    log:
        WORKDIR + "/logs/enrichment/go_enrichment_unique_degs/run_go_enrichment_unique_degs.log"
    threads:
        config["params"]["clusterProfiler"]["threads"]
    script:
        "../scripts/clusterprofiler_go_enrich_unique_degs.R"

# ==============================================================================================

rule compare_heat_stress_effect:
    input:
        dds = WORKDIR + "/results/deseq2/dds.rds",
        dds_subset = expand(
                        WORKDIR + "/results/deseq2/subset/dds_{analysis_name}.rds",
                        analysis_name = [a["name"] for a in config["params"]["deseq2"]["analyses"]]
                    ),
        rds_dir = expand(
                        WORKDIR + "/results/deseq2/diffexp/comparison_{analysis_name}",
                        analysis_name = [a["name"] for a in config["params"]["deseq2"]["analyses"]]
                    ),
        term_to_gene= WORKDIR + "/results/go_enrichment/mappings/term_to_gene.rds",
        term_to_name= WORKDIR + "/results/go_enrichment/mappings/term_to_name.rds",
    output:
        plots_dir = directory(WORKDIR + "/results/plots/enrichment/compare_heat_stress_effect"),
        res_dir = directory(WORKDIR + "/results/go_enrichment/compare_heat_stress_effect"),
        done = touch(WORKDIR + "/results/plots/enrichment/compare_heat_stress_effect/compare_heat_stress_effect.done")
    params:
        funcs_script = ENRICH_FUNCS_SCRIPT,
        padj_threshold = config["params"]["deseq2"]["padj_threshold"],
        lfc_threshold = config["params"]["deseq2"]["lfc_threshold"],
        p_cutoff = config["params"]["clusterProfiler"]["p_cutoff"],
        q_cutoff = config["params"]["clusterProfiler"]["q_cutoff"],
        minGSSize = config["params"]["clusterProfiler"]["minGSSize"],
        p_adj_method = config["params"]["clusterProfiler"]["p_adj_method"],
        n_category = config["params"]["clusterProfiler"]["n_category"],
    conda:
        "../envs/venn_enrich.yaml"
    log:
        WORKDIR + "/logs/enrichment/compare_heat_stress_effect/compare_heat_stress_effect.log"
    threads:
        config["params"]["clusterProfiler"]["threads"]
    script:
        "../scripts/compare_heat_stress_effect.R"
        
