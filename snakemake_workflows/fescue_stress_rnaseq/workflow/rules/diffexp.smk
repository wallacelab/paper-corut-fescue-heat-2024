# Define the path for common DE analysis functions script
DEA_FUNCS_SCRIPT= os.path.abspath("./scripts/deseq2_common_functions.R")

# ==================================================================================================================================================
#   Deseq2
# ==================================================================================================================================================

rule deseq2_init:
    input:
        txi = WORKDIR + "/results/tximport/txi.RDS",
        sample_metadata = WORKDIR + "/results/tables/sample_metadata.tsv",
    output:
        WORKDIR + "/results/deseq2/dds.rds",
    params:
        exclude= config["params"]["deseq2"].get("exclude", None),
    conda:
        "../envs/deseq2.yaml"
    threads:
        config["params"]["deseq2"]["threads"],
    log:
        WORKDIR + "/logs/deseq2/deseq2_init.log",
    script:
        "../scripts/deseq2_init.R"

# ==================================================================================================================================================
#   DEseQ2 Plot Dispersion Estimates
# ==================================================================================================================================================

rule deseq2_plot_dispersion:
    input:
        WORKDIR + "/results/deseq2/dds.rds",
    output:
        report(WORKDIR + "/results/plots/deseq2/dispersion_est/dispersion_estimates.png", 
                caption= "../report/deseq2-dispersion_est.rst",
                category= "DEseQ2",
                subcategory="Dispersion Estimates"),
    conda:
        "../envs/deseq2.yaml"
    log:
        WORKDIR + "/logs/deseq2/dispersion_est/dispersion_estimates.log",
    script:
        "../scripts/deseq2_plot_dispersion.R"

# ==================================================================================================================================================
#   Deseq2 Plot PCA
# ==================================================================================================================================================

rule deseq2_pca:
    input:
        WORKDIR + "/results/deseq2/dds.rds",
    output:
        report(WORKDIR + "/results/plots/deseq2/pca/pca.{variable}.png", 
                caption= "../report/deseq2-pca.rst",
                category= "DEseQ2",
                subcategory= "PCA"),
        report(WORKDIR + "/results/plots/deseq2/pca/pca.{variable}.combined.png", 
                caption= "../report/deseq2-pca.rst",
                category= "DEseQ2",
                subcategory= "PCA"),
    conda:
        "../envs/deseq2.yaml"
    log:
        WORKDIR + "/logs/deseq2/pca/pca.{variable}.log",
    script:
        "../scripts/deseq2_plot_pca.R"

# ==================================================================================================================================================
#   Deseq2 Subset Data for Analysis
# ==================================================================================================================================================

rule deseq2_subset:
    input:
        dds = WORKDIR + "/results/deseq2/dds.rds",
    output:
        vsd_subset = WORKDIR + "/results/deseq2/subset/vsd_{analysis_name}.rds",
        dds_subset = WORKDIR + "/results/deseq2/subset/dds_{analysis_name}.rds",
    params:
        analysis=lambda wildcards: next(a for a in config["params"]["deseq2"]["analyses"] if a["name"] == wildcards.analysis_name),
    conda:
        "../envs/deseq2.yaml"
    threads:
        config["params"]["deseq2"]["threads"],
    log:
        WORKDIR + "/logs/deseq2/subset/deseq2_analysis_{analysis_name}.log",
    script:
        "../scripts/deseq2_subset.R"

# ==================================================================================================================================================
#   Deseq2 Exploratory Plots
# ==================================================================================================================================================

rule deseq2_exploratory_plots:
    input:
        vsd_subset = WORKDIR + "/results/deseq2/subset/vsd_{analysis_name}.rds",
        dds_subset = WORKDIR + "/results/deseq2/subset/dds_{analysis_name}.rds",
    output:
        pca_combined = report(WORKDIR + "/results/plots/deseq2/comparison_{analysis_name}/eda/pca.{analysis_name}.combined.png", 
                        caption="../report/deseq2-exploratory.rst",
                        category= "DEseQ2-Exploratory",
                        subcategory= "{analysis_name}"),
        pca_pc1_pc2 = WORKDIR + "/results/plots/deseq2/comparison_{analysis_name}/eda/pca.{analysis_name}.pc1_pc2.png",
        eigencor = report(WORKDIR + "/results/plots/deseq2/comparison_{analysis_name}/eda/{analysis_name}.eigencorplot.png", 
                        caption="../report/deseq2-exploratory.rst",
                        category= "DEseQ2-Exploratory",
                        subcategory= "{analysis_name}"),
        cooks = report(WORKDIR + "/results/plots/deseq2/comparison_{analysis_name}/eda/{analysis_name}.cooks_boxplot.png", 
                        caption="../report/deseq2-exploratory.rst",
                        category= "DEseQ2-Exploratory",
                        subcategory= "{analysis_name}"),
        disp_est = report(WORKDIR + "/results/plots/deseq2/comparison_{analysis_name}/eda/{analysis_name}.disp_est.png", 
                        caption="../report/deseq2-exploratory.rst",
                        category= "DEseQ2-Exploratory",
                        subcategory= "{analysis_name}"),
        density = report(WORKDIR + "/results/plots/deseq2/comparison_{analysis_name}/eda/{analysis_name}.density_plot.png", 
                        caption="../report/deseq2-exploratory.rst",
                        category= "DEseQ2-Exploratory",
                        subcategory= "{analysis_name}"),
        dist_heatmap = report(WORKDIR + "/results/plots/deseq2/comparison_{analysis_name}/eda/{analysis_name}.sample_dist_heatmap.png", 
                        caption="../report/deseq2-exploratory.rst",
                        category= "DEseQ2-Exploratory",
                        subcategory= "{analysis_name}"),
        count_heatmap = report(WORKDIR + "/results/plots/deseq2/comparison_{analysis_name}/eda/{analysis_name}.count_matrix_heatmap.png", 
                        caption="../report/deseq2-exploratory.rst",
                        category= "DEseQ2-Exploratory",
                        subcategory= "{analysis_name}"),
        gene_cluster_heatmap = report(WORKDIR + "/results/plots/deseq2/comparison_{analysis_name}/eda/{analysis_name}.gene_cluster_heatmap.png", 
                        caption="../report/deseq2-exploratory.rst",
                        category= "DEseQ2-Exploratory",
                        subcategory= "{analysis_name}"),
        done = touch(WORKDIR + "/results/plots/deseq2/comparison_{analysis_name}/eda/deseq2_exploratory_plots.done"),
    conda:
        "../envs/deseq2_exploratory.yaml"
    params:
        funcs_script = DEA_FUNCS_SCRIPT,
        color_by = lambda wildcards: [analysis for analysis in config["params"]["deseq2"]["analyses"] if analysis["name"] == wildcards.analysis_name][0]["exploratory_plots"]["color_by"],
        shape_by = lambda wildcards: [analysis for analysis in config["params"]["deseq2"]["analyses"] if analysis["name"] == wildcards.analysis_name][0]["exploratory_plots"]["shape_by"],
        n_top = 500, # number of top genes to plot in PCA
    log:
        WORKDIR + "/logs/deseq2/exploratory_plots/deseq2_exploratory_plots_{analysis_name}.log",
    threads:
        config["params"]["deseq2"]["threads"],
    script:
        "../scripts/deseq2_exploratory_plots.R"

# ==================================================================================================================================================
#   Deseq2 Run Differential Expression
# ==================================================================================================================================================

rule deseq2_diffexp:
    input:
        vsd_subset = WORKDIR + "/results/deseq2/subset/vsd_{analysis_name}.rds",
        dds_subset = WORKDIR + "/results/deseq2/subset/dds_{analysis_name}.rds",
    output:
        tables_dir = directory(WORKDIR + "/results/tables/deseq2/comparison_{analysis_name}/diffexp"),
        rds_dir = directory(WORKDIR + "/results/deseq2/diffexp/comparison_{analysis_name}"),
    params:
        padj_threshold = config["params"]["deseq2"]["padj_threshold"],
        lfc_threshold = config["params"]["deseq2"]["lfc_threshold"],
        p_cutoff = config["params"]["clusterProfiler"]["p_cutoff"],
        q_cutoff = config["params"]["clusterProfiler"]["q_cutoff"],
        minGSSize = config["params"]["clusterProfiler"]["minGSSize"],
        p_adj_method = config["params"]["clusterProfiler"]["p_adj_method"],
        funcs_script = DEA_FUNCS_SCRIPT,
    log:
        WORKDIR + "/logs/deseq2/diffexp/deseq2_diffexp_{analysis_name}.log",
    conda:
        "../envs/deseq2_diffexp.yaml"
    threads:
        config["params"]["deseq2"]["threads"],
    script:
        "../scripts/deseq2_run_diffexp.R"

# ==================================================================================================================================================
# Deseq2 Plot P-value Histogram
# ==================================================================================================================================================

rule deseq2_plot_pval_hist:
    input:
        dds_subset = WORKDIR + "/results/deseq2/subset/dds_{analysis_name}.rds",
        rds_dir = WORKDIR + "/results/deseq2/diffexp/comparison_{analysis_name}",
    output:
        plots_dir = report(
                        directory(WORKDIR + "/results/plots/deseq2/comparison_{analysis_name}/diffexp/pval_hist"),
                        patterns=["{name}.png"],
                        caption="../report/deseq2-pval_hist.rst",
                        category="DEseQ2-Differential Expression",
                        subcategory="{analysis_name}"),
        done = touch(WORKDIR + "/results/plots/deseq2/comparison_{analysis_name}/diffexp/pval_hist/deseq2_plot_pval_hist.done"),
    params:
        analysis=lambda wildcards: next(a for a in config["params"]["deseq2"]["analyses"] if a["name"] == wildcards.analysis_name),
        padj_threshold = config["params"]["deseq2"]["padj_threshold"],
        lfc_threshold = config["params"]["deseq2"]["lfc_threshold"],
        funcs_script = DEA_FUNCS_SCRIPT,
    log:
        WORKDIR + "/logs/deseq2/diffexp/plot_pval_hist/deseq2_plot_pval_hist_{analysis_name}.log",
    conda:
        "../envs/deseq2_diffexp.yaml"
    threads:
        config["params"]["deseq2"]["threads"],
    script:
        "../scripts/deseq2_plot_pval_hist.R"

# ==================================================================================================================================================
#   Deseq2 Plot MA
# ==================================================================================================================================================

rule deseq2_plot_ma:
    input:
        dds_subset = WORKDIR + "/results/deseq2/subset/dds_{analysis_name}.rds",
        rds_dir = WORKDIR + "/results/deseq2/diffexp/comparison_{analysis_name}",
    output:
        plots_dir = report(
                        directory(WORKDIR + "/results/plots/deseq2/comparison_{analysis_name}/diffexp/ma"),
                        patterns=["{name}.png"],
                        caption="../report/deseq2-ma.rst",
                        category="DEseQ2-Differential Expression",
                        subcategory="{analysis_name}"),
        done = touch(WORKDIR + "/results/plots/deseq2/comparison_{analysis_name}/diffexp/ma/deseq2_plot_ma.done"),
    params:
        analysis=lambda wildcards: next(a for a in config["params"]["deseq2"]["analyses"] if a["name"] == wildcards.analysis_name),
        padj_threshold = config["params"]["deseq2"]["padj_threshold"],
        lfc_threshold = config["params"]["deseq2"]["lfc_threshold"],
        funcs_script = DEA_FUNCS_SCRIPT,
    log:
        WORKDIR + "/logs/deseq2/diffexp/plot_ma/deseq2_plot_ma_{analysis_name}.log",
    conda:
        "../envs/deseq2_diffexp.yaml"
    threads:
        config["params"]["deseq2"]["threads"],
    script:
        "../scripts/deseq2_plot_ma.R"

# ==================================================================================================================================================
#   Deseq2 Plot Volcano
# ==================================================================================================================================================

rule deseq2_plot_volcano:
    input:
        dds_subset = WORKDIR + "/results/deseq2/subset/dds_{analysis_name}.rds",
        rds_dir = WORKDIR + "/results/deseq2/diffexp/comparison_{analysis_name}",
    output:
        plots_dir = report(
                        directory(WORKDIR + "/results/plots/deseq2/comparison_{analysis_name}/diffexp/volcano"),
                        patterns=["{name}.png"],
                        caption="../report/deseq2-volcano.rst",
                        category="DEseQ2-Differential Expression",
                        subcategory="{analysis_name}"),
        done = touch(WORKDIR + "/results/plots/deseq2/comparison_{analysis_name}/diffexp/volcano/deseq2_plot_volcano.done"),
    params:
        analysis=lambda wildcards: next(a for a in config["params"]["deseq2"]["analyses"] if a["name"] == wildcards.analysis_name),
        padj_threshold = config["params"]["deseq2"]["padj_threshold"],
        lfc_threshold = config["params"]["deseq2"]["lfc_threshold"],
        funcs_script = DEA_FUNCS_SCRIPT,
    log:
        WORKDIR + "/logs/deseq2/diffexp/plot_volcano/deseq2_plot_volcano_{analysis_name}.log",
    conda:
        "../envs/deseq2_diffexp.yaml"
    threads:
        config["params"]["deseq2"]["threads"],
    script:
        "../scripts/deseq2_plot_volcano.R"

# ==================================================================================================================================================
#   Deseq2 Plot Top LFC Barplot
# ==================================================================================================================================================

rule deseq2_plot_lfc_bar:
    input:
        dds_subset = WORKDIR + "/results/deseq2/subset/dds_{analysis_name}.rds",
        rds_dir = WORKDIR + "/results/deseq2/diffexp/comparison_{analysis_name}",
    output:
        plots_dir = report(
                        directory(WORKDIR + "/results/plots/deseq2/comparison_{analysis_name}/diffexp/lfc_bar"),
                        patterns=["{name}.png"],
                        caption="../report/deseq2-lfc_bar.rst",
                        category="DEseQ2-Differential Expression",
                        subcategory="{analysis_name}"),
        done = touch(WORKDIR + "/results/plots/deseq2/comparison_{analysis_name}/diffexp/lfc_bar/deseq2_plot_lfc_bar.done"),
    params:
        analysis=lambda wildcards: next(a for a in config["params"]["deseq2"]["analyses"] if a["name"] == wildcards.analysis_name),
        padj_threshold = config["params"]["deseq2"]["padj_threshold"],
        lfc_threshold = config["params"]["deseq2"]["lfc_threshold"],
        top_n = 20,
        funcs_script = DEA_FUNCS_SCRIPT,
    log:
        WORKDIR + "/logs/deseq2/diffexp/plot_lfc_bar/deseq2_plot_lfc_bar_{analysis_name}.log",
    conda:
        "../envs/deseq2_diffexp.yaml"
    threads:
        config["params"]["deseq2"]["threads"],
    script:
        "../scripts/deseq2_plot_lfc_bar.R"

# ==================================================================================================================================================
#   Deseq2 Plot Expression Heatmap of Top DEGs
# ==================================================================================================================================================

rule deseq2_plot_exp_heatmap:
    input:
        dds_subset = WORKDIR + "/results/deseq2/subset/dds_{analysis_name}.rds",
        vsd_subset = WORKDIR + "/results/deseq2/subset/vsd_{analysis_name}.rds",
        rds_dir = WORKDIR + "/results/deseq2/diffexp/comparison_{analysis_name}",
    output:
        plots_dir = report(
                        directory(WORKDIR + "/results/plots/deseq2/comparison_{analysis_name}/diffexp/exp_heatmap"),
                        patterns=["{name}.png"],
                        caption="../report/deseq2-exp_heatmap.rst",
                        category="DEseQ2-Differential Expression",
                        subcategory="{analysis_name}"),
        done = touch(WORKDIR + "/results/plots/deseq2/comparison_{analysis_name}/diffexp/exp_heatmap/deseq2_plot_exp_heatmap.done"),
    params:
        analysis=lambda wildcards: next(a for a in config["params"]["deseq2"]["analyses"] if a["name"] == wildcards.analysis_name),
        padj_threshold = config["params"]["deseq2"]["padj_threshold"],
        lfc_threshold = config["params"]["deseq2"]["lfc_threshold"],
        top_n = 20,
        funcs_script = DEA_FUNCS_SCRIPT,
    log:
        WORKDIR + "/logs/deseq2/diffexp/plot_exp_heatmap/deseq2_plot_exp_heatmap_{analysis_name}.log",
    conda:
        "../envs/deseq2_diffexp.yaml"
    threads:
        config["params"]["deseq2"]["threads"],
    script:
        "../scripts/deseq2_plot_exp_heatmap.R"

# ==================================================================================================================================================
#   Deseq2 Plot Top DEG Barplot
# ==================================================================================================================================================

rule deseq2_plot_deg_bar:
    input:
        dds_subset = WORKDIR + "/results/deseq2/subset/dds_{analysis_name}.rds",
        rds_dir = WORKDIR + "/results/deseq2/diffexp/comparison_{analysis_name}",
    output:
        plots_dir = report(
                        directory(WORKDIR + "/results/plots/deseq2/comparison_{analysis_name}/diffexp/deg_bar"),
                        patterns=["{name}.png"],
                        caption="../report/deseq2-deg_bar.rst",
                        category="DEseQ2-Differential Expression",
                        subcategory="{analysis_name}"),
        done = touch(WORKDIR + "/results/plots/deseq2/comparison_{analysis_name}/diffexp/deg_bar/deseq2_plot_deg_bar.done"),
    params:
        analysis=lambda wildcards: next(a for a in config["params"]["deseq2"]["analyses"] if a["name"] == wildcards.analysis_name),
        padj_threshold = config["params"]["deseq2"]["padj_threshold"],
        lfc_threshold = config["params"]["deseq2"]["lfc_threshold"],
        funcs_script = DEA_FUNCS_SCRIPT,
    log:
        WORKDIR + "/logs/deseq2/diffexp/plot_deg_bar/deseq2_plot_deg_bar_{analysis_name}.log",
    conda:
        "../envs/deseq2_diffexp.yaml"
    threads:
        config["params"]["deseq2"]["threads"],
    script:
        "../scripts/deseq2_plot_deg_bar.R"

# ==================================================================================================================================================
#   Deseq2 Summarize Results
# ==================================================================================================================================================

rule deseq2_summarize_results:
    input:
        dds_subset = WORKDIR + "/results/deseq2/subset/dds_{analysis_name}.rds",
        rds_dir = WORKDIR + "/results/deseq2/diffexp/comparison_{analysis_name}",
        tcp_hits = config["ref"]["plant"]["tcp_hits"],
    output:
        tables_dir = report(    
                        directory(WORKDIR + "/results/tables/deseq2/comparison_{analysis_name}/summary"),
                        patterns=["{name}.csv"],
                        caption="../report/deseq2-summarize_results.rst",
                        category="DEseQ2-DE Summary",
                        subcategory="{analysis_name}"),
        done = touch(WORKDIR + "/results/tables/deseq2/comparison_{analysis_name}/summary/deseq2_summarize_results.done"),
    params:
        analysis=lambda wildcards: next(a for a in config["params"]["deseq2"]["analyses"] if a["name"] == wildcards.analysis_name),
        padj_threshold = config["params"]["deseq2"]["padj_threshold"],
        lfc_threshold = config["params"]["deseq2"]["lfc_threshold"],
    log:
        WORKDIR + "/logs/deseq2/diffexp/summarize/deseq2_summarize_results_{analysis_name}.log",
    conda:
        "../envs/deseq2_diffexp.yaml"
    threads:
        config["params"]["deseq2"]["threads"],
    script:
        "../scripts/deseq2_summarize_results.R"

# ==================================================================================================================================================
#   Deseq2 Venn Diagram
# ==================================================================================================================================================

rule deseq2_venn:
    input:
        dds_subset = WORKDIR + "/results/deseq2/subset/dds_{analysis_name}.rds",
        rds_dir = WORKDIR + "/results/deseq2/diffexp/comparison_{analysis_name}",
    output:
        plots_dir = report(
                        directory(WORKDIR + "/results/plots/deseq2/comparison_{analysis_name}/venn"),
                        patterns=["{name}.png"],
                        caption="../report/deseq2-venn.rst",
                        category="DEseQ2-Venn",
                        subcategory="{analysis_name}"),
        done = touch(WORKDIR + "/results/plots/deseq2/comparison_{analysis_name}/venn/deseq2_venn.done"),
    params:
        analysis=lambda wildcards: next(a for a in config["params"]["deseq2"]["analyses"] if a["name"] == wildcards.analysis_name),
        funcs_script = DEA_FUNCS_SCRIPT,
        padj_threshold = config["params"]["deseq2"]["padj_threshold"],
        lfc_threshold = config["params"]["deseq2"]["lfc_threshold"],
    log:
        WORKDIR + "/logs/deseq2/diffexp/venn/deseq2_venn_{analysis_name}.log",
    conda:
        "../envs/venndetail.yaml"
    threads:
        config["params"]["deseq2"]["threads"],
    script:
        "../scripts/deseq2_venn.R"

# =========== # =========== # =========== # =========== # =========== # =========== # ========== #

rule deseq2_venn_genotype:
    input:
        dds_subset = expand(
                        WORKDIR + "/results/deseq2/subset/dds_{analysis_name}.rds",
                        analysis_name = [a["name"] for a in config["params"]["deseq2"]["analyses"]]
                    ),
        rds_dir = expand(
                        WORKDIR + "/results/deseq2/diffexp/comparison_{analysis_name}",
                        analysis_name = [a["name"] for a in config["params"]["deseq2"]["analyses"]]
                    ),
    output:
        plots_dir = report(
                        directory(WORKDIR + "/results/plots/deseq2/venn/by_genotype"),
                        patterns=["{name}.png"],
                        caption="../report/deseq2-venn.rst",
                        category="DEseQ2-Venn",
                        subcategory="Venn by Genotype"),
        done = touch(WORKDIR + "/results/plots/deseq2/venn/by_genotype/deseq2_venn.by_genotype.done"),
    params:
        funcs_script = DEA_FUNCS_SCRIPT,
        padj_threshold = config["params"]["deseq2"]["padj_threshold"],
        lfc_threshold = config["params"]["deseq2"]["lfc_threshold"],
    log:
        WORKDIR + "/logs/deseq2/diffexp/venn/deseq2_venn.by_genotype.log",
    conda:
        "../envs/venndetail.yaml"
    threads:
        config["params"]["deseq2"]["threads"],
    script:
        "../scripts/deseq2_venn_genotype.R"

# ==================================================================================================================================================