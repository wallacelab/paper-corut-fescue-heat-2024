# ==================================================================================================================================================
#   Convert GO mappings to gmt
# ==================================================================================================================================================

rule convert_go_mappings_to_gmt:
    input:
        tcp_hits = config["ref"]["plant"]["tcp_hits"],
    output:
        gmt= "resources/ref/annotation/fescue.go_list.gmt"
    log:
       WORKDIR + "/logs/cemitool/convert_go_mappings_to_gmt.log"
    script:
        "../scripts/convert_go_mappings_to_gmt.py"

# ==================================================================================================================================================
#   Run CEMiTool
# ==================================================================================================================================================

rule cemitool:
    input:
        dds = WORKDIR + "/results/deseq2/dds.rds",
        dds_subset = WORKDIR + "/results/deseq2/subset/dds_{analysis_name}.rds",
        gmt = "resources/ref/annotation/fescue.go_list.gmt",
        interactions = config["ref"]["plant"]["interactions"],
    output:
        cemitool_outdir= report(
                            directory(WORKDIR + "/results/cemitool/{analysis_name}"),
                            patterns=["{name}.html"],
                            caption="../report/cemitool.rst",
                            category="CemiTool-Coexp. Network Analysis",
                            subcategory="{analysis_name}"),
    params:
        analysis= lambda wildcards: next(a for a in config["params"]["deseq2"]["analyses"] if a["name"] == wildcards.analysis_name),
        verbose= config["params"]["CEMiTool"]["verbose"],
        plot= config["params"]["CEMiTool"]["plot"],
        apply_vst= config["params"]["CEMiTool"]["apply_vst"],
        filter= config["params"]["CEMiTool"]["filter"],
        count_filter= config["params"]["CEMiTool"]["count_filter"]["activate"],
        min_counts= config["params"]["CEMiTool"]["count_filter"]["min_counts"],
        min_samples= config["params"]["CEMiTool"]["count_filter"]["min_samples"],
        network_type= config["params"]["CEMiTool"]["network_type"],
        cor_method= config["params"]["CEMiTool"]["cor_method"],
    conda:
        "../envs/cemitool.yaml"
    threads: 
        config["params"]["CEMiTool"]["threads"],
    log:
        WORKDIR + "/logs/cemitool/{analysis_name}/run_cemitool.log"
    script:
        "../scripts/cemitool.R"

# ==================================================================================================================================================

    
