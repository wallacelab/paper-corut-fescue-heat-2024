# =================================================================================================
#     Import Libraries
# =================================================================================================

import glob
import os
import json
import pandas as pd
import socket, platform
import math
from snakemake.io import expand
from snakemake.utils import R
from snakemake.utils import min_version
from snakemake.utils import validate
from snakemake.io import glob_wildcards
import re
from os.path import join, basename, dirname
import pathlib
from os import path
from datetime import datetime

# =================================================================================================
#     Set Workdir
# =================================================================================================

WORKDIR = config["workdir"]["path"]

# =================================================================================================
#     Sample and Phenotype Sheets + Wildcard Constraints
# =================================================================================================

# ===================== Samples ===================== #

# Read samples sheet (samples.tsv)
samples = (
    pd.read_csv(
        config["samples"],
        sep="\t",
        dtype=str,
        comment="#",
    )
    .set_index(["sample_name"], drop=False)
    .sort_index()
)
samples.index.names = ["sample_id"]

# ====================== Units ====================== #

# Read units sheet (units.tsv)
units = (
    pd.read_csv(
        config["units"], 
        sep="\t", 
        dtype=str, 
        comment="#",
    )
    .set_index(["sample_name", "unit_name"], drop=False)
    .sort_index()
)
units.index.names = ["sample_id", "unit_id"]

# Make sure indeces always str type
units.index = units.index.set_levels(
    [i.astype(str) for i in units.index.levels]
)

wildcard_constraints:
    sample="|".join(samples.index),
    unit="|".join(units["unit_name"]),

# =================================================================================================
#     Pipeline User Output
# =================================================================================================

author = "Adnan Kivanc Corut"
date_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
snake_version = snakemake.__version__
py_version= sys.version.split(' ')[0]
snakefile_path= workflow.snakefile
base_dir= workflow.basedir
work_dir= os.getcwd()
config_file= ", ".join(workflow.configfiles)

# Helpful messages
logger.info("# ================================================================================== #")
logger.info("")
logger.info(f"     Date:            {date_time}")
logger.info(f"     Author:          {author}")
logger.info(f"     Snakemake version:          {snake_version}")
logger.info(f"     Python version:             {py_version}")
logger.info("")
logger.info(f"     Snakefile:          {snakefile_path}")
logger.info(f"     Base directory:     {base_dir}")
logger.info(f"     Working directory:  {work_dir}")
logger.info(f"     Results directory:  {WORKDIR}")
logger.info(f"     Config files:       {config_file}")
logger.info("")
logger.info("# ================================================================================== #")
logger.info("")

# =================================================================================================
#     Common Helper Functions
# =================================================================================================

## Source: https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/blob/main/workflow/rules/common.smk
def get_fastqs(wildcards):
    """Get raw FASTQ files from unit sheet."""
    return units.loc[(wildcards.sample, wildcards.unit), "fq1"]

## Source: https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/blob/main/workflow/rules/common.smk
def get_trimmed(wildcards):
    """Get trimmed FASTQ files."""
    return expand(WORKDIR + "/results/trimmed/{sample}-{unit}.fastq.gz", **wildcards)

# Function to get the trimmed and rRNA removed FASTQ files from ribodetector
def get_rrna_removed(wildcards):
    """Get trimmed FASTQ files."""
    return expand(WORKDIR + "/results/ribodetector/{sample}-{unit}_1_trimmed.norrna.fq", **wildcards)

# Function to get the input files for salmon
def get_salmon_input(wildcards):
    if config["settings"]["rrna_removal"]["activate"]:
        return expand(WORKDIR + "/results/ribodetector/{sample}-{unit}_1_trimmed.norrna.fq", **wildcards)
    else:
        return expand(WORKDIR + "/results/trimmed/{sample}-{unit}.fastq.gz", **wildcards)

# Function to get salmon parameters from config file
def salmon_params(wildcards, input):
    """Get Salmon params."""
    extra = config["params"]["salmon"]["quant_extra"]
    b= config["params"]["salmon"]["num_bootstrap"]
    fld_mean= config["params"]["salmon"]["fld_mean"]
    fld_sd= config["params"]["salmon"]["fld_sd"]
    extra += (
        " --numBootstraps {bootstrap} " "--fldMean {fragment_len_mean} " "--fldSD {fragment_len_sd}" 
    ).format(bootstrap= b, fragment_len_mean= fld_mean, fragment_len_sd= fld_sd)
    return extra

# Get the average read length from the seqkit stats file
def get_avg_read_length(filename):
    # Read the file into a DataFrame
    stats = pd.read_csv(filename, sep='\t')

    # Get the 'avg_len' value from the first (and presumably only) row
    avg_len = round(stats['avg_len'].values[0])

    return avg_len

# =================================================================================================
#     Target Ouput Function
# =================================================================================================

def get_target_output(wildcards):
    """
    Get all requested inputs (target outputs) for rule all.
    """

    target_output = []
    #### QC #####
    target_output.extend(
        expand(
            WORKDIR + "/results/qc/multiqc.html"
        )
    ),
    target_output.extend(
        expand(
            WORKDIR + "/results/qc/multiqc_trim.3prime.html"
        )
    ),
    target_output.extend(
        expand(
            WORKDIR + "/results/qc/multiqc_ribodetector.html",
        )
    ),
    # #### Kraken #####
    # target_output.extend(
    #     expand(
    #         [
    #             WORKDIR + "/results/plots/kraken/final_kraken_report.plot_by_condition.pdf",
    #         ]
    #     )
    # ),
    # target_output.extend(
    #     expand(
    #         [
    #             WORKDIR + "/results/qc/multiqc_kraken.html",
    #         ]
    #     )
    # ),
    # if config["settings"]["kraken"]["use_custom_db"]:
    #     target_output.extend(
    #         expand(
    #             [
    #                 WORKDIR + "/results/qc/multiqc_kraken_{db_name}.html",
    #             ],
    #             db_name= config["params"]["kraken2"]["db_name"]
    #         )
    #     ),
    #     target_output.extend(
    #         expand(
    #             [
    #                 WORKDIR + "/results/kraken_{db_name}/krona/multi-krona.html",
    #             ],
    #             db_name= config["params"]["kraken2"]["db_name"]
    #         )
    #     ),
    #### Bowtie2 #####
    target_output.extend(
        expand(
            WORKDIR + "/results/qc/multiqc_bowtie2.html",
        )
    ),
    #### Salmon #####
    target_output.extend(
        expand(
            WORKDIR + "/results/qc/multiqc_salmon.html"
        )
    ),
    target_output.extend(
        expand(
            WORKDIR + "/results/salmon/{u.sample_name}-{u.unit_name}/quant_bootstraps.tsv",
            u=units.itertuples()
        )
    ),
    target_output.extend(
        expand(
            [   
                WORKDIR + "/results/plots/coverage/{u.sample_name}-{u.unit_name}.salmon_mapping.coverage.png",
                WORKDIR + "/results/plots/coverage/{u.sample_name}-{u.unit_name}.salmon_mapping.coverage_length_ranges.png",
            ],
            u=units.itertuples()
        )
    ),
    target_output.extend(
        expand(
            WORKDIR + "/results/tables/catchSalmon/catchSalmon.RDS"
        )
    ),
    #### DESeq2 #####
    target_output.extend(
        expand(
            WORKDIR + "/results/plots/deseq2/dispersion_est/dispersion_estimates.png",
        )
    ),
    # get all the variables to plot a PCA for
    pca_variables = list(config["params"]["deseq2"]["variables_of_interest"])
    target_output.extend(
        expand(
            [   
                WORKDIR + "/results/plots/deseq2/pca/pca.{variable}.png",
                WORKDIR + "/results/plots/deseq2/pca/pca.{variable}.combined.png",
            ],
            variable=pca_variables)
    )
    target_output.extend(
        expand(
            [
                WORKDIR + "/results/plots/deseq2/comparison_{analysis}/eda/pca.{analysis}.combined.png",
                WORKDIR + "/results/plots/deseq2/comparison_{analysis}/eda/{analysis}.eigencorplot.png",
                WORKDIR + "/results/plots/deseq2/comparison_{analysis}/eda/{analysis}.cooks_boxplot.png",
                WORKDIR + "/results/plots/deseq2/comparison_{analysis}/eda/{analysis}.disp_est.png",
                WORKDIR + "/results/plots/deseq2/comparison_{analysis}/eda/{analysis}.density_plot.png",
                WORKDIR + "/results/plots/deseq2/comparison_{analysis}/eda/{analysis}.sample_dist_heatmap.png",
                WORKDIR + "/results/plots/deseq2/comparison_{analysis}/eda/{analysis}.count_matrix_heatmap.png", 
                WORKDIR + "/results/plots/deseq2/comparison_{analysis}/eda/{analysis}.gene_cluster_heatmap.png",
                WORKDIR + "/results/plots/deseq2/comparison_{analysis}/eda/deseq2_exploratory_plots.done",
            ],
            analysis=[a['name'] for a in config["params"]["deseq2"]["analyses"]])
    )
    target_output.extend(
        expand(
            [
                WORKDIR + "/results/tables/deseq2/comparison_{analysis}/diffexp",
                WORKDIR + "/results/deseq2/diffexp/comparison_{analysis}",
            ],
            analysis=[a['name'] for a in config["params"]["deseq2"]["analyses"]])
    )
    target_output.extend(
        expand(
            [
                WORKDIR + "/results/plots/deseq2/comparison_{analysis}/diffexp/volcano/deseq2_plot_volcano.done",
                WORKDIR + "/results/plots/deseq2/comparison_{analysis}/diffexp/ma/deseq2_plot_ma.done",
                WORKDIR + "/results/plots/deseq2/comparison_{analysis}/diffexp/lfc_bar/deseq2_plot_lfc_bar.done",
                WORKDIR + "/results/plots/deseq2/comparison_{analysis}/diffexp/deg_bar/deseq2_plot_deg_bar.done",
                WORKDIR + "/results/plots/deseq2/comparison_{analysis}/diffexp/exp_heatmap/deseq2_plot_exp_heatmap.done",
            ],
            analysis=[a['name'] for a in config["params"]["deseq2"]["analyses"]])
    )
    target_output.extend(
        expand(
            [
                WORKDIR + "/results/tables/deseq2/comparison_{analysis}/summary/deseq2_summarize_results.done",
            ],
            analysis=[a['name'] for a in config["params"]["deseq2"]["analyses"]])
    )
    target_output.extend(
        expand(
            [
                WORKDIR + "/results/plots/deseq2/comparison_{analysis}/venn/deseq2_venn.done",
                WORKDIR + "/results/plots/deseq2/venn/by_genotype/deseq2_venn.by_genotype.done"
            ],
            analysis=[a['name'] for a in config["params"]["deseq2"]["analyses"]])
    )
    #### Functional Enrichment #####
    target_output.extend(
        expand(
            [
                WORKDIR + "/results/tables/enrichment/comparison_{analysis}/go",
                WORKDIR + "/results/plots/enrichment/comparison_{analysis}/go",
            ],
            analysis=[a['name'] for a in config["params"]["deseq2"]["analyses"]])
    )
    target_output.extend(
        expand(
            [
                WORKDIR + "/results/tables/enrichment/comparison_{analysis}/kegg",
                WORKDIR + "/results/plots/enrichment/comparison_{analysis}/kegg",
            ],
            analysis=[a['name'] for a in config["params"]["deseq2"]["analyses"]])
    )
    target_output.extend(
        expand(
            [
                WORKDIR + "/results/tables/enrichment/comparison_{analysis}/wikipathways",
                WORKDIR + "/results/plots/enrichment/comparison_{analysis}/wikipathways"
            ],
            analysis=[a['name'] for a in config["params"]["deseq2"]["analyses"]])
    )
    target_output.extend(
        expand(
            [
                WORKDIR + "/results/plots/enrichment/comparison_{analysis}/compareCluster",
                WORKDIR + "/results/go_enrichment/compareCluster/comparison_{analysis}"
            ],
            analysis=[a['name'] for a in config["params"]["deseq2"]["analyses"]])
    )
    target_output.extend(
        expand(
            WORKDIR + "/results/plots/enrichment/compareCluster_all/compareCluster_all.done"
        )
    )
    target_output.extend(
        expand(
            WORKDIR + "/results/plots/enrichment/go_enrichment_unique_degs/go_enrichment_unique_degs.done"
        )
    )
    target_output.extend(
        expand(
            WORKDIR + "/results/plots/enrichment/compare_heat_stress_effect/compare_heat_stress_effect.done"
        )
    )
    #### CEMiTool #####
    target_output.extend(
        expand(
            [
                WORKDIR + "/results/cemitool/{analysis}"
            ],
            analysis=[a['name'] for a in config["params"]["deseq2"]["analyses"]])   
    )

    return target_output