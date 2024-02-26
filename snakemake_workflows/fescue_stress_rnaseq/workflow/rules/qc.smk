# =================================================================================================
#     Generate QC Stats Using FastQC on Raw Reads
# =================================================================================================

rule fastqc:
    input:
        get_fastqs
    output:
        html= WORKDIR + "/results/qc/fastqc/{sample}-{unit}_fastqc.html",
        zip= WORKDIR + "/results/qc/fastqc/{sample}-{unit}_fastqc.zip",
    log:
        WORKDIR + "/logs/qc/fastqc/{sample}-{unit}.log",
    threads:
        config["params"]["fastqc"]["threads"]
    message: 
        "Performing quality control analysis using FastQC on the following file: {input}"
    wrapper:
        "v1.20.0/bio/fastqc"

# =================================================================================================
#     MultiQC
# =================================================================================================

rule multiqc:
    input:
        expand(
            WORKDIR + "/results/qc/fastqc/{u.sample_name}-{u.unit_name}_fastqc.zip",
            u=units.itertuples(),
        )
    output:
        report(
            WORKDIR + "/results/qc/multiqc.html",
            caption="../report/multiqc.rst",
            category="Quality Control",
        ),
    log:
        WORKDIR + "/logs/qc/multiqc/multiqc.log",
    params:
        extra= "-v -d --interactive"
    message: 
        "Performing MultiQC on the FastQC results..."
    wrapper:
        "v1.20.0/bio/multiqc"

# =================================================================================================
#     Generate QC Stats Using FastQC on Trimmed Reads
# =================================================================================================

rule fastqc_trim:
    input:
        get_trimmed,
    output:
        html= WORKDIR + "/results/qc/fastqc_trim/{sample}-{unit}_trimmed_fastqc.html",
        zip= WORKDIR + "/results/qc/fastqc_trim/{sample}-{unit}_trimmed_fastqc.zip",
    log:
        WORKDIR + "/logs/qc/fastqc/trim/{sample}-{unit}_trimmed.log",
    threads:
        config["params"]["fastqc"]["threads"]
    message: 
        "Performing quality control analysis using FastQC on trimmed fastq file: {input}"
    wrapper:
        "v1.20.0/bio/fastqc"

# =================================================================================================
#     MultiQC Trim
# =================================================================================================

rule multiqc_trim:
    input:
        expand(
            WORKDIR + "/results/qc/fastqc_trim/{u.sample_name}-{u.unit_name}_trimmed_fastqc.zip",
            u=units.itertuples(),
        ),
        expand(
            WORKDIR + "/results/trimmed/{u.sample_name}-{u.unit_name}.qc.txt",
            u=units.itertuples(),
        ),
    output:
        report(
            WORKDIR + "/results/qc/multiqc_trim.3prime.html",
            caption="../report/multiqc.rst",
            category="Quality Control",
        ),
    log:
        WORKDIR + "/logs/qc/multiqc/multiqc_trim.3prime.log",
    params:
        extra= "-v -d --interactive"
    message: 
        "Performing MultiQC on the FastQC of trimmed reads..."
    wrapper:
        "v1.20.0/bio/multiqc"

# =================================================================================================
#     Generate QC Stats Using Ribodetector Results
# =================================================================================================

rule fastqc_norrna:
    input:
        WORKDIR + "/results/ribodetector/{sample}-{unit}_1_trimmed.norrna.fq",
    output:
        html= WORKDIR + "/results/ribodetector/{sample}-{unit}_1_trimmed.norrna_fastqc.html",
        zip= WORKDIR + "/results/ribodetector/{sample}-{unit}_1_trimmed.norrna_fastqc.zip",
    log:
        WORKDIR + "/logs/qc/fastqc/fastqc_norrna/{sample}-{unit}_fastqc_norrna.log",
    threads:
        config["params"]["fastqc"]["threads"]
    message: 
        "Performing quality control analysis using FastQC on nonrrna fastq file: {input}"
    wrapper:
        "v1.20.0/bio/fastqc"

# =================================================================================================
#     MultiQC Ribodetector
# =================================================================================================

rule multiqc_ribodetector:
    input:
        expand(
            WORKDIR + "/results/ribodetector/{u.sample_name}-{u.unit_name}_1_trimmed.norrna_fastqc.zip",
            u=units.itertuples(),
        )
    output:
        report(
            WORKDIR + "/results/qc/multiqc_ribodetector.html",
            caption="../report/multiqc.rst",
            category="Quality Control",
        ),
    log:
        WORKDIR + "/logs/qc/multiqc/multiqc_ribodetector.log",
    params:
        extra= "-v -d --interactive"
    message: 
        "Performing MultiQC on the kallisto results..."
    wrapper:
        "v1.20.0/bio/multiqc"

# =================================================================================================
#     MultiQC STAR
# =================================================================================================

rule multiqc_star:
    input:
        expand(
            WORKDIR + "/results/star/{u.sample_name}-{u.unit_name}/Aligned.sortedByCoord.out.bam",
            u=units.itertuples(),
        ),
        expand(
            WORKDIR + "/results/star/{u.sample_name}-{u.unit_name}/Log.final.out",
            u=units.itertuples(),
        )
    output:
        report(
            WORKDIR + "/results/qc/multiqc_star.html",
            caption="../report/multiqc.rst",
            category="Quality Control",
        ),
    log:
        WORKDIR + "/logs/qc/multiqc/multiqc_star.log",
    params:
        extra= "-v -d --interactive"
    message: 
        "Performing MultiQC on STAR alignment results..."
    wrapper:
        "v1.20.0/bio/multiqc"

# =================================================================================================
#     MultiQC Kallisto
# =================================================================================================

rule multiqc_kallisto:
    input:
        expand(
            WORKDIR + "/logs/kallisto/quant/{u.sample_name}-{u.unit_name}.kallisto_quant.3prime.log",
            u=units.itertuples(),
        )
    output:
        report(
            WORKDIR + "/results/qc/multiqc_kallisto.3prime.html",
            caption="../report/multiqc.rst",
            category="Quality Control",
        ),
    log:
        WORKDIR + "/logs/qc/multiqc/multiqc_kallisto.3prime.log",
    params:
        extra= "-v -d --interactive"
    message: 
        "Performing MultiQC on the kallisto results..."
    wrapper:
        "v1.20.0/bio/multiqc"

# =================================================================================================
#     MultiQC Kallisto
# =================================================================================================

rule multiqc_salmon:
    input:
        expand(
            WORKDIR + "/results/salmon/{u.sample_name}-{u.unit_name}",
            u=units.itertuples(),
        )
    output:
        report(
            WORKDIR + "/results/qc/multiqc_salmon.html",
            caption="../report/multiqc.rst",
            category="Quality Control",
        ),
    log:
        WORKDIR + "/logs/qc/multiqc/multiqc_salmon.log",
    params:
        extra= "-v -d --interactive"
    message: 
        "Performing MultiQC on the Salmon results..."
    wrapper:
        "v1.20.0/bio/multiqc"

# =================================================================================================
#     MultiQC Kraken2
# =================================================================================================

rule multiqc_kraken:
    input:
        expand(
            WORKDIR + "/results/kraken/{u.sample_name}-{u.unit_name}.kraken.report",
            u=units.itertuples(),
        )
    output:
        report(
            WORKDIR + "/results/qc/multiqc_kraken.html",
            caption="../report/multiqc.rst",
            category="Quality Control",
        ),
    log:
        WORKDIR + "/logs/qc/multiqc/multiqc_kraken.log",
    params:
        extra= "-v -d --interactive"
    message: 
        "Performing MultiQC on the Kraken2 results..."
    wrapper:
        "v1.20.0/bio/multiqc"

rule multiqc_kraken_custom:
    input:
        expand(
            WORKDIR + "/results/kraken_{db_name}/{u.sample_name}-{u.unit_name}.kraken.report",
            u=units.itertuples(), 
            db_name=config["params"]["kraken2"]["db_name"]
        )
    output:
        report(
            WORKDIR + "/results/qc/multiqc_kraken_{db_name}.html",
            caption="../report/multiqc.rst",
            category="Quality Control",
        ),
    log:
        WORKDIR + "/logs/qc/multiqc/multiqc_kraken_{db_name}.log",
    params:
        extra= "-v -d --interactive"
    message: 
        "Performing MultiQC on the Kraken2 results..."
    wrapper:
        "v1.20.0/bio/multiqc"

# =================================================================================================
#      Qualimap
# =================================================================================================
rule qualimap:
    input:
        # BAM aligned
        bam= WORKDIR + "/results/coverage/{sample}-{unit}/{sample}-{unit}.salmon_mapping.sorted.bam",
    output:
        directory(WORKDIR + "/results/qc/qualimap/{sample}-{unit}"),
    log:
        WORKDIR + "/logs/qualimap/bamqc/{sample}-{unit}.log",
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=4096,
    wrapper:
        "v2.6.0/bio/qualimap/bamqc"

# =================================================================================================
#      preseq
# =================================================================================================

rule preseq_lc_extrap_bam:
    input:
        WORKDIR + "/results/coverage/{sample}-{unit}/{sample}-{unit}.salmon_mapping.sorted.bam",
    output:
        WORKDIR + "/results/test_bam/{sample}-{unit}.lc_extrap"
    params:
        "-v"   #optional parameters
    log:
        WORKDIR + "/logs/test_bam/{sample}-{unit}.log"
    wrapper:
        "v2.6.0/bio/preseq/lc_extrap"

# =================================================================================================
#      RSeQC
# =================================================================================================

rule rseqc_stat:
    input:
        WORKDIR + "/results/coverage/{sample}-{unit}/{sample}-{unit}.salmon_mapping.sorted.bam",
    output:
        WORKDIR + "/results/qc/rseqc/{sample}-{unit}.stats.txt",
    log:
        WORKDIR + "/logs/rseqc/rseqc_stat/{sample}-{unit}.log",
    conda:
        "../envs/rseqc.yaml"
    shell:
        "bam_stat.py -i {input} > {output} 2> {log}"

rule rseqc_readdup:
    input:
        WORKDIR + "/results/coverage/{sample}-{unit}/{sample}-{unit}.salmon_mapping.sorted.bam",
    output:
        WORKDIR + "/results/qc/rseqc/{sample}-{unit}.readdup.DupRate_plot.pdf",
    log:
         WORKDIR + "/logs/rseqc/rseqc_readdup/{sample}_{unit}.log",
    params:
        prefix=lambda w, output: output[0].replace(".DupRate_plot.pdf", ""),
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_duplication.py -i {input} -o {params.prefix} > {log} 2>&1"


rule rseqc_readgc:
    input:
        WORKDIR + "/results/coverage/{sample}-{unit}/{sample}-{unit}.salmon_mapping.sorted.bam",
    output:
        WORKDIR + "/results/qc/rseqc/{sample}-{unit}.readgc.GC_plot.pdf",
    log:
        WORKDIR + "/logs/rseqc/rseqc_readgc/{sample}_{unit}.log",
    params:
        prefix=lambda w, output: output[0].replace(".GC_plot.pdf", ""),
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_GC.py -i {input} -o {params.prefix} > {log} 2>&1"

rule rseqc_geneBody_coverage:
    input:
        bam=WORKDIR + "/results/coverage/{sample}-{unit}/{sample}-{unit}.salmon_mapping.sorted.bam",
        bai=WORKDIR + "/results/coverage/{sample}-{unit}/{sample}-{unit}.salmon_mapping.sorted.bam.bai",
        bed="/scratch/ac32082/03.FescueHeatStress/fescue_stress_scripts/nextflow_workflows/map_transcriptome/work/ff/31aa19b3c3f054b4ab13027a8834d1/results/bed/fescue_tcp.map_to_refv01.bed"
    output:
        WORKDIR + "/results/qc/rseqc/{sample}-{unit}.geneBody_coverage.png",
    params:
        prefix=lambda w, output: output[0].replace(".png", ""),
    log:
        WORKDIR + "/logs/rseqc/rseqc_geneBody_coverage/{sample}_{unit}.log",
    conda:
        "../envs/rseqc.yaml"
    shell:
        "geneBody_coverage.py -r {input.bed} -i {input.bam} -f png -o {params.prefix} > {log} 2>&1"

# =================================================================================================
#      MultiQC BAM
# =================================================================================================

rule multiqc_bam:
    input:
        expand(
            WORKDIR + "/results/qc/qualimap/{u.sample_name}-{u.unit_name}",
            u=units.itertuples(),
        ),
        expand(
            WORKDIR + "/results/test_bam/{u.sample_name}-{u.unit_name}.lc_extrap",
            u=units.itertuples(),
        ),
        expand(
            WORKDIR + "/results/qc/rseqc/{u.sample_name}-{u.unit_name}.stats.txt",
            u=units.itertuples(),
        ),
        expand(
            WORKDIR + "/results/qc/rseqc/{u.sample_name}-{u.unit_name}.readdup.DupRate_plot.pdf",
            u=units.itertuples(),
        ),
        expand(
            WORKDIR + "/results/qc/rseqc/{u.sample_name}-{u.unit_name}.readgc.GC_plot.pdf",
            u=units.itertuples(),
        ),
        # expand(
        #     WORKDIR + "/results/qc/rseqc/{u.sample_name}-{u.unit_name}.geneBody_coverage.png",
        #     u=units.itertuples(),
        # ),
    output:
        report(
            WORKDIR + "/results/qc/multiqc_bam.html",
            caption="../report/multiqc.rst",
            category="Quality Control",
        ),
    log:
        WORKDIR + "/logs/qc/multiqc/multiqc_bam.log",
    params:
        extra= "-v -d --interactive"
    message: 
        "Performing MultiQC on the BAM files.."
    wrapper:
        "v1.20.0/bio/multiqc"

# =================================================================================================
#      MultiQC Bowtie2
# =================================================================================================

rule multiqc_bowtie2:
    input:
        expand(
            WORKDIR + "/logs/bowtie2/{u.sample_name}-{u.unit_name}.log",
            u=units.itertuples(),
        ),
        expand(
            WORKDIR + "/results/bowtie2/{u.sample_name}-{u.unit_name}/{u.sample_name}-{u.unit_name}.bowtie2.align.sorted.bam.bai",
            u=units.itertuples(),
        )
    output:
        report(
            WORKDIR + "/results/qc/multiqc_bowtie2.html",
            caption="../report/multiqc.rst",
            category="Quality Control",
        ),
    log:
        WORKDIR + "/logs/qc/multiqc/multiqc_bowtie2.log",
    params:
        extra= "-v -d --interactive"
    message: 
        "Performing MultiQC on the Bowtie2 results..."
    wrapper:
        "v1.20.0/bio/multiqc"

# =================================================================================================