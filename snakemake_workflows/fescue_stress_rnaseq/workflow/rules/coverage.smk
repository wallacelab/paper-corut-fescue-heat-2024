# =================================================================================================
#     Convert SAM to BAM
# =================================================================================================

rule sam_to_bam:
    input:
        WORKDIR + "/results/salmon/{sample}-{unit}/{sample}-{unit}.salmon_mapping.sam"
    output:
        temp(WORKDIR + "/results/coverage/{sample}-{unit}/{sample}-{unit}.salmon_mapping.bam")
    conda:
        "../envs/coverage.yaml"
    threads: 8
    shell:
        "samtools view -@ {threads} -bS {input} > {output}"

# =================================================================================================
#     Sort BAM
# =================================================================================================

rule sort_bam:
    input:
        WORKDIR + "/results/coverage/{sample}-{unit}/{sample}-{unit}.salmon_mapping.bam"
    output:
        WORKDIR + "/results/coverage/{sample}-{unit}/{sample}-{unit}.salmon_mapping.sorted.bam"
    conda:
        "../envs/coverage.yaml"
    threads: 8
    shell:
        "samtools sort -@ {threads} {input} -o {output}"

# =================================================================================================
#     Index BAM
# =================================================================================================

rule index_bam:
    input:
        WORKDIR + "/results/coverage/{sample}-{unit}/{sample}-{unit}.salmon_mapping.sorted.bam"
    output:
        WORKDIR + "/results/coverage/{sample}-{unit}/{sample}-{unit}.salmon_mapping.sorted.bam.bai"
    conda:
        "../envs/coverage.yaml"
    threads: 8
    shell:
        "samtools index -@ 8 {input}"

# =================================================================================================
#    Calculate coverage
# =================================================================================================

rule calculate_coverage:
    input:
        WORKDIR + "/results/coverage/{sample}-{unit}/{sample}-{unit}.salmon_mapping.sorted.bam"
    output:
        WORKDIR + "/results/coverage/{sample}-{unit}.salmon_mapping.bedgraph"
    conda:
        "../envs/coverage.yaml"
    threads: 8
    shell:
        "bedtools genomecov -ibam {input} -bg -split > {output}"

# =================================================================================================
#    Plot coverage
# =================================================================================================

rule plot_coverage:
    input:
        WORKDIR + "/results/coverage/{sample}-{unit}.salmon_mapping.bedgraph"
    output:
        avg_plot = WORKDIR + "/results/plots/coverage/{sample}-{unit}.salmon_mapping.coverage.png",
        length_range_plot = WORKDIR + "/results/plots/coverage/{sample}-{unit}.salmon_mapping.coverage_length_ranges.png"
    params:
        normalize = True
    conda:
        "../envs/matplotlib.yaml"
    log:
        WORKDIR + "/logs/coverage/{sample}-{unit}.salmon_mapping.coverage.log"
    threads: 8
    script:
        "../scripts/plot_coverage.py"

# =================================================================================================