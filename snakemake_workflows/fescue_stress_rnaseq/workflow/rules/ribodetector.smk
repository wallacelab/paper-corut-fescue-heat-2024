# =================================================================================================
#     Remove rRNA reads from the trimmed reads using ribodetector
# =================================================================================================

rule ribodetector:
    input:
        r1= WORKDIR + "/results/trimmed/{sample}-{unit}.fastq.gz",
        stats= WORKDIR + "/results/stats/seqkit/{sample}-{unit}.stats.tsv",
    output:
        norrna= WORKDIR + "/results/ribodetector/{sample}-{unit}_1_trimmed.norrna.fq",
    params:
        memory= config["params"]["ribodetector"]["mem"],
        read_length= lambda wildcards, input: get_avg_read_length(input.stats),
    conda:
        "../envs/ribodetector.yaml"
    log:
        WORKDIR + "/logs/ribodetector/{sample}-{unit}.ribodetector.log"
    threads:
        config["params"]["ribodetector"]["threads"]
    shell:
        """
        ribodetector_cpu -t {threads} \
        -l {params.read_length} \
        -i {input.r1} \
        -e rrna \
        -o {output.norrna}
        """

# =================================================================================================
