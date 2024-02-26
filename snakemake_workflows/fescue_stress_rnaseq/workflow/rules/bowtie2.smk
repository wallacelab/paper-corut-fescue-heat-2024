# =================================================================================================
#     Bowtie2 Indexing
# =================================================================================================

rule bowtie2_build:
    input:
        ref="resources/ref/genome/plant.genome.fasta",
    output:
        multiext(
            "resources/ref/genome/bowtie2_index/plant.genome",
            ".1.bt2l",
            ".2.bt2l",
            ".3.bt2l",
            ".4.bt2l",
            ".rev.1.bt2l",
            ".rev.2.bt2l",
        ),
    log:
        WORKDIR + "/logs/bowtie2/bowtie2_build.log",
    params:
        extra="",  # optional parameters
    threads: 64
    wrapper:
        "v2.6.0/bio/bowtie2/build"

# =================================================================================================
#     Bowtie2 Alignment
# =================================================================================================

rule bowtie2_align:
    input:
        sample=get_rrna_removed,
        idx=multiext(
            "resources/ref/genome/bowtie2_index/plant.genome",
            ".1.bt2l",
            ".2.bt2l",
            ".3.bt2l",
            ".4.bt2l",
            ".rev.1.bt2l",
            ".rev.2.bt2l",
        ),
    output:
        temp(WORKDIR + "/results/bowtie2/{sample}-{unit}/{sample}-{unit}.bowtie2.align.bam")
    log:
        WORKDIR + "/logs/bowtie2/{sample}-{unit}.log",
    params:
        extra="--end-to-end --no-unal",  # optional parameters
    threads: 8  # Use at least two threads
    wrapper:
        "v2.6.0/bio/bowtie2/align"

# =================================================================================================
#     Samtools Sorting
# =================================================================================================

rule samtools_sort:
    input:
        WORKDIR + "/results/bowtie2/{sample}-{unit}/{sample}-{unit}.bowtie2.align.bam",
    output:
        WORKDIR + "/results/bowtie2/{sample}-{unit}/{sample}-{unit}.bowtie2.align.sorted.bam",
    log:
        WORKDIR + "/logs/samtools/{sample}-{unit}.log",
    params:
        extra="-m 12G",  # optional parameters
    threads: 8
    wrapper:
        "v2.6.0/bio/samtools/sort"

# =================================================================================================
#     Samtools Indexing
# =================================================================================================

rule samtools_index:
    input:
        WORKDIR + "/results/bowtie2/{sample}-{unit}/{sample}-{unit}.bowtie2.align.sorted.bam",
    output:
        WORKDIR + "/results/bowtie2/{sample}-{unit}/{sample}-{unit}.bowtie2.align.sorted.bam.bai",
    log:
        WORKDIR + "/logs/samtools/{sample}-{unit}.log",
    params:
        extra="-m 12G",  # optional parameters
    threads: 8
    wrapper:
        "v2.6.0/bio/samtools/index"

# =================================================================================================
