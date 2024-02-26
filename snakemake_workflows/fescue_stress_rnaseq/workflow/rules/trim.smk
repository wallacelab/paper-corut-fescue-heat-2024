# =================================================================================================
#     Create Symlinks
# =================================================================================================

rule create_symlink:
    input:
        get_fastqs,
    output:
        WORKDIR + "/results/reads/{sample}-{unit}_1.fastq.gz",
    message:
        "Creating symbolic links for fastq file..."
    log:
        WORKDIR + "/logs/create_symlink/{sample}-{unit}_1.fastq.gz.create_symlink.log",
    threads: 1
    shell:
        """
        echo Working on fastq file: {input}
        echo Symlinking {input} to {output}
        ln -rs {input} {output} 2> {log}
        """

# =================================================================================================
#     Trimming Single End Reads
# =================================================================================================

# rule trim_galore_se:
#     input:
#         WORKDIR + "/results/reads/{sample}-{unit}_1.fastq.gz"
#     output:
#         WORKDIR + "/results/trimmed/{sample}-{unit}_1_trimmed.fq.gz",
#         WORKDIR + "/results/trimmed/{sample}-{unit}_1.fastq.gz_trimming_report.txt",
#     params:
#         extra="--three_prime_clip_R1 12",
#     message: 
#         "Performing trimming using trim_galore on the following fastq file: {input}"
#     log:
#         WORKDIR + "/logs/trim_galore/{sample}-{unit}.log",
#     wrapper:
#         "v1.20.0/bio/trim_galore/se"

# =================================================================================================
#     Trimming Read 1 in QuantSeq FWD libraries
# =================================================================================================

# Adapted from: https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth/blob/main/workflow/rules/trim_3prime.smk
# https://faqs.lexogen.com/faq/what-is-the-adapter-sequence-i-need-to-use-for-t-1

# Rule cutadapt1 checks and removes poly-A tails and sequence qualilty score <20.
rule cutadapt1:
    input:
        WORKDIR + "/results/reads/{sample}-{unit}_1.fastq.gz",
    output:
        fastq= WORKDIR + "/results/trimmed/{sample}-{unit}.1.1.fastq.gz",
        qc= WORKDIR + "/results/trimmed/{sample}-{unit}.1.1.qc.txt",
    params:
        extra=lambda w: str(
            "-m 20 -O 20 -a " "polyA=A{20}" " -a " "QUALITY=G{20}" " -n 2"
        ),
    threads: config["params"]["cutadapt"]["threads"],
    log:
        WORKDIR + "/logs/cutadapt/{sample}-{unit}.1.1.log",
    wrapper:
        "v1.31.1/bio/cutadapt/se"

# Rule cutadapt2 checks if reads contains the polyA-stretch + adapter at the 3'
rule cutadapt2:
    input:
        fastq= WORKDIR + "/results/trimmed/{sample}-{unit}.1.1.fastq.gz",
    output:
        fastq= WORKDIR + "/results/trimmed/{sample}-{unit}.2.1.fastq.gz",
        qc= WORKDIR + "/results/trimmed/{sample}-{unit}.2.1.qc.txt",
    params:
        adapters=lambda w: str(
            '-m 20 -O 3 --nextseq-trim=10 -a "'
            'A{18}AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=3;max_error_rate=0.100000"'
            ""
        ),
    threads: config["params"]["cutadapt"]["threads"],
    log:
        WORKDIR + "/logs/cutadapt/{sample}-{unit}.2.1.log",
    wrapper:
        "v1.31.1/bio/cutadapt/se"

# Rule cutadapt3 removes final set of adapters from 5'end and discards reads based on trimmed read length
rule cutadapt3:
    input:
        fastq= WORKDIR + "/results/trimmed/{sample}-{unit}.2.1.fastq.gz",
    output:
        fastq= WORKDIR + "/results/trimmed/{sample}-{unit}.fastq.gz",
        qc= WORKDIR + "/results/trimmed/{sample}-{unit}.qc.txt",
    params:
        adapters=lambda w: str(
            '-m 20 -O 20 -g "'
            'r1adapter=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=20"'
            " --discard-trimmed"
        ),
    threads: config["params"]["cutadapt"]["threads"],
    log:
        WORKDIR + "/logs/cutadapt/{sample}-{unit}.log",
    wrapper:
        "v1.31.1/bio/cutadapt/se"

rule max_read_length:
    input:
        get_trimmed,
    output:
        WORKDIR + "/results/stats/max-read-length.json",
    log:
        WORKDIR + "/logs/max-read-length.log",
    conda:
        "../envs/pysam.yaml"
    script:
        "../scripts/get-max-read-length.py"

# =================================================================================================
#   Seqkit Stats of Trimmed Reads
# =================================================================================================

rule seqkit_stats:
    input:
        fastx= WORKDIR + "/results/trimmed/{sample}-{unit}.fastq.gz",
    output:
        stats= WORKDIR + "/results/stats/seqkit/{sample}-{unit}.stats.tsv",
    log:
        WORKDIR + "/logs/stats/{sample}-{unit}.seqkit_stats.log",
    params:
        extra="--all --tabular",
    threads: config["params"]["seqkit"]["threads"],
    wrapper:
        "v1.31.1-39-gb5b9878a/bio/seqkit/stats"

# =================================================================================================