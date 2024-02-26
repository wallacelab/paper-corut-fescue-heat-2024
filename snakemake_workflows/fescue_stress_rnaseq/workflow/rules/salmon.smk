# =================================================================================================
#     Building Salmon Index
# =================================================================================================

if not config["settings"]["salmon"]["salmon_sa_index"]: 
    # Standard index
    rule salmon_index:
        input: 
            "resources/ref/transcriptome/plant.transcriptome.fasta",
        output:
            directory(WORKDIR + "/results/salmon/transcriptome_index")
        log:
            WORKDIR + "/logs/salmon/salmon_index.log"
        conda:
            "../envs/salmon.yaml"
        threads: 
            config["params"]["salmon_index"]["threads"]
        params:
            # optional parameters
            extra= config["params"]["salmon_index"]["index_extra"]
        message: "Building Salmon index....."
        shell:
            """
            salmon index --threads {threads} -t {input} -i {output} {params.extra} > {log} 2>&1
            """


if config["settings"]["salmon"]["salmon_sa_index"]: 
    # Decoy-aware index
    rule salmon_index_decoy:
        input:
            transcriptome= "resources/ref/transcriptome/plant.transcriptome.fasta",
            genome= "resources/ref/genome/plant.genome.fasta",
        output:
            gentrome= WORKDIR + "/results/salmon/decoy/gentrome.fa",
            decoy= WORKDIR + "/results/salmon/decoy/decoys.txt",
            bak= WORKDIR + "/results/salmon/decoy/decoys.txt.bak"
        conda:
            "../envs/salmon.yaml"
        threads:
            config["params"]["salmon_index"]["threads"]
        log:
            WORKDIR + "/logs/salmon/salmon_index_decoy.log"
        shell:
            """
            grep "^>" {input.genome} | cut -d " " -f 1 > {output.decoy}
            sed -i.bak -e 's/>//g' {output.decoy}
            cat {input.transcriptome} {input.genome} > {output.gentrome}
            """

    rule salmon_sa_index:
        input:
            gentrome= WORKDIR + "/results/salmon/decoy/gentrome.fa",
            decoy= WORKDIR + "/results/salmon/decoy/decoys.txt"
        output:
            directory(WORKDIR + "/results/salmon/gentrome_index")
        log:
            WORKDIR + "/logs/salmon/salmon_sa_index.log"
        conda:
            "../envs/salmon.yaml"
        threads:
            config["params"]["salmon_index"]["threads"]
        shell:
            """
            salmon index -p {threads} -t {input.gentrome} -d {input.decoy} -i {output} 2>&1 {log}
            """

# =================================================================================================
#     Quantifying Transcripts with Salmon
# =================================================================================================

if not config["settings"]["salmon"]["salmon_sa_index"]:
    rule salmon_quant:
        input:
            fq = get_salmon_input,
            index = WORKDIR + "/results/salmon/transcriptome_index"
        output:
            dir = directory(WORKDIR + "/results/salmon/{sample}-{unit}"),
            quant = WORKDIR + "/results/salmon/{sample}-{unit}/quant.sf",
            lib = WORKDIR + "/results/salmon/{sample}-{unit}/lib_format_counts.json",
            sam = WORKDIR + "/results/salmon/{sample}-{unit}/{sample}-{unit}.salmon_mapping.sam"
        log:
            WORKDIR + "/logs/salmon/salmon_quant/{sample}-{unit}.log"
        params:
            out_prefix= lambda w, output: os.path.dirname(output.quant),
            # optional parameters
            lib_type= config["params"]["salmon"]["lib_type"],
            extra= salmon_params
        conda:
            "../envs/salmon.yaml"
        threads: 
            config["params"]["salmon"]["threads"]
        message: "Quantifying transcripts with Salmon....."
        shell:
            """
            salmon quant -p {threads} -i {input.index} \
            -l {params.lib_type} -r {input.fq} \
            {params.extra} -o {params.out_prefix} \
            --writeMappings={output.sam} 2>&1 {log}
            """

if config["settings"]["salmon"]["salmon_sa_index"]:
    rule salmon_quant:
        input:
            fq = get_salmon_input,
            index = WORKDIR + "/results/salmon/gentrome_index"
        output:
            dir = directory(WORKDIR + "/results/salmon/{sample}-{unit}"),
            quant = WORKDIR + "/results/salmon/{sample}-{unit}/quant.sf",
            lib = WORKDIR + "/results/salmon/{sample}-{unit}/lib_format_counts.json"
        log:
            WORKDIR + "/logs/salmon/salmon_quant_sa/{sample}-{unit}.log"
        params:
            out_prefix= lambda w, output: os.path.dirname(output.quant),
            # optional parameters
            lib_type= config["params"]["salmon"]["lib_type"],
            extra= salmon_params
        conda:
            "../envs/salmon.yaml"
        threads: 
            config["params"]["salmon"]["threads"]
        message: "Quantifying transcripts with Salmon....."
        shell:
            """
            salmon quant -p {threads} -i {input.index} \
            -l {params.lib_type} {params.extra} \
            -r {input.fq} -o {params.out_prefix} 2>&1 {log}
            """

rule convert_bootsraps_to_tsv:
    input:
        WORKDIR + "/results/salmon/{sample}-{unit}",
    output:
        tsv= WORKDIR + "/results/salmon/{sample}-{unit}/quant_bootstraps.tsv"
    params:
        out_prefix= lambda w, output: os.path.dirname(output.tsv),
    conda:
        "../envs/salmon.yaml"
    shell:
        """
        python scripts/convert_bootstraps_to_tsv.py {input} {params.out_prefix}
        """

# =================================================================================================