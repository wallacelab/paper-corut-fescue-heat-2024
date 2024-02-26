# =======================================================================================================
#    Creating a symbolic link for the reference transcriptome fasta
# =======================================================================================================

rule plant_transcriptome_symlink:
    input:
        transcriptome = config["ref"]["plant"]["transcriptome"]
    output:
        "resources/ref/transcriptome/plant.transcriptome.fasta",
    message:
        "Creating symbolic link for reference transcriptome fasta..."
    threads: 1
    shell:
        """
        ln -s {input.transcriptome} {output}
        """

# =======================================================================================================
#    Creating a symbolic link for the reference genome fasta
# =======================================================================================================

rule plant_genome_symlink:
    input:
        genome = config["ref"]["plant"]["genome"]
    output:
        "resources/ref/genome/plant.genome.fasta",
    message:
        "Creating symbolic link for reference genome fasta..."
    threads: 1
    shell:
        """
        ln -s {input.genome} {output}
        """

# =======================================================================================================
#    Creating a symbolic link for the reference transcriptome fasta
# =======================================================================================================

rule endophyte_genome_symlink:
    input:
        genome = config["ref"]["endophyte"]["genome"]
    output:
        "resources/ref/genome/endophyte.genome.fasta",
    message:
        "Creating symbolic link for reference genome fasta..."
    threads: 1
    shell:
        """
        ln -s {input.genome} {output}
        """

# =======================================================================================================
#    Creating a symbolic link for the reference transcriptome fasta
# =======================================================================================================

rule chloroplast_genome_symlink:
    input:
        chloroplast = config["ref"]["chloroplast"]["genome"]
    output:
        "resources/ref/chloroplast/chloroplast.genome.fasta",
    message:
        "Creating symbolic link for chloroplast genome fasta..."
    threads: 1
    shell:
        """
        ln -s {input.chloroplast} {output}
        """

# =======================================================================================================
#    Parsing tall fescue transcripts information(annotation) file
# =======================================================================================================

rule parse_annotation:
    input:
        tcp_anno = config["ref"]["plant"]["annotation"]
    output:
        t2g = "resources/ref/annotation/fescue.transcriptome_info.txt",
        go_list = "resources/ref/annotation/fescue.go_list.txt",
    message:
        "Parsing transcriptome annotation file..."
    log:
        WORKDIR + "/logs/ref/parse_annotation.log"
    script:
        "../scripts/parse_annotation_file.py"
    
# =======================================================================================================

