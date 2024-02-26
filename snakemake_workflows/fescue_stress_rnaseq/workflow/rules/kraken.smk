# =================================================================================================
#     Building Kraken Database
# =================================================================================================

rule kraken_build:
    input:
        fasta= "resources/ref/genome/endophyte.genome.fasta",
    output:
        directory(WORKDIR + "/results/kraken/kraken_db"),
    conda:
        "../envs/kraken.yaml"
    threads: 
        config["params"]["kraken2"]["threads"]
    log:
        log1 = WORKDIR + "/logs/kraken/kraken_build/download-taxonomy.log",
        log2 = WORKDIR + "/logs/kraken/kraken_build/add-to-library.log",
        log3 = WORKDIR + "/logs/kraken/kraken_build/build.log"
    shell:
        """
        kraken2-build --threads {threads} --download-taxonomy --db {output} 2> {log.log1}
        kraken2-build --threads {threads} --add-to-library {input.fasta} --db {output} 2> {log.log2}
        kraken2-build --threads {threads} --build --db {output} 2> {log.log3}
        """

# =================================================================================================
#     Running Kraken2 with the Specified Database
# =================================================================================================

rule run_kraken:
    input:
        reads= WORKDIR + "/results/trimmed/{sample}-{unit}_1_trimmed.fq.gz",
        db= WORKDIR + "/results/kraken/kraken_db"
    output:
        kraken_report = WORKDIR + "/results/kraken/{sample}-{unit}.kraken.report",
        kraken_out = temp(WORKDIR + "/results/kraken/{sample}-{unit}.kraken")
    params:
        extra = config["params"]["kraken2"]["extra"]
    conda:
        "../envs/kraken.yaml"
    threads: 
        config["params"]["kraken2"]["threads"]
    log:
        WORKDIR + "/logs/kraken/run_kraken/{sample}-{unit}.kraken.log"
    shell:
        """
        kraken2 --db {input.db} --threads {threads} --output {output.kraken_out} \
        --report {output.kraken_report} --use-names {params.extra} {input.reads} >> {log}
        """

# =================================================================================================
#     Combining Multiple Kraken Reports into a Combined Report File
# =================================================================================================

rule kreport_to_mpa:
    input:
        WORKDIR + "/results/kraken/{sample}-{unit}.kraken.report",
    output:
        WORKDIR + "/results/kraken/mpa_report/{sample}-{unit}.mpa.report.txt"
    conda:
        "../envs/kraken.yaml"
    log:
        WORKDIR + "/logs/kraken/kreport_to_mpa/{sample}-{unit}.kreport_to_mpa.log"
    shell:
        """
        kreport2mpa.py --display-header --no-intermediate-ranks -r {input} -o {output} 2> {log}
        """

rule combine_mpa:
    input:
        expand(
            WORKDIR + "/results/kraken/mpa_report/{u.sample_name}-{u.unit_name}.mpa.report.txt",
            u=units.itertuples(),
        )
    output:
        WORKDIR + "/results/kraken/mpa_report/combined.mpa.report"
    conda:
        "../envs/kraken.yaml"
    log:
        WORKDIR + "/logs/kraken/combine_mpa/combine_mpa.log"
    shell:
        """
        combine_mpa.py -i {input} -o {output}  2> {log}
        """

rule parse_mpa_report:
    input:
        WORKDIR + "/results/kraken/mpa_report/combined.mpa.report"
    output:
        WORKDIR + "/results/kraken/final_kraken_report.txt"
    conda:
        "../envs/pandas.yaml"
    log:
        WORKDIR + "/logs/kraken/parse_mpa_report/parse_mpa_report.log"
    script:
        "../scripts/parse_kraken_mpa_report.py"

rule plot_final_kraken_report:
    input:
        WORKDIR + "/results/kraken/final_kraken_report.txt",
        config["samples"]
    output:
        plot_by_condition = report(
            WORKDIR + "/results/plots/kraken/final_kraken_report.plot_by_condition.pdf",
            caption="../report/plot_final_kraken_report.rst",
            category="Kraken2 Report",
        ),
        # plot_by_genotype = report(
        #     WORKDIR + "/results/plots/kraken/final_kraken_report.plot_by_genotype.pdf",
        #     caption="../report/plot_final_kraken_report.rst",
        #     category="Kraken2 Report",
        # )
    conda:
        "../envs/matplotlib.yaml"
    log:
        WORKDIR + "/logs/plots/kraken/plot_final_kraken_report.log"
    script:
        "../scripts/plot_kraken_final_report.py"

# =================================================================================================
#     Running Kraken2 with a Custom Database
#     ## Specify the Database Name and Path in the Config File  
# =================================================================================================

rule run_kraken_custom:
    input:
        reads= WORKDIR + "/results/trimmed/{sample}-{unit}_1_trimmed.fq.gz",
        db= config["params"]["kraken2"]["db_path"]
    output:
        kraken_report = WORKDIR + "/results/kraken_{db_name}/{sample}-{unit}.kraken.report",
        kraken_out = temp(WORKDIR + "/results/kraken_{db_name}/{sample}-{unit}.kraken")
    params:
        extra = config["params"]["kraken2"]["extra"]
    conda:
        "../envs/kraken.yaml"
    threads: 
        config["params"]["kraken2"]["threads"]
    log:
        WORKDIR + "/logs/kraken/run_kraken_custom/kraken_{db_name}/{sample}-{unit}.kraken.log"
    shell:
        """
        kraken2 --db {input.db} --threads {threads} --quick --output {output.kraken_out} \
        --report {output.kraken_report} --use-names {params.extra} {input.reads} >> {log}
        """

rule krona:
    input:
        expand(WORKDIR + "/results/kraken_{db_name}/{u.sample_name}-{u.unit_name}.kraken.report",
            u=units.itertuples(),
            db_name=config["params"]["kraken2"]["db_name"]
        )
    output:
        WORKDIR + "/results/kraken_{db_name}/krona/multi-krona.html"
    conda:
        "../envs/krona.yaml"
    log:
        WORKDIR + "/logs/kraken/krona/kraken_{db_name}/multi-krona.log"
    shell:
        """
        ktUpdateTaxonomy.sh
        ktImportTaxonomy -t 5 -m 3 -o {output} {input} 2> {log}
        """
        
# =================================================================================================