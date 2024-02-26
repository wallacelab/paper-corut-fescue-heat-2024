# Write me a snakemake rule for catchSalmon.R

rule catch_salmon:
    input:
        salmon_quant = expand(WORKDIR + "/results/salmon/{u.sample_name}-{u.unit_name}", 
                               u=units.itertuples()),
        sample_metadata = WORKDIR + "/results/tables/sample_metadata.tsv",
    output:
        catch = WORKDIR + "/results/tables/catchSalmon/catchSalmon.RDS",
    log:
        log = WORKDIR + "/logs/catch_salmon/catchSalmon.log",
    conda:
        "../envs/edger.yaml"
    script:
        "../scripts/catchSalmon.R"