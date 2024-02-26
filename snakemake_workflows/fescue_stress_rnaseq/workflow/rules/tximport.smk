# =================================================================================================
#       Run tximport
# =================================================================================================

rule tximport:
    input:
        quant = expand(WORKDIR + "/results/salmon/{u.sample_name}-{u.unit_name}/quant.sf", 
                        u=units.itertuples()),
        metadata = WORKDIR + "/results/tables/sample_metadata.tsv",
    output:
        txi = WORKDIR + "/results/tximport/txi.RDS"
    conda:
        "../envs/tximport.yaml"
    threads: 
        config["params"]["tximport"]["threads"]
    log:
        WORKDIR + "/logs/tximport/tximport_salmon.log"
    script:
        "../scripts/tximport.R"

# =================================================================================================