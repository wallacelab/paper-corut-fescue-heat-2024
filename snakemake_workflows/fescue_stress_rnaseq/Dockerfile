FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="aee6eb98569d32744a70275062bd225f1331cc3f97a4ad0d30e4d4de32ac34d1"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: envs/kallisto.yaml
#   prefix: /conda-envs/b52ef4b47b1b3dd64c41bb3b8c581f55
#   name: kallisto
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - kallisto =0.48.0
RUN mkdir -p /conda-envs/b52ef4b47b1b3dd64c41bb3b8c581f55
COPY envs/kallisto.yaml /conda-envs/b52ef4b47b1b3dd64c41bb3b8c581f55/environment.yaml

# Conda environment:
#   source: envs/kraken.yaml
#   prefix: /conda-envs/aefba476cc36b155b9cee0e333ea284d
#   name: kraken2
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - kraken2=2.1.2
#     - krakentools=1.2
#     - libiconv=1.16
RUN mkdir -p /conda-envs/aefba476cc36b155b9cee0e333ea284d
COPY envs/kraken.yaml /conda-envs/aefba476cc36b155b9cee0e333ea284d/environment.yaml

# Conda environment:
#   source: envs/krona.yaml
#   prefix: /conda-envs/b6feeadd3ddf7f56165510dde54a3e3a
#   name: krona
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - krona=2.8.1
RUN mkdir -p /conda-envs/b6feeadd3ddf7f56165510dde54a3e3a
COPY envs/krona.yaml /conda-envs/b6feeadd3ddf7f56165510dde54a3e3a/environment.yaml

# Conda environment:
#   source: envs/matplotlib.yaml
#   prefix: /conda-envs/446f8633fb47ed71fe64ac01165c8869
#   name: plotting
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - matplotlib-base=3.2.2
#     - numpy=1.20.2
#     - python=3.7.0
#     - scipy=1.6.2
#     - seaborn=0.11.1
#     - seaborn-base=0.11.1
RUN mkdir -p /conda-envs/446f8633fb47ed71fe64ac01165c8869
COPY envs/matplotlib.yaml /conda-envs/446f8633fb47ed71fe64ac01165c8869/environment.yaml

# Conda environment:
#   source: envs/pandas.yaml
#   prefix: /conda-envs/3a10f4576e6a0dae688721cbab5194c2
#   name: pandas
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - pandas=1.5.2
#     - numpy=1.23.4
RUN mkdir -p /conda-envs/3a10f4576e6a0dae688721cbab5194c2
COPY envs/pandas.yaml /conda-envs/3a10f4576e6a0dae688721cbab5194c2/environment.yaml

# Conda environment:
#   source: envs/ribodetector.yaml
#   prefix: /conda-envs/8a08cae83873cfff945e9aa7090d6e05
#   name: ribodetector
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - ribodetector=0.2.7
#     - python=3.8.0
RUN mkdir -p /conda-envs/8a08cae83873cfff945e9aa7090d6e05
COPY envs/ribodetector.yaml /conda-envs/8a08cae83873cfff945e9aa7090d6e05/environment.yaml

# Conda environment:
#   source: envs/salmon.yaml
#   prefix: /conda-envs/efc608477657ed0b799eeaca3cb96d30
#   name: salmon
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - salmon=1.9.0
RUN mkdir -p /conda-envs/efc608477657ed0b799eeaca3cb96d30
COPY envs/salmon.yaml /conda-envs/efc608477657ed0b799eeaca3cb96d30/environment.yaml

# Conda environment:
#   source: envs/sleuth.yaml
#   prefix: /conda-envs/0ec91b423f78cc265ab50144e07c9393
#   name: sleuth
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - r-base=4.2.2
#     - r-sleuth=0.30.1
#     - bioconductor-rhdf5lib=1.20
#     - bioconductor-rhdf5=2.42
#     - r-gridextra=2.3
#     - r-pheatmap=1.0.12
#     - r-tidyverse=1.3
#     - r-ggpubr=0.4
#     - bioconductor-limma=3.54
#     - r-annotables=0.1.90
#     - openssl=1.1
RUN mkdir -p /conda-envs/0ec91b423f78cc265ab50144e07c9393
COPY envs/sleuth.yaml /conda-envs/0ec91b423f78cc265ab50144e07c9393/environment.yaml

# Conda environment:
#   source: envs/wasabi.yaml
#   prefix: /conda-envs/978ae1e41af0109313edea093d2ec6ec
#   name: wasabi
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - r-tidyverse=1.3
#     - r-wasabi=1.0.1
RUN mkdir -p /conda-envs/978ae1e41af0109313edea093d2ec6ec
COPY envs/wasabi.yaml /conda-envs/978ae1e41af0109313edea093d2ec6ec/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v1.20.0/bio/fastqc/environment.yaml
#   prefix: /conda-envs/8678da093d0504bc9ac36dfd48d4c62e
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - fastqc =0.11.9
RUN mkdir -p /conda-envs/8678da093d0504bc9ac36dfd48d4c62e
ADD https://github.com/snakemake/snakemake-wrappers/raw/v1.20.0/bio/fastqc/environment.yaml /conda-envs/8678da093d0504bc9ac36dfd48d4c62e/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v1.20.0/bio/multiqc/environment.yaml
#   prefix: /conda-envs/eebc5adeb5e51053b4346bc959a40041
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - multiqc =1.13
RUN mkdir -p /conda-envs/eebc5adeb5e51053b4346bc959a40041
ADD https://github.com/snakemake/snakemake-wrappers/raw/v1.20.0/bio/multiqc/environment.yaml /conda-envs/eebc5adeb5e51053b4346bc959a40041/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v1.20.0/bio/trim_galore/se/environment.yaml
#   prefix: /conda-envs/a225483da95f0f7f6981c658f7dd702b
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - trim-galore ==0.6.6
RUN mkdir -p /conda-envs/a225483da95f0f7f6981c658f7dd702b
ADD https://github.com/snakemake/snakemake-wrappers/raw/v1.20.0/bio/trim_galore/se/environment.yaml /conda-envs/a225483da95f0f7f6981c658f7dd702b/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/b52ef4b47b1b3dd64c41bb3b8c581f55 --file /conda-envs/b52ef4b47b1b3dd64c41bb3b8c581f55/environment.yaml && \
    mamba env create --prefix /conda-envs/aefba476cc36b155b9cee0e333ea284d --file /conda-envs/aefba476cc36b155b9cee0e333ea284d/environment.yaml && \
    mamba env create --prefix /conda-envs/b6feeadd3ddf7f56165510dde54a3e3a --file /conda-envs/b6feeadd3ddf7f56165510dde54a3e3a/environment.yaml && \
    mamba env create --prefix /conda-envs/446f8633fb47ed71fe64ac01165c8869 --file /conda-envs/446f8633fb47ed71fe64ac01165c8869/environment.yaml && \
    mamba env create --prefix /conda-envs/3a10f4576e6a0dae688721cbab5194c2 --file /conda-envs/3a10f4576e6a0dae688721cbab5194c2/environment.yaml && \
    mamba env create --prefix /conda-envs/8a08cae83873cfff945e9aa7090d6e05 --file /conda-envs/8a08cae83873cfff945e9aa7090d6e05/environment.yaml && \
    mamba env create --prefix /conda-envs/efc608477657ed0b799eeaca3cb96d30 --file /conda-envs/efc608477657ed0b799eeaca3cb96d30/environment.yaml && \
    mamba env create --prefix /conda-envs/0ec91b423f78cc265ab50144e07c9393 --file /conda-envs/0ec91b423f78cc265ab50144e07c9393/environment.yaml && \
    mamba env create --prefix /conda-envs/978ae1e41af0109313edea093d2ec6ec --file /conda-envs/978ae1e41af0109313edea093d2ec6ec/environment.yaml && \
    mamba env create --prefix /conda-envs/8678da093d0504bc9ac36dfd48d4c62e --file /conda-envs/8678da093d0504bc9ac36dfd48d4c62e/environment.yaml && \
    mamba env create --prefix /conda-envs/eebc5adeb5e51053b4346bc959a40041 --file /conda-envs/eebc5adeb5e51053b4346bc959a40041/environment.yaml && \
    mamba env create --prefix /conda-envs/a225483da95f0f7f6981c658f7dd702b --file /conda-envs/a225483da95f0f7f6981c658f7dd702b/environment.yaml && \
    mamba clean --all -y
