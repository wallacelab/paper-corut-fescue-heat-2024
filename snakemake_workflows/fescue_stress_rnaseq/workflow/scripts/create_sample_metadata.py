# Purpose: Create sample metadata file for downstream analysis
import sys
import csv
import os, argparse, glob, shutil
import pandas as pd
import numpy as np

sys.stderr = open(snakemake.log[0], "w") # Redirect stderr to log file

sample_metadata = snakemake.params.units[["sample_name", "unit_name"]].merge(snakemake.params.samples, on="sample_name") # Merge sample and unit metadata

sample_metadata["sample_name"] = sample_metadata.apply(
    lambda row: "{}-{}".format(row["sample_name"], row["unit_name"]), axis=1
)   # Create sample name column

# Create path column for salmon results
sample_metadata["path"] = list(snakemake.input.wasabi_output)

# Delete columns that are not needed
del sample_metadata["unit_name"]

# Rename the sample names column
sample_metadata.rename(columns={'sample_name': 'sample'}, inplace=True)

sample_metadata.to_csv(snakemake.output[0], sep="\t", index=False)  # Write the sample metadata file