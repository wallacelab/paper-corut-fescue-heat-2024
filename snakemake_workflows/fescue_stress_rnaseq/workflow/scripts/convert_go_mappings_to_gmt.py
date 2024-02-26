# Purpose: Convert GO mappings to GMT file for GSEA analysis

import sys
import csv
import os, argparse, glob, shutil
import pandas as pd

sys.stderr = sys.stdout = open(snakemake.log[0], "w") # Redirect stderr to log file

# Read annotation file into dataframe with only the columns we need
df = pd.read_csv(snakemake.input["tcp_hits"], sep=",", usecols=["TFQuery", "BP descriptionMF"])

# Rename columns and switch order
df.columns = ["GeneID", "GOID"]

# Switch order of columns
df = df[["GOID", "GeneID"]]

# Remove rows if GOID is NaN
df = df.dropna(subset=['GOID'])

print(df.head())

# Group by GO ID
grouped = df.groupby('GOID')['GeneID'].apply(list)

print(grouped.head())

# Write to GMT file
with open(snakemake.output["gmt"], 'w') as f:
    for go_id, genes in grouped.items():
        f.write(f"{go_id}\t{go_id}\t" + "\t".join(genes) + "\n")
