import csv
import os, argparse, glob, shutil
import pandas as pd
import numpy as np
import sys

# logging
with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f

    kraken_report = pd.read_csv(snakemake.input[0], sep="\t")

    kraken_report_edited = kraken_report.drop([0,1,2,3,4,5,6])
    kraken_report_edited["#Classification"] = kraken_report_edited["#Classification"].apply(lambda x: x.split('|')[-1])
    kraken_report_edited.rename(columns={"#Classification": "Sample_name"},inplace=True)

    kraken_report_edited_transposed = kraken_report_edited.T
    kraken_report_edited_transposed.reset_index(level=0, inplace=True)
    new_header = kraken_report_edited_transposed.iloc[0] #grab the first row for the header
    kraken_report_edited_transposed = kraken_report_edited_transposed[1:] #take the data less the header row
    kraken_report_edited_transposed.columns = new_header #set the header row as the df header
    kraken_report_edited_transposed["Sample_name"] = kraken_report_edited_transposed["Sample_name"].apply(lambda x: x.split('.kraken.report')[0])
    kraken_report_edited_transposed.rename(columns={"s__Epichloe_festucae": "Epichloe_festucae_counts"},inplace=True)

    kraken_report_edited_transposed.to_csv(snakemake.output[0], index=False, sep="\t")