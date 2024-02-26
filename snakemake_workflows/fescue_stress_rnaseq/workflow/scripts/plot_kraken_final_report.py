import csv
import sys
import os, argparse, glob, shutil
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy as sp
from scipy import stats
import matplotlib.ticker as ticker
from matplotlib.pyplot import hist

# logging
with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f

    kraken_final_report = pd.read_csv(snakemake.input[0], 
                                    sep="\t")

    units = pd.read_csv(snakemake.input[1], 
                                    sep="\t")

    kraken_final_report_edited = pd.concat([kraken_final_report, units], axis=1)
    print(kraken_final_report_edited)
    
    kraken_final_report_edited_sort = kraken_final_report_edited.sort_values(by=['treatment', 'genotype'])

    fig_dims = (16, 8)
    fig, ax = plt.subplots(figsize=fig_dims)
    sns_plot_by_cond = sns.barplot(x ='Sample_name', y = 'Epichloe_festucae_counts', hue="treatment", 
                                   data = kraken_final_report_edited_sort, 
                                   dodge=False)
    
    sns_plot_by_cond.set_xticklabels(sns_plot_by_cond.get_xticklabels(), 
                                     rotation=40, ha="right", fontsize=6)
    
    sns_plot_by_cond.set_xlabel("Sample Name", fontsize=18, weight='bold')
    
    sns_plot_by_cond.set_ylabel("Epichloe festucae Counts (Kraken2)", 
                                fontsize=18, weight='bold')
    
    sns_plot_by_cond.legend(loc='upper right', fontsize=16, 
                            title = "Treatment", title_fontsize= 20)
    fig.tight_layout()
    
    ## Save the plot
    sns_plot_by_cond.figure.savefig(snakemake.output.plot_by_condition, dpi=600)

    # fig_dims = (16, 8)
    # fig, ax = plt.subplots(figsize=fig_dims)
    # sns_plot_by_tiss = sns.barplot(x ='Sample_name', y = 'Epichloe_festucae_counts', hue="genotype", data = kraken_final_report_edited)
    # sns_plot_by_tiss.set_xticklabels(sns_plot_by_tiss.get_xticklabels(), rotation=40, ha="right", fontsize=7)
    # sns_plot_by_cond.set_xlabel("Sample Name",fontsize=18, weight='bold')
    # sns_plot_by_cond.set_ylabel("Epichloe festucae Counts (Kraken2)",fontsize=18, weight='bold')
    # sns_plot_by_tiss.legend(loc='center right', fontsize=16, title = "Genotype", title_fontsize= 20)
    # fig.tight_layout()
    
    # ## Save the plot
    # sns_plot_by_tiss.figure.savefig(snakemake.output.plot_by_tissue, dpi=600)