import csv
import sys
import os, argparse, glob, shutil
import pandas as pd
import numpy as np

# logging
sys.stderr = open(snakemake.log[0], "w")

# read in tall fescue annotation file
tall_fescue_anno = pd.read_csv(snakemake.input.tcp_anno, 
                               sep=",", 
                               usecols=["TFQuery","TAIR", "BP GO", "MF GO"], 
                               dtype=object, 
                               na_values=[], 
                               keep_default_na=False)

# rename columns
tall_fescue_anno.rename(columns={   
    'TFQuery': 'target_id', 
    'TAIR': 'TAIR_hit',
    'BP GO': 'bp_GO',
    'MF GO': 'mf_GO',},
                        inplace=True)

tall_fescue_anno["bp_GO"]= tall_fescue_anno["bp_GO"].apply(lambda x: x.split(' ')) # split GO terms
tall_fescue_anno["mf_GO"]= tall_fescue_anno["mf_GO"].apply(lambda x: x.split(' ')) # split GO terms

tall_fescue_bp_GO_list = tall_fescue_anno[["target_id", "bp_GO"]] # create a new dataframe with only GO terms
tall_fescue_bp_GO_list = tall_fescue_bp_GO_list.explode('bp_GO') # explode the GO terms
tall_fescue_bp_GO_list = tall_fescue_bp_GO_list.replace('NA', np.NaN) # replace NA with NaN
tall_fescue_bp_GO_list.dropna(inplace=True) # drop NaN
tall_fescue_bp_GO_list.rename(columns={
    'target_id': 'GeneID', 
    "bp_GO": 'GOID'}, 
                              inplace=True) # rename columns

tall_fescue_bp_GO_list = tall_fescue_bp_GO_list[['GOID','GeneID']] # reorder columns

tall_fescue_mf_GO_list = tall_fescue_anno[["target_id", "mf_GO"]] # create a new dataframe with only GO terms
tall_fescue_mf_GO_list = tall_fescue_mf_GO_list.explode('mf_GO') # explode the GO terms
tall_fescue_mf_GO_list = tall_fescue_mf_GO_list.replace('NA',np.NaN) # replace NA with NaN
tall_fescue_mf_GO_list.dropna(inplace=True) # drop NaN
tall_fescue_mf_GO_list.rename(columns={
    'target_id': 'GeneID', 
    "mf_GO": 'GOID'}, 
                              inplace=True) # rename columns

tall_fescue_mf_GO_list = tall_fescue_mf_GO_list[['GOID','GeneID']] # reorder columns

tall_fescue_anno = tall_fescue_anno.drop(['bp_GO', 'mf_GO'], 
                                         axis=1) # drop GO terms from tall_fescue_anno


tall_fescue_anno.to_csv(snakemake.output.t2g, 
                        index=False, 
                        sep="\t") # write tall fescue annotation file

tall_fescue_go_list = pd.concat([tall_fescue_bp_GO_list,tall_fescue_mf_GO_list],
                                ignore_index=True) # combine bp and mf GO terms
tall_fescue_go_list.sort_values(['GeneID'], 
                                inplace=True) # sort by GeneID
tall_fescue_go_list.to_csv(snakemake.output.go_list, 
                           index=False, 
                           sep="\t") # write tall fescue GO list
# tall_fescue_bp_GO_list.to_csv(snakemake.output.bp_GO, 
#                               index=False, 
#                               sep="\t") # write tall fescue bp GO list
# tall_fescue_mf_GO_list.to_csv(snakemake.output.mf_GO, 
#                               index=False, 
#                               sep="\t") # write tall fescue mf GO list