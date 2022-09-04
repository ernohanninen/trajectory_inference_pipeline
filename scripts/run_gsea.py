#!/usr/bin/env python3
"""
Title: run_gsea.py
Date: 2021-08-08
Author: Erno Hänninen
Description:
   This script runs GSEA to a ranked list.
   The list from slingshot + tradeseq is ranked using values from Wald test (highest first)
   The list from palantir + spearman correlation is ranked using spearman correlation (highest value first)
   Output is GSEA plots
List of functions:
    -    
List of "non standard" modules:
    -    scanpy, pandas, gseapy
Error handling:
    -
Procedure:
    - This script is called from processes.nf
"""

import scanpy as sc
import argparse
import os
import pandas as pd
import re
from gseapy.plot import gseaplot
import gseapy as gsea

#Creating the parser and reading the arguments
parser = argparse.ArgumentParser(description = "")
parser.add_argument("ranked_list", help = "ranked list")
parser.add_argument("dataset", help = "dataset name")
parser.add_argument("gene_set_path", help = "path to gene set")
parser.add_argument("ti_method", help = "Trajectory method")
args = parser.parse_args()
#Reading the arguments from parser to variables
ranked_list, dataset, gene_set_path, ti_method = args.ranked_list, args.dataset, args.gene_set_path, args.ti_method

#Extract the trajectory name from the ranked_list path name
trajectory_name = ranked_list.split("_")[-1].replace(".rnk", "").replace("_ranked", "")
print("TRAJECTORY : ", trajectory_name)

#Read the ranked list from the file 
ranked_list = pd.read_csv(ranked_list, sep = " ", index_col=0, names=["#ranked"])
ranked_list = ranked_list.iloc[:,0] #Convert the table to right format

print("PRINT RANKED LIST")
print(ranked_list.head())

#Specify the output dir
out_dir = dataset + "_gsea_" + trajectory_name


#Run the gsea
gsea_res = gsea.prerank(rnk=ranked_list, gene_sets=gene_set_path, permutation_num=1000, outdir = out_dir, format="png", seed=6, graph_num=6, min_size=15)
#Plotting
gsea_res.res2d.sort_index().head()
terms = gsea_res.res2d.index
gseaplot(rank_metric=gsea_res.ranking, term=terms[0], **gsea_res.results[terms[0]])

print("READY")