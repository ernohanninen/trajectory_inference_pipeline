#!/usr/bin/env python3
"""
Title: run_paga.py
Date: 2021-08-09
Author: Erno Hänninen
Description:
   This scripts runs PAGA + diffusion pseudotime trajectory inference
   This is a script where the parameters (clusters, start cluster) are configured for vec_smc_per dataset
List of functions:
    -    
List of "non standard" modules:
    -    scanpy, pandas, numpy, fa2, matplotlib
Error handling:
    -
Procedure:
    - This script is called from processes.nf
"""
#Load python packages
import numpy as np
import pandas as pd
import argparse
import scanpy as sc
import fa2
import matplotlib.pyplot as pl
from matplotlib import rcParams
import matplotlib
#Plotting settings
matplotlib.use('Agg')
sc.settings.set_figure_params(dpi=80, frameon=False, figsize=(5, 5), facecolor='white')  # low dpi (dots per inch) yields small inline figures

sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)

#Initialize the parser and read the input arguments
parser = argparse.ArgumentParser(description = "")
parser.add_argument("input_data", help = "input data")
parser.add_argument("dataset", help = "dataset name")
args = parser.parse_args()
#Get the arguments to variables
input_data, dataset = args.input_data, args.dataset


#Load the data
adata = sc.read(input_data)

#Plot umap
sc.pl.umap(adata, color = "new_leiden_1_harmony", legend_loc = "on data",  save="_new_leiden_1_harmony.png")
#Run paga
sc.tl.paga(adata, groups='new_leiden_1_harmony')
sc.pl.paga(adata, color=['new_leiden_1_harmony'], save="_graph.png") #Plot paga
sc.tl.draw_graph(adata, init_pos='paga') #Draw force directed graph
sc.pl.draw_graph(adata, color=['new_leiden_1_harmony'], legend_loc='on data', save="_1.png") #Plot force directed graph
adata.uns['iroot'] = np.flatnonzero(adata.obs['new_leiden_1_harmony']  == 6)[0] #Set start cluster
sc.tl.dpt(adata) #Run diffusion pseudotime function
sc.pl.draw_graph(adata, color=['new_leiden_1_harmony', 'dpt_pseudotime'], legend_loc='on data', save="_2.png") #plot force directed graph
sc.tl.umap(adata, init_pos='paga') #run umap using paga
sc.pl.umap(adata, color=['new_leiden_1_harmony', 'dpt_pseudotime'], legend_loc='on data', save="_pseudotime.png") #Plot umap