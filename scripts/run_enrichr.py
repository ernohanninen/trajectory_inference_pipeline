#!/usr/bin/env python3
"""
Title: run_enrichr.py
Date: 2021-08-05
Author: Erno Hänninen
Description:
   This script runs Enrichr to a list of genes.
   The genes are the top ranked and significant genes from tradeSeq associationTest / spearman correlation
   Output is Enrichr plots
List of functions:
    -    
List of "non standard" modules:
    -    scanpy, pandas, gseapy
Error handling:
    -
Procedure:
    - This script is called from processes.nf
"""
#Import packages
import scanpy as sc
import os
import pandas as pd
import re
from gseapy.plot import gseaplot
import gseapy as gsea
import argparse

#Initialize parser and read the arguments
parser = argparse.ArgumentParser(description = "")

parser.add_argument("gene_list", help = "gene list")
parser.add_argument("input_data", help = "input data")
parser.add_argument("dataset", help = "dataset name")
parser.add_argument("gene_set_path", help = "path to gene set")
parser.add_argument("ti_method", help = "trajectory method")
parser.add_argument("num_genes", help = "number of genes to run the enrichr with")

args = parser.parse_args()
#Read the arguments to variables
gene_list_path, input_data, dataset, gene_set_path, ti_method, num_genes = args.gene_list, args.input_data, args.dataset, args.gene_set_path, args.ti_method, args.num_genes

#Read the h5ad file
adata = sc.read(input_data)

#Get the trajectory name from the file path
trajectory_name = gene_list_path.split("_")[-1].replace(".txt", "").replace("_ranked", "")
print("TRAJECTORY : ", trajectory_name)

#Read the genes (adata.var.index) from file to list
gene_list = []
with open(gene_list_path, "r") as file:
    for line in file:
        line = line.replace("\n", "")
        gene_list.append(line)

#Get the user specified number of top genes for enrichr
gene_list = gene_list[0:int(num_genes)]
print(len(gene_list))
print(gene_list[0:10])

#Specify the output dir
out_dir = dataset + "_enrichr_" + trajectory_name

#Using all genes from adata.raw as background genes
background_set = adata.raw.var.index.to_list()

#Run the enrichr
enr = gsea.enrichr(gene_list=gene_list,
                gene_sets=gene_set_path,
                background = background_set,
                outdir = out_dir)
