#!/usr/bin/env python3
"""
Title: run_enrichr.py
Date: 2021-08-08
Author: Erno Hänninen
Description:
   This script runs Palantir trajectory inference algorithm for input data.
   The palantir outputs gene expression trends for each lineage. 
   This script runs spearman correlation to test whether the gene expression trend correlates with lineage probability.
   Spearman correlation outputs correlation table, which is ordered so that highest correlation is first.
   The script prepares the input for GSEA and Enrichr.
   For GSEA it creates a ranked list with gene name and correlation metric, gene with the highest correlation is first.
   Significant and positively correlated genes are used as input for Enrichr.
   This script returns the, plots, computed diffusion components, correlation tables, ranked lists and gene lists for each lineage

List of functions:
    -    
List of "non standard" modules:
    -    scanpy, pandas, gseapy, scipy, numpy, matplotlib, seaborn, palantir
Error handling:
    -
Procedure:
    - This script is called from processes.nf
"""

import argparse
import palantir
import scanpy as sc
import numpy as np
import pandas as pd
import os
import scipy 
from gseapy.plot import gseaplot
import gseapy as gsea
# Plotting 
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
from scipy.sparse import csr_matrix
# Inline plotting
#%matplotlib inline

sns.set_style('ticks')
matplotlib.rcParams['figure.figsize'] = [4, 4]
matplotlib.rcParams['figure.dpi'] = 100
matplotlib.rcParams['image.cmap'] = 'Spectral_r'
warnings.filterwarnings(action="ignore", module="matplotlib", message="findfont")

# Set random seed
np.random.seed(10)

#initialize parser and get the input arguments
parser = argparse.ArgumentParser(description = "")
parser.add_argument("input_data", help = "input data")
parser.add_argument("dataset", help = "dataset name")
parser.add_argument("start_cell", help = "dataset name")
parser.add_argument("palantir_clusters", help = "dataset name")


#parser.add_argument("ti_method", help = "Trajectory method")
args = parser.parse_args()
#Read the arguments to variables
input_data, dataset, start_cell, clusters = args.input_data, args.dataset, args.start_cell, args.palantir_clusters

#Import the data
adata = sc.read(input_data)

n_diffusion_components = 10

#Run diffusion map
pca_projections = pd.DataFrame(adata.obsm['X_pca'], index=adata.obs_names)
dm_res = palantir.utils.run_diffusion_maps(pca_projections, n_components=n_diffusion_components)

#Estimating low dimensional embeddeing of the data
ms_data = palantir.utils.determine_multiscale_space(dm_res)

##Plot umap color by sample_id
#sc.pl.embedding(adata, basis='umap', color="sample_id", size = 1000000 / adata.n_obs, save="_sample_id.png")

#Plot umap color by new_leiden_1_harmony
sc.pl.embedding(adata, basis='umap', color=clusters, size = 1000000 / adata.n_obs, save="_" + clusters + ".png")

#Convert adata.X to sparse matrix, this is needed for magic imputations
adata.X = csr_matrix(adata.X)
#Run magic imputations
adata.layers['MAGIC_imputed_data'] = palantir.utils.run_magic_imputation(adata, dm_res)

#Visualize diffusion components
umap = pd.DataFrame(adata.obsm['X_umap'], index=adata.obs_names)
palantir.plot.plot_diffusion_components(umap, dm_res)
plt.savefig("figures/_diffusion_components.png")

#Write diffusion components to file
print(dm_res["EigenVectors"])
dm_res["EigenVectors"].to_csv(dataset+"_diff_components.csv", sep = "\t")

#Running palantir
pr_res = palantir.core.run_palantir(ms_data, start_cell, num_waypoints=1000)

#Start cell
palantir.plot.highlight_cells_on_tsne(umap, start_cell)
plt.savefig("figures/start_cell.png")


#Terminal states
palantir.plot.highlight_cells_on_tsne(umap, pr_res.branch_probs.columns.values)
plt.savefig("figures/terminal_states.png")

#Rename the terminal states
terminal_states = pr_res.branch_probs.columns.values.tolist() #Extract the terminal states to a list
for i in range(len(pr_res.branch_probs.columns.values)):
    pr_res.branch_probs.rename(columns = {pr_res.branch_probs.columns.values[i] : i+1}, inplace = True)
    
print(pr_res.branch_probs.columns.values)

#Visualizing the results
#First plot is the pseudotime, second differentiation potential and
#third is the terminal states probabilities
palantir.plot.plot_palantir_results(pr_res, umap)
plt.savefig("figures/pseudotime.png")

#Get highly variable genes from the dataset 
gene_list = adata.var.index
print("GENE LIST for gene trends")
print(len(gene_list))

#Computing gene trends
#The function returns a dictionary of gene expression trends and standard deviations for each branch
#Using this function with variable genes extracted above
imp_df = pd.DataFrame(adata[:, gene_list].layers['MAGIC_imputed_data'],
                     index=adata.obs_names, columns=gene_list)

gene_trends = palantir.presults.compute_gene_trends( pr_res, imp_df.loc[:, gene_list])

#Run spearman correlation to test if gene trend correlates with lineage probability, this way we identified the genes associated with each lineage 
correlation_dict = {} #Dictinary where the terminal state is stored as a key and correlation table as a value
terminal_states = pr_res.branch_probs.columns.values.tolist() #Extract the terminal states to a list
for terminal_state in terminal_states: #Loop thru the terminal states
    print(terminal_state)
    res_list = [] #List where the spearman results are stored
    lineage_probability = gene_trends[terminal_state]["trends"].columns.values.tolist() #Get the lineage probability from gene_trends table
    for i in range(gene_trends[terminal_state]["trends"].shape[0]): #Loop over the gene_trends table
        gene_trend = gene_trends[terminal_state]["trends"].iloc[i, :].values.tolist() #Extract one row at a time
        gene = gene_trends[terminal_state]["trends"].index.values[i] #Get the row index (gene name)
        correlation, pvalue = scipy.stats.spearmanr(gene_trend, lineage_probability) #Run the spearman correlation
        res_list.append([gene,correlation,pvalue]) #Store the results to a list
    
    df = pd.DataFrame(res_list, columns=["gene", 'correlation','pvalue'])#Store the results from the list to a dataframe
    df = df.set_index('gene') #Change the gene name as index
    df.index.names = [None] #Remove the index column name
    df = df.sort_values(by = ["correlation"], ascending=False) #Sort values by correlation  
    #df = df.loc[df['pvalue'] <= 0.05] #Get only the significant genes
    print("Correlation table")
    print(df.head())
    df.to_csv("correlation_"+ dataset + "_lineage_" + str(terminal_state) + ".csv", sep="\t") #Write correlation results to csv

    ranked_list = df.iloc[:,0]
    ranked_list.to_csv(dataset + "_ranked_lineage" + str(terminal_state)+ ".rnk", sep = " ", header=False)
    print("Ranked list")
    print(ranked_list)


    gene_list = df.loc[df["correlation"] > 0]
    gene_list = gene_list.loc[gene_list["pvalue"] <= 0.05]
    gene_list = gene_list.index.values.tolist()
    print("Gene list")
    print(len(gene_list))
    print(gene_list[0:10])

    with open(dataset + "_geneList_lineage" + str(terminal_state) + ".txt", "w") as file:
        for gene in gene_list:
            file.write(gene + "\n")
    
    print("\n")

