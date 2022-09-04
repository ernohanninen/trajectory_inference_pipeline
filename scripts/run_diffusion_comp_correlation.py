"""
- This script run's spearman correlation to identify the genes which associate to specific component and runs GSEA and ENRICHR to the genes
"""
from gseapy.plot import gseaplot
import gseapy as gsea
import scanpy as sc
import numpy as np
import pandas as pd
import os
import scipy
# Plotting
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
import argparse

#initialize parser and get the input arguments
parser = argparse.ArgumentParser(description = "")

parser.add_argument("input_data", help = "input data")
parser.add_argument("dataset", help = "dataset name")
parser.add_argument("gene_set", help = "Gene set for GSEA and Enrichr")
parser.add_argument("num_genes", help = "Number of genes")
parser.add_argument("diffusion_components", help = "diffusion components")


#parser.add_argument("ti_method", help = "Trajectory method")
args = parser.parse_args()

#Read the arguments to variables
input_data, dataset, gene_set, num_genes,diffusion_components = args.input_data, args.dataset, args.gene_set, args.num_genes, args.diffusion_components

#Read the input data to adata object
adata = sc.read(input_data)

components = pd.read_csv(diffusion_components, sep="\t", index_col=0)
print(components.head())

genes = adata.raw.var.gene_ids.index.values.tolist() #Get the gene names
correlation_dict = {} #Dictionary where the terminal state is stored as a key and correlation table $

#The loop loops as many times as shape[1] of dmres["eigenvectors"]
#Hence the gene should always correspond the row
for i in range(components.shape[1]): #Loop over the dm_res components
    print("COMPONENT : ", i)
    counter = 0 #Gene index
    res_list = [] #List where the spearman results are stored
    #Prepare the input for spearman correlation
    component_values = components.iloc[:,i].tolist()  #Get the dm_res component values to a list

    for j in range(adata.raw.X.shape[1]): #Loop over expression matrix
        #Prepare input for separman correlation
        expression_matrix = adata.raw.X.getcol(j).transpose() #Extracts column from adata.raw.X expression matrix by index and transposes it
        expression_list = expression_matrix.toarray().tolist()[0] #Converts the matrix extracted above to a list

        correlation, pvalue = scipy.stats.spearmanr(component_values, expression_list) #Run the spearman correlation
        gene = genes[counter] #Get the gene from the list using index
        counter += 1
        #if counter % 1500 == 0:
            #print(expression_list[0:5])
            #print("FOLLOW PROGRESS : ", counter)
        res_list.append([gene,correlation,pvalue]) #Store the results to a list


    df = pd.DataFrame(res_list, columns=["gene", 'correlation','pvalue'])#Store the results from the list to a dataframe
    df = df.set_index('gene') #Change the gene name as index
    df.index.names = [None] #Remove the index column name
    df = df.sort_values(by = ["correlation"], ascending=False) #Sort values by correlation
    print("Length : ", len(df))
    os.mkdir("component_" + str(i))
    df.to_csv("component_" + str(i) + "/component_" + str(i)+ "correlation.csv", sep="\t") #Write correlation results to csv
    correlation_dict.update( {str(i) : df} ) #Stores the terminal state and correlation table to a dictionary
    print(df.head())

#Input directory for downstream analysis
directory = "diffusion_comp_correlation/" + dataset

#Run enrichr for diffusion components
for key, value in correlation_dict.items():
    component_name = "component_" + key
    spearman_res = value

    print(spearman_res)
    #Remove NA values
    ranked_list = spearman_res.iloc[:,0].dropna()

    #Get significant genes
    spearman_res = spearman_res.loc[spearman_res['pvalue'] <= 0.05]
    #Be sure that we use only genes with positive correlation
    spearman_res = spearman_res.loc[spearman_res['correlation'] > 0]

    #Convert the genes to list, use only most correlating gene
    gene_list = spearman_res.index.values.tolist()[0:int(num_genes)]
    print(gene_list[0:10])
    #COnvert background genes to set
    background_set = adata.raw.var.index.to_list()


    #Specify the output dir
    out_dir=component_name + "/"+ component_name + "enrichr"

    #Run the enrichr
    enr = gsea.enrichr(gene_list=gene_list,
                    gene_sets=gene_set,
                    background = background_set,
                    outdir = out_dir)

#Run GSEA for diffusion components
for key, value in correlation_dict.items():
    component_name = "component_" + key
    spearman_res = value

    #Ranked list contains the genes and correlation as a ranking metric, highest correlation first
    #Remove NaN values
    ranked_list = spearman_res.iloc[:,0].dropna()
    print(ranked_list)
    print(len(ranked_list))

    #Specify the output dir
    out_dir=component_name + "/"+ component_name + "gsea"

    #Run gsea
    gsea_res = gsea.prerank(rnk=ranked_list,
                    gene_sets=gene_set,
                    permutation_num=1000, # reduce number to speed up testing
                    outdir = out_dir, format="png", seed=6,
                    graph_num=6, min_size=15)

    #Plot the results
    gsea_res.res2d.sort_index().head()
    terms = gsea_res.res2d.index
