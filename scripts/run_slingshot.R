"Title: run_slingshot.R
Date: 2021-08-06
Author: Erno Hänninen
Description:
    This script runs slingshot for the input dataset and plots the trajectories and pseudotime.
    In addition the script extracts the variable genes from sce_object (same as adata.var.index)
List of functions:
    -    
List of non standard modules:
    -    stringr, zellkonverter, SingleCellExperiment,
         Seurat, SeuratDisk, RColorBrewer,
         traviz, slingshot, devtools
Error handling:
    -
Procedure:
    - This script is called from processes.nf"

library(devtools)
devtools::install_github("mojaveazure/seurat-disk")
devtools::install_github("theislab/zellkonverter")
library(SeuratDisk)
library(zellkonverter)
library(stringr)
library(SingleCellExperiment)
library(Seurat)
library(RColorBrewer)
library(grDevices) #Base package
library(traviz)
library(slingshot)

#Input arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[[1]]
dataset <- args[[2]]
start_cluster <- args[[3]]
extend_value <- args[[4]]
clusters <- args[[5]]
#Extract the input path from file_path
file_name <- str_split(basename(input_file), ".h5ad")[[1]][1] #Get file name without name extension


#Read the h5ad
#The sce_object contains only the variable features
sce_object <- readH5AD(input_file)

#If h5seurat file doesn't exist convert h5ad file to h5seurat
if(!file.exists(paste0(paste0(dirname(input_file),"/"), paste0(file_name,".h5seurat")))){
  Convert(input_file, dest = "h5seurat",  overwrite=TRUE) #Convert h5ad to seurat
  #Read the h5seurat file
  seurat_object <- LoadH5Seurat(paste0(dirname(input_file), "/", file_name, ".h5seurat"), assays="RNA")
}else{ #Else the h5seurat file already exists, no conversion needed
  seurat_object <- LoadH5Seurat(paste0(dirname(input_file), "/", file_name, ".h5seurat"), assays="RNA") #read the file
}

#Saving the converted file to the work directory, so that we can use it in other scripts
SaveH5Seurat(seurat_object)
#This data contains all the observations
raw_data <- as.SingleCellExperiment(seurat_object)


#RUNNING SLINGSHOT
#Output is a SlingshotDataSet object 
results <- slingshot(sce_object, reducedDim = 'X_umap', clusterLabels=colData(sce_object)[[clusters]], start.clus=start_cluster, extend=extend_value)
#Saving the slignshot result to file
saveRDS(results, file = paste0(dataset, "_slingshot_results.rds"))


###################### PLOTTING ##########################
#Plot the trajecotries, color with user specified clusters
plotGeneCount(SlingshotDataSet(results), clusters=colData(sce_object)[[clusters]])

#Visualizing the inferred lineage for the single-trajectory data with points colored by the pseudotime
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(results$slingPseudotime_1, breaks=100)]
plot(reducedDims(results)$X_umap, col = plotcol, pch=16)
lines(SlingshotDataSet(results), lwd=2, col='black')

#Extract the variable gene names from the object
genes <- rownames(sce_object)
print(length(genes))
print(typeof(genes))
#Write the gene names to table
#These genes are used for tradeSeq fitGAM function
writeLines(genes, "genes.txt")
