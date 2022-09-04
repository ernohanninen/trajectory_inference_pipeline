"Title: run_tradeSeq.R
Date: 2021-08-10
Author: Erno Hänninen
Description:
    - Runs tradeSeq association test function for slingshot result
    - For each lineage writes ranked list and gene list to file
    
List of functions:
    -    
List of non standard modules:
    -     tidyr, SingleCellExperiment,
         Seurat, SeuratDisk,
         slingshot, dplyr, tradeSeq
Error handling:
    -
Procedure:
    - This script is called from processes.nf"


library(tidyr)
library(SingleCellExperiment)
library(tradeSeq)
library(Seurat)
library(SeuratDisk)
library(grDevices) #Base package
library(slingshot)
library(dplyr)

#Get the input arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[[1]]
trajectory_results <- args[[2]]
dataset <- args[[3]]
var_genes_path <- args[[4]]

#Read the original dataset and convert to SCE dataset
data <- LoadH5Seurat(input_file, assays="RNA")
raw_data <- as.SingleCellExperiment(data)
print(raw_data)

#Read the slingshot results stored in RDS object
results <- readRDS(trajectory_results)

#Extract input for fitGam function
counts <- assays(raw_data)$counts %>% as.matrix() # extract count matrix
sling_pseudotime <- slingPseudotime(results, na = FALSE) #Extract the Slingshot pseudotime values
sling_weights = slingCurveWeights(results) #Extract Slingshot weights

#Get the genes from file to use as an input for fitGAM function
hvg_genes <- readLines(var_genes_path)
#hvg_genes <- hvg_genes[1:10]


#fitGAM function fits the NB-GAM model for each gene (A negative binomial generalized additive model)
#Returns singleCellExperiment object where tradeSeq results are stored
fitGam_res <- fitGAM(counts=counts,pseudotime=sling_pseudotime, cellWeights=sling_weights, genes = hvg_genes,sce=TRUE)
saveRDS(fitGam_res, file=paste0(dataset, "_fitGam.rds")) #Save the fitgam results to file

#associationTest checks wheter gene expression is associated with a particular lineage
#argument lineages=TRUE, test for all lineages independently
#Returns a matrix with wald statistic, the degrees of freedom and the (unadjusted) p-value for each gene
assoRes <- associationTest(fitGam_res, lineages=TRUE)


#slingLineages returns a list of lineages, represented by ordered sets of clusters.
lin <- slingLineages(results)

#Get the genes associated to each lineage
#Filter the dataframes and make a new df for each lineage
for(i in 1:length(lin)){
    print(paste0("LINEAGE: ", lin))
    #Remove na values of pvalue column  
    assign("temp", assoRes %>% drop_na(paste0("pvalue_", i)))
    #adjust p-values and create a new vector containing the adjusted p-values
    assign("p_val_adj", p.adjust(temp[[paste0("pvalue_", i)]], "fdr"))
    #Add the new vector to the temp df
    temp[paste0("p_val_adj_", i)]<- p_val_adj
    #Extract the lineage specific values to a new df
    temp <- temp[, c(paste0("waldStat_", i), paste0("df_",i), paste0("pvalue_",i), paste0("p_val_adj_",i))]
    #Order the df based on waldStat column
    temp <- temp[order(temp[paste0("waldStat_", i)], decreasing=TRUE),] 
  
   
    if(!is.null(temp)){
        print("Output from association test")
        print(temp)
        #Write lineage specific associationTest results to file 
        write.table(temp,paste0(dataset,"_associationTest_", paste0("lineage_",i), ".csv"),sep="\t", col.names = NA, row.names = TRUE)
        
        #Prepare the ranked list for GSEA
        #Gene name as row index, waldStat as a value
        ranked_list <- temp[,paste0("waldStat_", i), drop=FALSE]
        names(ranked_list)<-NULL
        print("RANKED LIST")
        print(ranked_list)
        #Write the ranked list to file
        write.table(ranked_list,paste0(dataset,"_ranked_", paste0("lineage",i), ".rnk"))
  
        #Prepare the gene list for enrichr
        #Use only significant genes
        filtered_temp <- temp[which(temp[,paste0("p_val_adj_", i)]<=0.05),]
        #Get the row names (gene names) from the filtered table
        gene_list <- rownames(filtered_temp)
        print("GENE LIST")
        print(gene_list)
        #Write genes line by line to file
        writeLines(gene_list, paste0(dataset,"_geneList_", paste0("lineage",i), ".txt"))

        #Rename the the dataframe
        assign(paste0("lineage", i), temp)
    }   
}

