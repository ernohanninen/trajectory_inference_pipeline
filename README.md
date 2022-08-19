# Pipeline for identifying the underlying biological processes of the inferred trajectory
This is a trajectory inference pipeline, which first constructs a trajectory with either Palantir or Slingshot, then identifies the genes associated to each inferred trajectory lineage, and finally runs GSEA and Enrichr to the lineage associated genes. 

## Description of the pipeline
The pipeline is a nextflow (21.10.6) pipeline that combines scripts written in R (4.1.3) and Python (3.8.13). The trajectory inference methods of pipeline are choosed based on a comparison study of different trajectory methods. Slingshot (R package, 2.2.0), Palantir (Python package, 1.1), Paga from Scanpy package (Python package, 1.9.1)  and Monocle3 (R package, 1.0.0) trajectory methods was compared to see, which of these works best with our data. Out of the methods Slingshot, Palantir and Paga was included to this pipeline, but the downstream analysis in the pipeline is limited only for the trajectories inferred from Slingshot and Palantir methods. 

For the inferred Slingshot and Palantir trajectories the pipeline uses two different approaches to identify the genes associated to particular trajectory lineage. For the inferred Slingshot trajectory associationTest function from tradeSeq (R package, 1.8.0) is used to identify the genes associated to each trajectory lineage. For the inferred Palantir trajectory spearmanr function from  scipy (Python package, 1.8.1) is used to run Spearman correlation to identify the genes that correlates with specific lineage. To identify the biological processes underlying the trajectory lineage the pipeline contains scripts for GSEA and Enrichr from gseapy (Python package, 0.10.8). 

## Usage of the pipeline
Setup the pipeline:
1. Open terminal
2. Clone the repository to your computer
```
git clone https://github.com/ernohanninen/trajectory_inference_pipeline.git
```
3. Navigate to the pipeline environments directory
```
cd trajectory_inference_pipeline/environments/
```
4. Setup the conda environment for R scripts
```
conda env create -f Renv.yml
```
5. Setup the conda environment for Python scripts
```
conda env create -f PYenv.yml
```
6. Navigate back to the trajectory_inference_pipeline directory
```
cd ..
```
7. Open the pipeline.config file in a text editor. See details from the "Update configuration file" -section, how to fill in the pipeline.config file.
```
nano pipeline.config
```
8. Run the pipeline
```
nextflow main.nf -c pipeline.config
```

The results of the analysis are stored in a user specified output file. 

## Update configuration file
Before every analysis the configuration file (pipeline.config) needs to be updated. The pipeline.config file controls which scripts are executed during the analysis. 
 - params.input_data: Submit a h5ad file. See below the requirements for the input h5ad file.
 - params.dataset: dataset name, this information is used to neame the output directory.
 - params.output_folder: folder where the output is stored.
 - params.skip_slingshot: if false the Slingshot trajectory inference method is executed.
 - params.skip_palantir: if false the Palantir trajectory inference method is executed.
 - params.skip_paga: if false PAGA the trajectory inference method is executed.
 - params.slingshot_start_cluster: the cluster where the slingshot trajectory origins.
 - params.slingshot_extend_value: specifies the method how to handle root and leaf clusters when constructing the slingshot trajectory curves. Possible values: "y", "n" and "pc1". More information see extend from: https://rdrr.io/github/kstreet13/slingshot/man/slingParams.html
 - params.sligshot_clusters: clusters to use to construct and visualize the slingshot trajectory. These should be clusters from the h5ad data.
 - params.palantir_start_cell: Start cell for the Palantir trajectories.
 - params.palantir_clusters: clusters to be used to visualize the Palantir results. These should be clusters from the h5ad data.
 - params.skip_tradeSeq: if skip_slingshot is false and skip_tradeSeq is false the script which runs fitGAM and associationTest functions from tradeSeq package is executed.
- params.skip_gsea = if skip_palantir and skip_gsea is false, GSEA analysis is performad for Palantir lineages. If skip_slingshot, skip_tradeSeq and skip_gsea is false GSEA analysis is performed for Slingshot trajectories.
- params.skip_enrichr = if skip_palantir and skip_enrichr is false, Enrichr analysis is performad for Palantir lineages. If skip_slingshot, skip_tradeSeq and skip_enrichr is false Enrichr analysis is performed for Slingshot trajectories.
- params.gene_set: specifies the gene set to be used in the analysis. The gene sets, should be loded in gene_sets folder.
- params.ernichr_num_genes: specifies the number of most associated genes to run the Enrichr with.

## Requirements for h5ad input file
- Anndata (adata) object with variable features in adata.var.
- adata.raw contains the raw counts.
- adata.var.index should return a list of variable genes.
- The input data needs to be dimensionality reduced (PCA and UMAP) so that in adata.obsm there is column named X_pca and X_umap.
- The input data needs to be clustered so that in adata.obs there is a column for cluster(s). Use on of these clusters in pipeline.config file as an argument for params.sligshot_clusters or params.slingshot_clusters.
- The adata should have neigbourhood graph computed.

## Description of the pipeline output files
#### Slingshot and tradeSeq
run_slingshot.R
- AnnData.h5Seurat contains file with seurat object, which is converted from h5ad file. 
- genes.txt file contains the variable genes. (adata.var.index)
- Rplots.pdf contains the trajectories plotted in a UMAP. One plot colored by user specified clusters and the other plot colored by pseudotime.
-  {dataset}_slingshot_results.rds contains the result from slingshot function.

run_tradeSeq.R script returns
- {dataset}_fitGam.rds contains the output from the fitGAM function. 
- {dataset}_associationTest_lineage_{num}.csv: contains the results from associationTest.
- {dataset}_geneList_lineage_{num}.txt: contains the top lineage associated and significant genes for Enrichr.
- {dataset}_ranked_lineage_{num}.rnk: contains all the genes from assiciationTest (genes with NA values have been filtered), even the non-significant, gene with highest Wald stat value is an top of the list. This is a ranked list for GSEA.

#### Palantir and Spearman correlation
run_palantir.py
- Figures of diffusion components, trajectories, start cell and terminal state
- correlation_{dataset}_lineage_{num}.csv: correlation table from spearman correlation. Gene with highest correlation is first.
- {dataset}_geneList_lineage_{num}.txt: Contains the top lineage correlating and significant genes for Enrichr.
- {dataset}_ranked_lineage_{num}.rnk: contains all the genes from Spearman correlation, even the non-significant, gene with highest correlation value is an top of the list. This is a ranked list for GSEA.

#### GSEA and Enrichr
- run_gsea.py outputs GSEA plots
- run_enrichr.py outputs Enrichment plots
