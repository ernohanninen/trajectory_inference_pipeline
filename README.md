# Pipeline for identifying the underlying biological processes of the inferred trajectory
This is a trajectory inference pipeline, which first constructs a trajectory with either Palantir or Slingshot, then identifies the genes associated with each inferred trajectory lineage, and finally runs GSEA and Enrichr to the lineage-associated genes. 

## Description of the pipeline
The pipeline is a Nextflow (21.10.6) based pipeline that combines scripts written in R (4.1.3) and Python (3.8.13). The trajectory inference methods of the pipeline were chosen based on a comparison study of different trajectory methods. Slingshot (R package, 2.2.0), Palantir (Python package, 1.1.0), Paga from Scanpy (Python package, 1.9.1), and Monocle3 (R package, 1.0.0) trajectory methods was compared to see, which of these works best with the data used. Out of the methods Slingshot, Palantir and Paga were included in this pipeline, but the downstream analysis in the pipeline is limited only to the trajectories inferred from Slingshot and Palantir methods. The data I used in this study didn’t have developmental trajectories, instead with trajectory inference methods I was trying to identify continuous cell states.

For the inferred Slingshot and Palantir trajectories, the pipeline uses two different approaches to identify the genes associated with a particular trajectory lineage. For the inferred Slingshot trajectory, associationTest function from tradeSeq (R package, 1.8.0) is used to identify the genes associated with a trajectory lineage. For the inferred Palantir trajectory, spearmanr function from  scipy (Python package, 1.8.1) is used to run Spearman correlation to identify the genes that correlate with a trajectory lineage. To identify the biological processes underlying the lineage the pipeline contains scripts for GSEA and Enrichr analysis from gseapy (Python package, 0.10.8). 

The pipeline is developed on Rocky Linux 8.6 and Conda 4.13.0 was used to manage the packages.

## Usage of the pipeline
To run the pipeline Conda and Nextflow are required. 

Setup the pipeline:
1. Open terminal:
2. Clone the repository to your computer:
```
git clone https://github.com/ernohanninen/trajectory_inference_pipeline.git
```
3. Navigate to the pipeline environments directory:
```
cd trajectory_inference_pipeline/environments/
```
4. Setup the conda environment for R scripts:
```
conda env create -p ./Renv -f Renv.yml
```
5. Setup the conda environment for Python scripts:
```
conda env create -p ./PYenv -f PYenv.yml
```
6. Navigate back to the trajectory_inference_pipeline directory:
```
cd ..
```
7. Create a directory where to store the input data
```
mkdir data
```
8. Open the pipeline.config file in a text editor. See details from the "Update configuration file" -section, how to update the pipeline.config file:
```
nano pipeline.config
```
9. Run the pipeline:
```
nextflow main.nf -c pipeline.config
```

The results of the analysis are stored in a user specified output file. 

## Update configuration file
Before every analysis the configuration file (pipeline.config) needs to be updated. The configuration file controls which scripts are executed during the analysis. 
 - params.input_data: Submit a h5ad file. See below the requirements for the input h5ad file.
 - params.dataset: dataset name, this information is used to name the output directory.
 - params.output_folder: folder where the output of the pipeline is stored.
 - params.skip_slingshot: if false Slingshot trajectory inference method is executed.
 - params.skip_palantir: if false Palantir trajectory inference method is executed.
 - params.skip_paga: if false PAGA trajectory inference method is executed.
 - params.slingshot_start_cluster: the cluster where the slingshot trajectory origins.
 - params.slingshot_extend_value: specifies the method how to handle the root and leaf clusters when constructing the Slingshot trajectory curves. Possible values: "y", "n" and "pc1". More information see "extend" from: https://rdrr.io/github/kstreet13/slingshot/man/slingParams.html
 - params.sligshot_clusters: clusters to be used to construct and visualize the Slingshot trajectory. These needs to be clusters from the h5ad data.
 - params.palantir_start_cell: Start cell for the Palantir trajectories.
 - params.palantir_clusters: clusters to be used to visualize the Palantir results. These needs to be clusters from the h5ad data.
 - params.skip_tradeSeq: if skip_slingshot is false and skip_tradeSeq is false, the script which runs fitGAM and associationTest functions from tradeSeq package is executed.
- params.skip_gsea = if skip_palantir and skip_gsea are false, GSEA analysis is performed for the Palantir trajectory. If skip_slingshot, skip_tradeSeq and skip_gsea are false, GSEA analysis is performed for the Slingshot trajectory.
- params.skip_enrichr = if skip_palantir and skip_enrichr are false, Enrichr analysis is performed for Palantir lineages. If skip_slingshot, skip_tradeSeq and skip_enrichr are false, Enrichr analysis is performed for Slingshot trajectories.
- params.gene_set: specifies the path to the gene set used in the Enrichr. The gene sets, should be loaded in gene_sets folder. See the available gene sets: https://maayanlab.cloud/Enrichr/#libraries
- params.enrichr_num_genes: specifies the number of most associated significant genes to run the Enrichr with.
- params.skip_diffusion_comp_correlation: if skip_palantir and skip_diffusion_comp_correlation are false spearman correlation is exectued to the diffusion components. In addition the scripts runs GSEA and Enrichr to the result from the spearman correlation.
params.diff_comp_gene_set: specifies the path to the gene set used in the Enrichr and GSEA for the genes correlating with diffusion component. The gene sets, should be loaded in gene_sets folder. See the available gene sets: https://maayanlab.cloud/Enrichr/#libraries
params.enrichr_diff_comp_num_genes: specifies the number of most associated significant genes to run the Enrichr with for the genes correlating with diffusion component.

## Notes of the pipeline
- Palantir diffusion components are computed with all observations in the data, making the GSEA function for diffusion components computationally demanding. If the dataset is big the GSEA function requires a lot of computational resources. 
- Even though the tradeSeq fitGAM function is executed only with the variable genes it also requires a lot of computational resources.

## Requirements for h5ad input file
- Anndata (adata) object with variable features in adata.var.
- adata.raw contains the raw counts.
- adata.var.index contains a list of variable genes.
- The input data needs to be dimensionality reduced (PCA and UMAP), so that in adata.obsm there are columns named X_pca and X_umap.
- The input data needs to be clustered, so that in adata.obs there is a column for cluster(s). Use one of these clusters in pipeline.config file as an argument for params.sligshot_clusters or params.slingshot_clusters.
- The adata should have neighborhood graph computed.

## Description of the pipeline output files
#### Slingshot and tradeSeq
run_slingshot.R script returns
- AnnData.h5Seurat contains file with Seurat object, which is converted from h5ad file. 
- genes.txt file contains the variable genes. (adata.var.index)
- Rplots.pdf contains the trajectories plotted in a UMAP. One plot colored by user specified clusters and the other plot colored by pseudotime.
- {dataset}_slingshot_results.rds contains the result from slingshot function.

run_tradeSeq.R script returns
- {dataset}_fitGam.rds contains the output from the fitGAM function. 
- {dataset}_associationTest_lineage_{num}.csv: contains the results from associationTest.
- {dataset}_geneList_lineage_{num}.txt: contains the top lineage associated significant genes for Enrichr.
- {dataset}_ranked_lineage_{num}.rnk: contains all the genes from assiciationTest (genes with NA values have been filtered away), even the non-significant. Gene with highest Wald stat value is an top of the list. This is the ranked list for GSEA.

#### Palantir and Spearman correlation
run_palantir.py script returns
- Figures of diffusion components, trajectories, start cell and terminal state
- File for the diffusion components
- correlation_{dataset}_lineage_{num}.csv: correlation table from spearman correlation. Gene with highest correlation is first.
- {dataset}_geneList_lineage_{num}.txt: Contains the top lineage correlating significant genes for Enrichr.
- {dataset}_ranked_lineage_{num}.rnk: contains all the genes from Spearman correlation, even the non-significant. Gene with highest correlation value is an top of the list. This is the ranked list for GSEA.

#### GSEA and Enrichr
- run_gsea.py outputs GSEA plots
- run_enrichr.py outputs Enrichment plots

#### Diffusion component
run_diffusion_comp_correlation.py script returns
- spearman correlation results for each diffusion component
- GSEA and ENRICHR plots
