# Pipeline for identifying the underlying biological processes of the inferred trajectory
This is a trajectory inference pipeline, which first constructs a trajectory with either Palantir or Slingshot, then identifies the genes associated to each inferred trajectory lineage, and finally runs GSEA and Enrichr to the lineage associated genes. 

## Description of the pipeline
The pipeline is a nextflow (21.10.6) pipeline that combines scripts written in R (4.1.3) and Python (3.8.13). The trajectory inference methods of pipeline are choosed based on a comparison study of different trajectory methods. Slingshot (R package, 2.2.0), Palantir (Python package, 1.1), Paga from Scanpy package (Python package, 1.9.1)  and Monocle3 (R package, 1.0.0) trajectory methods was compared to see, which of these works best with our data. Out of the methods Slingshot, Palantir and Paga was included to this pipeline, but the downstream analysis in the pipeline is limited only for the trajectories inferred from Slingshot and Palantir methods. 

For the inferred Slingshot and Palantir trajectories the pipeline uses two different approaches to identify the genes associated to particular trajectory lineage. For the inferred Slingshot trajectory associationTest function from tradeSeq (R package, 1.8.0) is used to identify the genes associated to each trajectory lineage. For the inferred Palantir trajectory spearmanr function from  scipy (Python package, 1.8.1) is used to run Spearman correlation to identify the genes that correlates with specific lineage. To identify the biological processes underlying the trajectory lineage the pipeline contains scripts for GSEA and Enrichr from gseapy (Python package, 0.10.8). 

## Usage of the pipeline
Setup the pipeline:
1. Open terminal
2. Clone the repository to your computer
3. Navigate to the pipeline environments directory
4. Setup the conda environments
5. Navigate back to the pipeline directory
6. Open the pipeline.config file in a text editor. See details from pipeline.config file section, how to fill in the pipeline.config file.
```
nano pipeline.config
```
7. Run the pipeline
```
nextflow main.nf -c pipeline.config
```

The results of the analysis are stored in a user specified output file. 

## Description of the pipeline output
#### Slingshot and tradeSeq
run_slingshot.R script returns AnnData.h5Seurat, genes.txt, Rplots.pdf and dataset_slingshot_results.rds. AnnData.h5Seurat contains file with seurat object, which is converted from h5ad file. genes.txt file contains the variable genes. Rplots.pdf contains the trajectories plotted in a UMAP. One plot colored by user specified clusters and the other plot colored by pseudotime. {dataset}_slingshot_results.rds contains the result from slingshot function.

run_tradeSeq.r script returns
{dataset}_fitGam.rds contains the output from the fitGAM function. 
{dataset}_associationTest_lineage_{num}.csv contains the results from associationTest
