# Pipeline for identifying the underlying biological processes of the inferred trajectory
This is a trajectory inference pipeline, which first constructs a trajectory with either Palantir or Slingshot, then identifies the genes associated to each inferred trajectory lineage, and finally runs GSEA and Enrichr to the lineage associated genes. 

## Description
The pipeline is a nextflow (21.10.6) pipeline that combines scripts written in R (4.1.3) and Python (3.8.13). The trajectory inference methods of pipeline are choosed based on a comparison study of different trajectory methods. Slingshot (R package, 2.2.0), Palantir (Python package, 1.1), Paga from Scanpy package (Python package, 1.9.1)  and Monocle3 (R package, 1.0.0) trajectory methods was compared to see, which of these works best with our data. Out of the methods Slingshot, Palantir and Paga was included to this pipeline, but the downstream analysis in the pipeline is limited only for the trajectories inferred from Slingshot and Palantir methods. 

For the inferred Slingshot and Palantir trajectories the pipeline uses two different approaches to identify the genes associated to particular trajectory lineage. For the inferred Slingshot trajectory associationTest function from tradeSeq (R package, 1.8.0) is used to identify the genes associated to each trajectory lineage. For the inferred Palantir trajectory spearmanr function from  scipy (Python package, 1.8.1) is used to run Spearman correlation to identify the genes that correlates with specific lineage. To identify the biological processes underlying the trajectory lineage the pipeline contains scripts for GSEA and Enrichr from gseapy (Python package, 0.10.8). 

## Usage of the pipeline

## Results from the comparison study
