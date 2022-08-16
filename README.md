# trajectory_inference_pipeline: Pipeline for identifying the underlying biological processes of the inferred trajectory
This is a trajectory inference pipeline, which first constructs a trajectory with either Palantir or Slingshot, then finds the genes associated to each inferred trajectory lineage, and finally runs GSEA and Enrichr to the lineage associated genes. 
