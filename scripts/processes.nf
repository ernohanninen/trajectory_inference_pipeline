#!/usr/bin/env nextflow
 
//This script calls the scripts

process run_slingshot {
    publishDir "$params.output_folder/slingshot_${dataset}/slingshot_trajectories", mode:'copy', overwrite: true
    conda "$projectDir/environments/Renv"
    input:
        path(input_data)
        val(dataset)
        val(start_cluster)
        val(extend_value)
        val(slingshot_clusters)
    output:
        path("Rplots.pdf")
        path("${dataset}_slingshot_results.rds", emit: slingshot_results)
        path("genes.txt", emit: var_genes)
        path("*.h5Seurat", emit: h5seurat_file)


    script:
    """
     Rscript $projectDir/scripts/run_slingshot.R \\
        $params.input_data \\
        $params.dataset \\
        $params.slingshot_start_cluster \\
        $params.slingshot_extend_value \\
        $params.slingshot_clusters \\
        . > run_slingshot.txt 2>&1
    """
}


process run_palantir{
    publishDir "$params.output_folder/palantir_${dataset}/palantir_trajectories", mode:'copy', overwrite: true

    conda "$projectDir/environments/PYenv"

    input:
        path(input_data)
        val(dataset)
        val(palantir_start_cell)
        val(palantir_clusters)

    output:
        path("figures/*")
        path("*_diff_components.csv", emit:palantir_diff_comp)
        path("correlation_*", emit: palantir_tra_correlation)
        path("*.rnk", emit:palantir_ranked_list)
        path("*lineage*.txt", emit:palantir_gene_list)

    script:
    """
    python3 $projectDir/scripts/run_palantir.py \\
        $params.input_data \\
        $params.dataset \\
        $params.palantir_start_cell \\
        $params.palantir_clusters \\
         > run_palantir.txt 2>&1
    """
}

process run_paga{
    publishDir "$params.output_folder/paga_${dataset}/paga_trajectories", mode:'copy', overwrite: true
    
    conda "$projectDir/environments/PYenv"
    input:
        path(input_data)
        val(dataset)

    output:
        path("figures/*")
        
    script:
    """
    python3 $projectDir/scripts/run_paga.py \\
        $params.input_data \\
        $params.dataset \\
         > run_paga.txt 2>&1
    """
}

process run_tradeSeq{
    publishDir "$params.output_folder/slingshot_${dataset}/tradeSeq", mode:'copy', overwrite: true
    conda "$projectDir/environments/Renv"
     
    input:
        path(h5seurat_file)
        path(trajectory_results)
        val(dataset)
        path(var_genes)

    output:
        path("${dataset}_fitGam.rds")
        path("${dataset}_associationTest_lineage_*", emit: associationTest_results)
        path("*.rnk", emit:tradeSeq_ranked_list)
        path("*lineage*.txt", emit:tradeSeq_gene_list)

 
    script:
    """
     Rscript $projectDir/scripts/run_tradeSeq.R \\
        $h5seurat_file \\
        $trajectory_results \\
        $params.dataset \\
        $var_genes \\
        . > run_tradeSeq.txt 2>&1
    """
}

process run_diffusion_comp_correlation{
    publishDir "$params.output_folder/palantir_${dataset}/diffusion_component_correlations", mode:'copy', overwrite: true
    conda "$projectDir/environments/PYenv"

    input:
        path(input_data)
        val(dataset)
        path(gene_set)
        val(num_genes)
        path(diffusion_components)

    output:
        path("component*/")

    script:
    """
     python3 $projectDir/scripts/run_diffusion_comp_correlation.py $params.input_data $params.dataset $params.diff_comp_gene_set $params.enrichr_diff_comp_num_genes $diffusion_components > run_diffusion_comp_correlation.txt 2>&1
    """
}

process run_gsea{
    publishDir "$params.output_folder/${ti_method}_${dataset}/gsea", mode:'copy', overwrite: true
    conda "$projectDir/environments/PYenv"

    input:
        each path(ranked_list)
        val(dataset)
        path(gene_set)
        val(ti_method)

    output:
        path("${dataset}*")

 
    script:
    """
     python3 $projectDir/scripts/run_gsea.py \\
        $ranked_list \\
        $params.dataset \\
        $params.gene_set \\
        $ti_method \\
         > run_gsea.txt 2>&1
    """
}



process run_enrichr{
    publishDir "$params.output_folder/${ti_method}_${dataset}/enrichr", mode:'copy', overwrite: true
    conda "$projectDir/environments/PYenv"
    input:
        each path(gene_list)
        path(input_data)
        val(dataset)
        path(gene_set)
        val(ti_method)
        val(num_genes)

    output:
        path("${dataset}*")

 
    script:
    """
     python3 $projectDir/scripts/run_enrichr.py \\
        $gene_list \\
        $params.input_data \\
        $params.dataset \\
        $params.gene_set \\
        $ti_method \\
        $params.enrichr_num_genes \\
         > run_enrichr.txt 2>&1
    """
}






