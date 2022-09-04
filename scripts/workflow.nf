//This script works between main.nf and processes.nf scripts

nextflow.enable.dsl=2
script_folder = "$baseDir/scripts"
include {run_slingshot} from "$script_folder/processes.nf"
include {run_tradeSeq} from "$script_folder/processes.nf"
include {run_gsea} from "$script_folder/processes.nf"
include {run_enrichr} from "$script_folder/processes.nf"
include {run_palantir} from "$script_folder/processes.nf"
include {run_paga} from "$script_folder/processes.nf"
include {run_diffusion_comp_correlation} from "$script_folder/processes.nf"


workflow wf_slingshot{
    take:
        input_data_path
        dataset
        start_cluster
        extend_value
        slingshot_clusters
    main:
        run_slingshot(input_data_path, dataset, start_cluster, extend_value, slingshot_clusters)
    emit:
        slingshot_results = run_slingshot.out.slingshot_results
        var_genes = run_slingshot.out.var_genes
        h5seurat_file = run_slingshot.out.h5seurat_file
}


workflow wf_palantir{
    take:
        input_data_path
        dataset
        palantir_start_cell
        palantir_clusters
    main:
        run_palantir(input_data_path, dataset, palantir_start_cell, palantir_clusters)
    emit:
        diffusion_components = run_palantir.out.palantir_diff_comp
        trajectory_correlation = run_palantir.out.palantir_tra_correlation
        ranked_list = run_palantir.out.palantir_ranked_list
        gene_list = run_palantir.out.palantir_gene_list
}

workflow wf_paga{
    take:
        input_data_path
        dataset
    main:
        run_paga(input_data_path, dataset) 
}


workflow wf_tradeseq{
    take:
        h5seurat_file
	trajectory_results
        dataset
        var_genes       
    main:
        run_tradeSeq(h5seurat_file,trajectory_results, dataset, var_genes)
    emit:
        tradSeq_results = run_tradeSeq.out.associationTest_results
        ranked_list = run_tradeSeq.out.tradeSeq_ranked_list
        gene_list = run_tradeSeq.out.tradeSeq_gene_list     
}

workflow wf_diffusion_comp_correlation{
    take:
        input_data
        dataset
        gene_set
        num_genes
        diffusion_components
    main:
        run_diffusion_comp_correlation(input_data, dataset, gene_set, num_genes, diffusion_components)
}

workflow wf_gsea{
    take:
        ranked_list
        dataset
        gene_set
        ti_method
    main:
        run_gsea(ranked_list, dataset, gene_set, ti_method)

}

workflow wf_enrichr{
    take:
        gene_list
        input_data
        dataset
        gene_set
        ti_method
        num_genes
    main:
        run_enrichr(gene_list, input_data, dataset, gene_set, ti_method, num_genes)
}


