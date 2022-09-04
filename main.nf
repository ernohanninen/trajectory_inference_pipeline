//This script is connected to pipeline.config file and controls which processes the program executes
//This scripts calls the workflow.nf file, which controls the workflow by calling the processes to run from processes.nf file. The processes.nf file calls the actual scripts.

nextflow.enable.dsl=2

script_folder = "$baseDir/scripts"

include {wf_slingshot} from "$script_folder/workflow.nf"
include {wf_palantir} from "$script_folder/workflow.nf"
include {wf_paga} from "$script_folder/workflow.nf"
include {wf_tradeseq} from "$script_folder/workflow.nf"
include {wf_diffusion_comp_correlation} from "$script_folder/workflow.nf"
include {wf_gsea} from "$script_folder/workflow.nf"
include {wf_enrichr} from "$script_folder/workflow.nf"


workflow {
    
    if(!params.skip_slingshot){
        wf_slingshot(params.input_data, params.dataset, params.slingshot_start_cluster, params.slingshot_extend_value, params.slingshot_clusters)       
    }    
    if(!params.skip_slingshot && !params.skip_tradeSeq){
        wf_tradeseq(wf_slingshot.out.h5seurat_file, wf_slingshot.out.slingshot_results, params.dataset, wf_slingshot.out.var_genes)
    }
    if(!params.skip_slingshot && !params.skip_tradeSeq && !params.skip_gsea){
        wf_gsea(wf_tradeseq.out.ranked_list, params.dataset, params.gene_set,"slingshot")
    }
    if(!params.skip_slingshot && !params.skip_tradeSeq && !params.skip_enrichr){
        wf_enrichr(wf_tradeseq.out.gene_list, params.input_data, params.dataset, params.gene_set, "slingshot", params.enrichr_num_genes)
    }
    if(!params.skip_palantir){
        wf_palantir(params.input_data, params.dataset, params.palantir_start_cell, params.palantir_clusters)
    }
    if(!params.skip_palantir && !params.skip_diffusion_comp_correlation){
        wf_diffusion_comp_correlation(params.input_data, params.dataset, params.diff_comp_gene_set, params.enrichr_diff_comp_num_genes,  wf_palantir.out.diffusion_components)
    }
    if(!params.skip_palantir && !params.skip_gsea){
        wf_gsea(wf_palantir.out.ranked_list, params.dataset, params.gene_set, "palantir")
    }
    if(!params.skip_palantir && !params.skip_enrichr){
        wf_enrichr(wf_palantir.out.gene_list, params.input_data, params.dataset, params.gene_set, "palantir", params.enrichr_num_genes)
    }
    if(!params.skip_paga){
        wf_paga(params.input_data, params.dataset)
    }
}
