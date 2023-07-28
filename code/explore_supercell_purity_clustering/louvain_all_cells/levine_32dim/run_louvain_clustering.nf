#!/usr/bin/env nextflow

process cluster_all_louvain {
    cpus "$params.max_cores"
    memory '120 GB'
    time '30 hours'
    module 'R/4.2.3'
    publishDir "$params.publish_dir/all_louvain", mode: 'copy'

    input:
    val(k_val)
    
    output:
    path(outdir)
    
    script:
    outdir = "k" + k_val
    
    """
    mkdir $outdir

    cluster_all_louvain.R \
        $params.raw_data \
        $params.marker_info_file \
        $k_val \
        $outdir \
        $params.max_cores \
        $params.benchmark_ntimes 
    
    """
}

workflow {
    louvain_params = Channel
        .fromPath(params.louvain_config)
        .splitCsv(header:true)
        .map {louvain_config -> louvain_config.get("k")}
        
    // Max walltime set to 48 hours, the maximum allowed in Milton.
    louvain_clusters = cluster_all_louvain(louvain_params)  
}