#!/usr/bin/env nextflow

process run_louvain_supercell {
    cpus 8
    memory '50 GB'
    time '12 hours'
    module 'R/4.2.3'
    publishDir "${projectDir}/out_supercell", mode: 'copy'

    input:
    tuple val(row), val(gamma)

    output:
    path(out_dur)
    path(out_dat)
    
    script:
    raw_data = "${projectDir}/data/supercell_mat_gamma${gamma}.csv"
    out_dat = "cluster_louvain_supercell_gamma${gamma}_k${row.k}.csv"
    out_dur = "duration_louvain_supercell_gamma${gamma}_k${row.k}.txt"
    
    """

    run_louvain_supercell.R \
        $raw_data \
        ${row.k} \
        ${gamma}
    """
}

process run_louvain_singlecell {
    cpus 8
    memory '200 GB'
    time '48 hours'
    module 'R/4.2.3'
    publishDir "${projectDir}/out_singlecell", mode: 'copy'

    input:
    val(row)

    output:
    path(out_dur)
    path(out_dat)
    
    script:
    raw_data = "${projectDir}/data/cell_dat_asinh.csv"
    out_dur = "duration_louvain_singlecell_k${row.k}.txt"
    out_dat = "cluster_louvain_singlecell_k${row.k}.csv"
    
    """

    run_louvain_singlecell.R \
        $raw_data \
        ${row.k}
    """
}

workflow {
    params_ch = Channel
        .fromPath("${projectDir}/data/params.csv")
        .splitCsv(header:true)
        
    // Max walltime set to 48 hours, the maximum allowed in Milton.
    run_louvain_singlecell(params_ch)

    gamma_ch = Channel.from(["10", "20", "30", "40"])
    // gamma_ch = Channel.from(["10"])
    combined_ch = params_ch.combine(gamma_ch)
    
    run_louvain_supercell(combined_ch)

}