#!/usr/bin/env nextflow

process run_cytofruv_supercell {
    cpus 12
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
    panel_data = "${projectDir}/data/Panel.csv"
    meta_data = "${projectDir}/data/Metadata.csv"
    out_dur = "cytofruv_supercell_gam${gamma}_setting${row.setting}.txt"
    out_dat = "dat_cytofruv_supercell_gam${gamma}_setting${row.setting}.csv"
    
    """

    run_cytofruv_supercell.R \
        $raw_data \
        $panel_data \
        ${row.meta} \
        ${row.k} \
        ${row.setting} \
        $meta_data \
        ${gamma}
    """
}

process run_cytofruv_singlecell {
    cpus 12
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
    panel_data = "${projectDir}/data/Panel.csv"
    out_dur = "cytofruv_singlecell_setting${row.setting}.txt"
    out_dat = "dat_cytofruv_singlecell_setting${row.setting}.csv"
    
    """

    run_cytofruv_singlecell.R \
        $raw_data \
        $panel_data\
        ${row.meta} \
        ${row.k} \
        ${row.setting} 
    """
}

workflow {
    params_ch = Channel
        .fromPath("${projectDir}/data/params.csv")
        .splitCsv(header:true)
        
    // Max walltime set to 48 hours, the maximum allowed in Milton.
    run_cytofruv_singlecell(params_ch)

    gamma_ch = Channel.from(["10", "20", "30", "40"])
    // gamma_ch = Channel.from(["10"])
    combined_ch = params_ch.combine(gamma_ch)
    
    run_cytofruv_supercell(combined_ch)

}