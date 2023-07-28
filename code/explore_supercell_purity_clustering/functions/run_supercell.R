#' Run SuperCellCyto
#'
#' Run SuperCellCyto on a cytometry data using a set of gamma values,
#' then write out the supercell_expression_matrix and supercell_cell_map out
#' as csv files, the supercell_object as a RDS file, and the amount of time
#' taken to run supercell as a text file.
#'
#' @details
#' Few things that are specific for our benchmarking:
#'
#' 1. The markers that are fed into SuperCellCyto are all cell type markers that
#' have been undergone arc-sinh transformation with co-factor = 5.
#' 2. For each gamma value, SuperCellCyto is ran 2x, and the amount of time
#' taken to do each run of SuperCellCyto is recorded using the tictoc package.
#' 3. Number of workers for the Multicoreparam is set to the default value set
#' by BiocParallel. At the time of benchmarking, this was set to parallel::detectCores() - 2.
#'
#' For each iteration, BPPARAM containing the MulticoreParam object is deleted
#' and gc is called. Otherwise we will run into memory issues.
#'
#' @param cell_dat_file character.
#' The name of the csv file containing the cytometry data.
#' @param marker_info_file character. The name of the csv file containing the marker information for
#' the cytometry data.
#' This csv file must contain marker_name (name of the marker), and
#' marker_class (type of marker - is the marker used to characterise cell type (type) or
#' cell state (state)?).
#' @param gamma_values vector of numeric. The list gamma values to run SuperCellCyto with.
#' @param save_dir character. Directory where to store the output files.
#'
run_supercell <- function(cell_dat_file,
                          marker_info_file,
                          gamma_values,
                          save_dir,
                          n_parallel_worker) {
    # Read data
    dt_full <- fread(cell_dat_file)
    dt_full$Sample <- as.character(dt_full$Sample)

    # TODO: remove me. For testing only
    # dt_full <- head(dt_full, n=200)

    # Select markers to use
    marker_info <- fread(marker_info_file)
    markers_to_use <- paste0(
        marker_info[marker_info$marker_class == "type"]$marker_name,
        "_asinh_cf5"
    )

    # How many times to run the method?
    n_rep <- 2

    tic.clearlog()

    for (gamma_value in gamma_values) {
        for (rep in seq(n_rep)) {
            message(paste("Running supercell gamma", gamma_value, ", rep", rep))

            BPPARAM <- MulticoreParam(workers = n_parallel_worker, tasks = length(unique(dt_full$Sample)))

            tic(paste0("supercell_gamma", gamma_value, "_rep", rep))
            res <- runSuperCellCyto(
                dt = dt_full,
                markers = markers_to_use,
                sample_colname = "Sample",
                cell_id_colname = "CellId",
                gam = gamma_value,
                k_knn = 5,
                BPPARAM = BPPARAM,
                load_balancing = TRUE
            )
            toc(log = TRUE, quiet = TRUE)

            rm(BPPARAM)
            gc()
        }

        # Write out the results
        fwrite(
            x = res$supercell_expression_matrix,
            file = paste0(save_dir, "/supercellExpMat_gamma", gamma_value, ".csv")
        )
        fwrite(
            x = res$supercell_cell_map,
            file = paste0(save_dir, "/supercellCellMap_gamma", gamma_value, ".csv")
        )
        saveRDS(
            object = res$supercell_object,
            file = paste0(save_dir, "/supercellObj_gamma", gamma_value, ".rds")
        )
    }

    # Export the supercell durations
    writeLines(
        text = unlist(tic.log(format = TRUE)),
        con = paste0(save_dir, "/supercell_runtime.txt")
    )

    tic.clearlog()
}
