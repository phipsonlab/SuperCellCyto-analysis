#' Prepare Levine_32dim data
#'
#' Download Levine_32dim data from HDCytoData package
#' and convert the cells, markers' info, and sample info
#' into a data.frame, and save them out as csv files.
#'
#' Require: HDCytoData and data.table packages.
#'
#' @param savedir character. Directory where to store the generated csv files.
#'
#' @author Givanna Putri
#'
prepare_data_levine <- function(savedir) {
    se <- Levine_32dim_SE()

    # get marker info
    marker_info <- data.table(as.data.frame(colData(se)))

    # actual data
    dat <- data.table(assays(se)$exprs)

    #  need to transform the data using asinh
    markers <- marker_info[marker_info$marker_class != "none"]$marker_name
    cf <- 5
    dat_trans <- asinh(dat[, markers, with = F] / cf)
    names(dat_trans) <- paste0(names(dat_trans), "_asinh_cf", cf)
    dat <- cbind(dat, dat_trans)

    dat$Gated_Population <- rowData(se)$population_id
    dat$Sample <- rowData(se)$patient_id

    # Assign cell ID
    dat$CellId <- paste0("Cell_", seq_len(nrow(dat)))

    # Export data
    fwrite(marker_info, paste0(savedir, "/levine_32_asinh_markers_info.csv"))
    fwrite(dat, paste0(savedir, "/levine_32_asinh.csv"))
}

#' Prepare Levine_32dim data
#'
#' Download Levine_32dim data from HDCytoData package
#' and convert the cells, markers' info, and sample info
#' into a data.frame, and save them out as csv files.
#'
#' Require: HDCytoData and data.table packages.
#'
#' @param savedir character. Directory where to store the generated csv files.
#'
#' @author Givanna Putri
#'
prepare_data_samusik <- function(savedir) {
    se <- Samusik_all_SE()

    # get marker info
    marker_info <- data.table(as.data.frame(colData(se)))

    # actual data
    dat <- data.table(assays(se)$exprs)

    #  need to transform the data using asinh
    markers <- marker_info[marker_info$marker_class != "none"]$marker_name
    cf <- 5
    dat_trans <- asinh(dat[, markers, with = F] / cf)
    names(dat_trans) <- paste0(names(dat_trans), "_asinh_cf", cf)
    dat <- cbind(dat, dat_trans)

    dat$Gated_Population <- rowData(se)$population_id
    dat$Sample <- rowData(se)$sample_id

    # Assign cell ID
    dat$CellId <- paste0("Cell_", seq_len(nrow(dat)))

    # Export data
    fwrite(marker_info, paste0(savedir, "/samusik_all_asinh_markers_info.csv"))
    fwrite(dat, paste0(savedir, "/samusik_all_asinh.csv"))
}
