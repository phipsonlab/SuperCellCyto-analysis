#' Cluster supercells using Louvain
#'
#' @param supercell_exp_mat_file
#' @param markers_to_use
#' @param k_vals
#' @param save_dir
#' @param seed
#'
#' @return
#' @export
#'
#' @examples
cluster_supercell_louvain <- function(supercell_exp_mat_file,
                                      markers_to_use,
                                      k_vals,
                                      save_dir,
                                      seed = 42) {
    exp_mt <- fread(supercell_exp_mat_file)

    exp_mt_toCluster <- exp_mt[, markers_to_use, with = F]

    # Subset for testing
    # exp_mt_toCluster <- rbind(head(exp_mt_toCluster, 100), tail(exp_mt_toCluster, 200))

    # How many times to run the clustering?
    n_rep <- 2

    tic.clearlog()

    for (k in k_vals) {
        for (rep in seq(n_rep)) {
            message(paste("Louvain clustering supercell k", k, ", rep", rep))

            set.seed(seed)

            tic(paste0("louvain_supercell_k", k, "_rep", rep))

            clusters <- clusterRows(
                x = exp_mt_toCluster,
                BLUSPARAM = NNGraphParam(
                    k = k,
                    cluster.fun = "louvain",
                    BPPARAM = MulticoreParam(RNGseed = seed)
                )
            )

            toc(log = TRUE, quiet = TRUE)

            clust_only <- data.table(louvain_cluster = clusters)

            fwrite(
                x = clust_only,
                file = here(save_dir, paste0("louvain_k", k, "_rep", rep, "_supercell_clusters.csv"))
            )
        }
    }

    # Export the clustering durations
    writeLines(
        text = unlist(tic.log(format = TRUE)),
        con = paste0(save_dir, "/louvain_supercell_runtime.txt")
    )

    tic.clearlog()
}


#' Cluster supercell using FlowSOM
#'
#' Cluster supercell using the FlowSOM algorithm.
#' Note, please run `set.seed` before running this function to ensure reproducibility.
#'
#' @return
#' @export
#'
#' @examples
cluster_supercell_flowsom <- function(supercell_exp_mat_file,
                                      markers_to_use,
                                      clust_configs,
                                      save_dir,
                                      seed = 42) {


    exp_mt <- fread(supercell_exp_mat_file)

    # Subset for testing
    # exp_mt_toCluster <- rbind(head(exp_mt_toCluster, 100), tail(exp_mt_toCluster, 200))

    # How many times to run the clustering?
    n_rep <- 2

    tic.clearlog()

    for (i in seq(length(clust_configs$metacluster))) {
        metacluster <- clust_configs$metacluster[i]
        grid_size <- clust_configs$grid_size[i]

        for (rep in seq(n_rep)) {

            message(paste("Flowsom clustering supercell meta", metacluster, "grid", grid_size, "rep", rep))

            tic(paste0("flowsom_supercell_meta", metacluster, "_grid", grid_size, "_rep", rep))

            clusters <- run.flowsom(
                dat=exp_mt,
                use.cols=markers_to_use,
                meta.k=metacluster,
                xdim=grid_size,
                ydim=grid_size,
                meta.seed=seed,
                clust.seed=seed
            )[, c("FlowSOM_cluster", "FlowSOM_metacluster"), with=F]

            toc(log = TRUE, quiet = TRUE)

        }
        fwrite(
            x = clusters,
            file = here(
                save_dir,
                paste0("flowsom_meta", metacluster, "_grid", grid_size, "_supercell_clusters.csv")
            )
        )
    }


    # Export the clustering durations
    writeLines(
        text = unlist(tic.log(format = TRUE)),
        con = here(save_dir, "flowsom_supercell_runtime.txt")
    )

    tic.clearlog()
}

#' Cluster all cells using FlowSOM
#'
#' @param cell_dat_file
#' @param markers_to_use
#' @param clust_configs
#' @param save_dir
#' @param seed
#'
#' @return
#' @export
#'
#' @examples
cluster_allcell_flowsom <- function(cell_dat_file,
                                      markers_to_use,
                                      clust_configs,
                                      save_dir,
                                      seed = 42) {


    exp_mt <- fread(cell_dat_file)

    # Subset for testing
    # exp_mt_toCluster <- rbind(head(exp_mt_toCluster, 100), tail(exp_mt_toCluster, 200))

    # How many times to run the clustering?
    n_rep <- 2

    tic.clearlog()

    for (i in seq(length(clust_configs$metacluster))) {
        metacluster <- clust_configs$metacluster[i]
        grid_size <- clust_configs$grid_size[i]

        for (rep in seq(n_rep)) {

            message(paste("Flowsom clustering all cells meta", metacluster, "grid", grid_size, "rep", rep))

            tic(paste0("flowsom_allcell_meta", metacluster, "_grid", grid_size, "_rep", rep))

            clusters <- run.flowsom(
                dat=exp_mt,
                use.cols=markers_to_use,
                meta.k=metacluster,
                xdim=grid_size,
                ydim=grid_size,
                meta.seed=seed,
                clust.seed=seed
            )[, c("FlowSOM_cluster", "FlowSOM_metacluster"), with=F]

            toc(log = TRUE, quiet = TRUE)

        }
        fwrite(
            x = clusters,
            file = here(
                save_dir,
                paste0("flowsom_meta", metacluster, "_grid", grid_size, "_allcell_clusters.csv")
            )
        )
    }


    # Export the clustering durations
    writeLines(
        text = unlist(tic.log(format = TRUE)),
        con = here(save_dir, "flowsom_allcell_runtime.txt")
    )

    tic.clearlog()
}
