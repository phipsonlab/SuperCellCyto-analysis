calc_ari_supercell_truth_louvain <- function(ground_truth_file,
                                     k_values,
                                     gamma_values,
                                     supercell_res_dir,
                                     clust_res_dir) {

    # For testing only
    # k_values <- seq(10, 30, 5)
    # gamma_values <- seq(5, 50, 5)
    # ground_truth_file <- "data/clustering_benchmarking/v2_pure_R_pipeline/levine_32dim/levine_32_asinh.csv"
    # supercell_res_dir <- here("output", "clustering_benchmarking", "20230511", "levine_32dim", "supercell_runs")
    # clust_res_dir <- here("output", "clustering_benchmarking", "20230511", "levine_32dim", "louvain_supercell_runs")

    ground_truth <- fread(ground_truth_file)[, c("CellId", "Gated_Population")]
    ground_truth <- ground_truth[ground_truth$Gated_Population != 'unassigned']

    ari_per_gamma <- lapply(gamma_values, function(gam) {
        ari_per_k <- sapply(k_values, function(k) {
            supercell_clusters <- data.frame(
                clusters = fread(
                    here(clust_res_dir,
                         paste0("gamma_", gam),
                         paste0("louvain_k", k, "_rep1_supercell_clusters.csv"))
                )$louvain_cluster,
                SuperCellID = fread(
                    here(supercell_res_dir, paste0("supercellExpMat_gamma", gam, ".csv"))
                )$SuperCellId
            )

            supercell_cell_map <- fread(here(supercell_res_dir, paste0("supercellCellMap_gamma", gam, ".csv")))

            # TODO: convert to merge.data.table for clarity
            all_data <- ground_truth[supercell_cell_map, on=c("CellId"), nomatch = 0]
            all_data <- all_data[supercell_clusters, on=c("SuperCellID"), nomatch = 0]

            return(ARI(all_data$clusters, all_data$Gated_Population))
        })
        return(data.frame(k = k_values, ari = ari_per_k, gamma = gam))
    })
    ari_per_gamma <- rbindlist(ari_per_gamma)
    return(ari_per_gamma)
}


calc_ari_supercell_allcell_louvain <- function(ground_truth_file,
                                               k_values,
                                               gamma_values,
                                               supercell_res_dir,
                                               clust_supercell_res_dir,
                                               clust_allcell_res_dir) {

    # For testing only
    # k_values <- seq(10, 30, 5)
    # gamma_values <- seq(5, 50, 5)
    # ground_truth_file <- "data/clustering_benchmarking/v2_pure_R_pipeline/levine_32dim/levine_32_asinh.csv"
    # supercell_res_dir <- here("output", "clustering_benchmarking", "20230511", "levine_32dim", "supercell_runs")
    # clust_supercell_res_dir <- here("output", "clustering_benchmarking", "20230511", "levine_32dim", "louvain_supercell_runs")
    # clust_allcell_res_dir <- here("output", "clustering_benchmarking", "20230511", "levine_32dim", "all_louvain")

    cell_ids <- fread(ground_truth_file)[, c("CellId")]

    ari_per_k <- lapply(k_values, function(k) {
        allcell_clusters <- fread(here(clust_allcell_res_dir,
                                       paste0("k", k),
                                       "louvain_clusters_allcells.csv"))
        setnames(allcell_clusters, "louvain_cluster", "allcell_clust")
        allcell_clusters$CellId <- cell_ids


        ari_per_gamma <- sapply(gamma_values, function(gam) {
            supercell_clusters <- data.frame(
                supercell_clust = fread(
                    here(clust_supercell_res_dir,
                         paste0("gamma_", gam),
                         paste0("louvain_k", k, "_rep1_supercell_clusters.csv"))
                )$louvain_cluster,
                SuperCellID = fread(
                    here(supercell_res_dir, paste0("supercellExpMat_gamma", gam, ".csv"))
                )$SuperCellId
            )

            supercell_cell_map <- fread(here(supercell_res_dir, paste0("supercellCellMap_gamma", gam, ".csv")))

            all_data <- merge.data.table(
                x = merge.data.table(
                    x = allcell_clusters,
                    y = supercell_cell_map,
                    by = "CellId"
                ),
                y = supercell_clusters,
                by = "SuperCellID"
            )
            return(ARI(all_data$allcell_clust, all_data$supercell_clust))
        })
        return(data.frame(k = k, ari = ari_per_gamma, gamma = gamma_values))
    })
    ari_per_k <- rbindlist(ari_per_k)

    return(ari_per_k)


}


calc_ari_supercell_truth_flowsom <- function(ground_truth_file,
                                             fsom_config,
                                             gamma_values,
                                             supercell_res_dir,
                                             clust_res_dir) {

    # For testing only
    # metaclusters <- seq(15, 30, by = 5)
    # grid_size <- seq(10, 14, by = 1)
    # fsom_config <- as.list(expand.grid(metacluster = metaclusters, grid_size = grid_size))
    # gamma_values <- seq(5, 50, 5)
    # ground_truth_file <- "data/clustering_benchmarking/v2_pure_R_pipeline/levine_32dim/levine_32_asinh.csv"
    # supercell_res_dir <- here("output", "clustering_benchmarking", "20230511", "levine_32dim", "supercell_runs")
    # clust_res_dir <- here("output", "clustering_benchmarking", "20230511", "levine_32dim", "flowsom_supercell_runs")

    ground_truth <- fread(ground_truth_file)[, c("CellId", "Gated_Population")]
    ground_truth <- ground_truth[ground_truth$Gated_Population != 'unassigned']

    ari_per_gamma <- lapply(gamma_values, function(gam) {
        ari_per_conf <- sapply(seq_len(length(fsom_config$metacluster)), function(i) {

            supercell_clusters <- data.frame(
                clusters = fread(
                    here(clust_res_dir,
                         paste0("gamma_", gam),
                         paste0(
                             "flowsom_meta", fsom_config$metacluster[i],
                             "_grid", fsom_config$grid_size[i],
                             "_supercell_clusters.csv"))
                )$FlowSOM_metacluster,
                SuperCellID = fread(
                    here(supercell_res_dir, paste0("supercellExpMat_gamma", gam, ".csv"))
                )$SuperCellId
            )

            supercell_cell_map <- fread(here(supercell_res_dir, paste0("supercellCellMap_gamma", gam, ".csv")))

            # TODO: convert to merge.data.table for clarity
            all_data <- ground_truth[supercell_cell_map, on=c("CellId"), nomatch = 0]
            all_data <- all_data[supercell_clusters, on=c("SuperCellID"), nomatch = 0]

            return(ARI(all_data$clusters, all_data$Gated_Population))
        })
        return(data.frame(
            n_metaclust = fsom_config$metacluster,
            grid_size = fsom_config$grid_size,
            ari = ari_per_conf,
            gamma = gam))
    })
    ari_per_gamma <- rbindlist(ari_per_gamma)
    return(ari_per_gamma)
}


calc_ari_supercell_allcell_flowsom <- function(ground_truth_file,
                                               fsom_config,
                                               gamma_values,
                                               supercell_res_dir,
                                               clust_supercell_res_dir,
                                               clust_allcell_res_dir) {

    # For testing only
    # metaclusters <- seq(15, 30, by = 5)
    # grid_size <- seq(10, 14, by = 1)
    # fsom_config <- as.list(expand.grid(metacluster = metaclusters, grid_size = grid_size))
    # gamma_values <- seq(5, 50, 5)
    # ground_truth_file <- "data/clustering_benchmarking/v2_pure_R_pipeline/levine_32dim/levine_32_asinh.csv"
    # supercell_res_dir <- here("output", "clustering_benchmarking", "20230511", "levine_32dim", "supercell_runs")
    # clust_supercell_res_dir <- here("output", "clustering_benchmarking", "20230511", "levine_32dim", "flowsom_supercell_runs")
    # clust_allcell_res_dir <- here("output", "clustering_benchmarking", "20230509", "levine_32dim", "flowsom_allcells")

    cell_ids <- fread(ground_truth_file)[, c("CellId")]

    ari_per_conf <- lapply(seq_len(length(fsom_config$metacluster)), function(i) {
        allcell_clusters <- fread(here(clust_allcell_res_dir,
                                       paste0("flowsom_meta", fsom_config$metacluster[i],
                                              "_grid", fsom_config$grid_size[i],
                                              "_allcell_clusters.csv")))
        setnames(allcell_clusters, "FlowSOM_metacluster", "allcell_clust")
        allcell_clusters$CellId <- cell_ids
        allcell_clusters$FlowSOM_cluster <- NULL

        ari_per_gamma <- sapply(gamma_values, function(gam) {
            supercell_clusters <- data.frame(
                supercell_clust = fread(
                    here(clust_supercell_res_dir,
                         paste0("gamma_", gam),
                         paste0(
                             "flowsom_meta", fsom_config$metacluster[i],
                             "_grid", fsom_config$grid_size[i],
                             "_supercell_clusters.csv"))
                )$FlowSOM_metacluster,
                SuperCellID = fread(
                    here(supercell_res_dir, paste0("supercellExpMat_gamma", gam, ".csv"))
                )$SuperCellId
            )

            supercell_cell_map <- fread(here(supercell_res_dir, paste0("supercellCellMap_gamma", gam, ".csv")))

            all_data <- merge.data.table(
                x = merge.data.table(
                    x = allcell_clusters,
                    y = supercell_cell_map,
                    by = "CellId"
                ),
                y = supercell_clusters,
                by = "SuperCellID"
            )
            return(ARI(all_data$allcell_clust, all_data$supercell_clust))
        })
        return(data.frame(
            n_metaclust = fsom_config$metacluster[i],
            grid_size = fsom_config$grid_size[i],
            ari = ari_per_gamma,
            gamma = gamma_values))
    })
    ari_per_conf <- rbindlist(ari_per_conf)

    return(ari_per_conf)


}

calc_supercell_purity <- function(ground_truth_file,
                                  gamma_values,
                                  supercell_res_dir) {
    # For testing only
    # gamma_values <- seq(5, 50, 5)
    # ground_truth_file <- "data/clustering_benchmarking/v2_pure_R_pipeline/levine_32dim/levine_32_asinh.csv"
    # supercell_res_dir <- here("output", "clustering_benchmarking", "20230511", "levine_32dim", "supercell_runs")

    ground_truth <- fread(ground_truth_file)[, c("CellId", "Gated_Population")]
    ground_truth <- ground_truth[ground_truth$Gated_Population != 'unassigned']

    purity_per_gamma <- lapply(gamma_values, function(gam) {
        supercell_cell_map <- fread(here(supercell_res_dir, paste0("supercellCellMap_gamma", gam, ".csv")))

        all_data <- merge.data.table(
            x = supercell_cell_map,
            y = ground_truth,
            by = "CellId"
        )

        purity_scores <- data.frame(supercell_purity(all_data$Gated_Population, all_data$SuperCellID))
        names(purity_scores) <- "purity"
        purity_scores$SuperCellID <- rownames(purity_scores)
        rownames(purity_scores) <- NULL
        purity_scores <- purity_scores[, c("SuperCellID", "purity")]
        purity_scores$gamma <- gam

        return(purity_scores)
    })
    purity_per_gamma <- rbindlist(purity_per_gamma)
    return(purity_per_gamma)
}






