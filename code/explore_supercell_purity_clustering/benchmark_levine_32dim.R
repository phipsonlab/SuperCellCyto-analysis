#' ---- Benchmark supercell on Levine 32dim data ----
#' Essentially:
#' 1. Run supercell on various gamma values (5-50 increment by 5).
#' 2. Cluster supercell using Louvain with k values set to 10-30 incremented by 5.
#' 3. Cluster supercell using FlowSOM.
#' 4. Cluster all cells using both Louvain and FlowSOM. Note, the former has to be
#' ran on a HPC system as it requires enormous amount of RAM.
#' 5. Calculate ARI between the clustering of supercell and ground truth,
#' and between the clustering of supercells and all cells.
#' 6. Calculate the purity of supercells.
#'
#' Before you begin, please open the SuperCellCyto_analysis.Rproj file.
#' This is so your working directory is set to the project root.

library(HDCytoData)
library(data.table)
library(SuperCellCyto)
library(BiocParallel)
library(tictoc)
library(bluster)
library(Spectre, include.only=c("run.flowsom")) # only need the run flowsom function
library(here)
library(SuperCell) # for calculating purity

# ---- 0.1 Setup functions ----
# The functions to do the benchmarking are all stored here:
benchmark_func_dir <- here("code", "explore_supercell_purity_clustering", "functions")
for (f in list.files(benchmark_func_dir)) {
    source(here(benchmark_func_dir, f))
}

# ---- 0.2 Setup some global configs ----
gamma_values <- seq(5, 50, 5)

set.seed(42)

# ---- 0.3 Setup required directories ----

# Setup all the directories.
# Main directories
directories <- list(
    data_savedir = here("data", "explore_supercell_purity_clustering", "levine_32dim"),
    res_savedir = here("output", "explore_supercell_purity_clustering", format(Sys.Date(), "%Y%m%d"), "levine_32dim")
)

# Supercell directories and all cells clustering directories
directories <- append(x = directories, values = list(
    supercell_res_savedir = here(directories$res_savedir, "supercell_runs"),
    flowsom_allcells_savedir = here(directories$res_savedir, "flowsom_allcells"),
    eval_savedir = here(directories$res_savedir, "evaluation")
))


# Supercell louvain clustering directories
directories <- append(
    x = directories,
    values = get_supercell_clustering_res_directories(
        root_dir = directories$res_savedir,
        gamma_values = gamma_values,
        clust_algo = "louvain"
    )
)

# Supercell flowsom clustering directories
directories <- append(
    x = directories,
    values = get_supercell_clustering_res_directories(
        root_dir = directories$res_savedir,
        gamma_values = gamma_values,
        clust_algo = "flowsom"
    )
)

# Create them all!
for (d in directories) {
    dir.create(d, recursive = TRUE)
}

# ---- 0.4 Prepare raw data ----
prepare_data_levine(directories$data_savedir)


# ---- 1. Run supercell ----
run_supercell(
    cell_dat_file = here(directories$data_savedir, "levine_32_asinh.csv"),
    marker_info_file = here(directories$data_savedir, "levine_32_asinh_markers_info.csv"),
    gamma_values = gamma_values,
    save_dir = directories$supercell_res_savedir,
    n_parallel_worker = parallel::detectCores() - 1
)

# ---- 2. Cluster supercell with Louvain ----

# Setup configs
marker_info <- fread(here(directories$data_savedir, "levine_32_asinh_markers_info.csv"))
markers_for_clustering <- paste0(
    marker_info[marker_info$marker_class == 'type']$marker_name,
    "_asinh_cf5"
)

k_vals <- seq(10, 30, 5)

for (gamma_value in gamma_values) {

    message(paste("Louvain clustering supercell gamma", gamma_value))

    cluster_supercell_louvain(
        supercell_exp_mat_file = here(
            directories$supercell_res_savedir,
            paste0("supercellExpMat_gamma", gamma_value, ".csv")
        ),
        markers_to_use = markers_for_clustering,
        k_vals = k_vals,
        seed = 42,
        save_dir = directories[[paste0("supercell_gamma", gamma_value, "_louvain_savedir")]]
    )
}

# ---- 3. Cluster supercell with FlowSOM ----

# Setup configs
marker_info <- fread(here(directories$data_savedir, "levine_32_asinh_markers_info.csv"))
markers_for_clustering <- paste0(
    marker_info[marker_info$marker_class == 'type']$marker_name,
    "_asinh_cf5"
)
metaclusters <- seq(15, 30, by = 5)
grid_size <- seq(10, 14, by = 1)

fsom_config <- as.list(expand.grid(metacluster = metaclusters, grid_size = grid_size))

for (gamma_value in gamma_values) {

    message(paste("FlowSOM clustering supercell gamma", gamma_value))

    cluster_supercell_flowsom(
        supercell_exp_mat_file = here(
            directories$supercell_res_savedir,
            paste0("supercellExpMat_gamma", gamma_value, ".csv")
        ),
        markers_to_use = markers_for_clustering,
        clust_configs = fsom_config,
        seed = 42,
        save_dir = directories[[paste0("supercell_gamma", gamma_value, "_flowsom_savedir")]]
    )
}

# ---- 4. Cluster all cells with FlowSOM ----
# Refer to the nextflow pipeline for the Louvain clustering.

# Note, these configs are the same as the configs used to cluster the supercells
# above. We're restating them here so they are explicit.
marker_info <- fread(here(directories$data_savedir, "levine_32_asinh_markers_info.csv"))
markers_for_clustering <- paste0(
    marker_info[marker_info$marker_class == 'type']$marker_name,
    "_asinh_cf5"
)
metaclusters <- seq(15, 30, by = 5)
grid_size <- seq(10, 14, by = 1)

fsom_config <- as.list(expand.grid(metacluster = metaclusters, grid_size = grid_size))

cluster_allcell_flowsom(
    cell_dat_file = here(directories$data_savedir, "levine_32_asinh.csv"),
    markers_to_use = markers_for_clustering,
    clust_configs = fsom_config,
    save_dir = directories$flowsom_allcells_savedir
)

# ---- 5. Evaluate supercells ----
# Note, these configs are the same as the configs used to cluster the supercells
# above. We're restating them here so they are explicit.
metaclusters <- seq(15, 30, by = 5)
grid_size <- seq(10, 14, by = 1)
fsom_config <- as.list(expand.grid(metacluster = metaclusters, grid_size = grid_size))
k_vals <- seq(10, 30, 5)

fwrite(
    x = calc_ari_supercell_allcell_flowsom(
        ground_truth_file = here(directories$data_savedir, "levine_32_asinh.csv"),
        fsom_config = fsom_config,
        gamma_values = gamma_values,
        supercell_res_dir = directories$supercell_res_savedir,
        clust_supercell_res_dir = here(directories$res_savedir, "flowsom_supercell_runs"),
        clust_allcell_res_dir = here("output", "explore_supercell_purity_clustering", "20230509", "levine_32dim", "flowsom_allcells")
    ),
    file = here(directories$eval_savedir, "flowsom_ari_vs_all.csv")
)

fwrite(
    x = calc_ari_supercell_truth_flowsom(
        ground_truth_file = here(directories$data_savedir, "levine_32_asinh.csv"),
        fsom_config = fsom_config,
        gamma_values = gamma_values,
        supercell_res_dir = directories$supercell_res_savedir,
        clust_res_dir = here(directories$res_savedir, "flowsom_supercell_runs")
    ),
    file = here(directories$eval_savedir, "flowsom_ari_vs_truth.csv")
)

fwrite(
    x = calc_ari_supercell_allcell_louvain(
        ground_truth_file = here(directories$data_savedir, "levine_32_asinh.csv"),
        k_values = k_vals,
        gamma_values = gamma_values,
        supercell_res_dir = directories$supercell_res_savedir,
        clust_supercell_res_dir = here(directories$res_savedir, "louvain_supercell_runs"),
        clust_allcell_res_dir = here(directories$res_savedir, "all_louvain")
    ),
    file = here(directories$eval_savedir, "louvain_ari_vs_all.csv")
)

fwrite(
    x = calc_ari_supercell_truth_louvain(
        ground_truth_file = here(directories$data_savedir, "levine_32_asinh.csv"),
        k_values = k_vals,
        gamma_values = gamma_values,
        supercell_res_dir = directories$supercell_res_savedir,
        clust_res_dir = here(directories$res_savedir, "louvain_supercell_runs")
    ),
    file = here(directories$eval_savedir, "louvain_ari_vs_truth.csv")
)

fwrite(
    x = calc_supercell_purity(
        ground_truth_file = here(directories$data_savedir, "levine_32_asinh.csv"),
        gamma_values = gamma_values,
        supercell_res_dir = directories$supercell_res_savedir
    ),
    file = here(directories$eval_savedir, "supercell_purities.csv")
)


