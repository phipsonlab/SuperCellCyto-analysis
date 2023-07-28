#' Script to run supercell on Oetjen B cell panel data
#' and to cluster the supercells using Louvain clustering.

library(data.table)
library(SuperCellCyto)
library(BiocParallel)
library(tictoc)
library(here)
library(bluster)

# ---- 1. Prepare data ----
cell_dat <- fread(here("output", "oetjen_b_cell_panel",
                       "20230404_cytoexplorer_gating",
                       "2023-04-04_single_live_cells.csv"))

out_dir <- here("output", "oetjen_b_cell_panel", format(Sys.Date(), "%Y%m%d"))

dir.create(out_dir, recursive = TRUE)

# Asinh transformation with co-factor 150
markers <- c(
  "CD21", "CD40", "CD10", "CD138", "CD80", "CD27", "CD19", "CD45",
  "PD_1", "CD20", "CD38", "CD86"
)

cf <- 150
asinh_colnames <- paste0(markers, "_asinh_cf", cf)
cell_dat[, (asinh_colnames) := lapply(.SD, function(x) asinh(x / cf)), .SDcols = markers]

# Assign the cell id
cell_dat$CellId <- paste0("Cell_", seq_len(nrow(cell_dat)))

# For latter quick access
fwrite(cell_dat, here(out_dir, "cell_dat_asinh.csv"))

# ---- 2. Run supercell ----
n_rep <- 2

# Run Supercell

tic.clearlog()

# We will set gamma to 20 for now.
for (x in seq_len(n_rep)) {
  message(paste0("supercell_run_", x))

  # Set up the parallel execution object
  BPPARAM <- MulticoreParam(workers = parallel::detectCores(), tasks = length(unique(cell_dat$sample)))

  tic(paste0("supercell_run_", x))
  supercell_obj <- runSuperCellCyto(
    dt = cell_dat,
    markers = paste0(markers, "_asinh_cf", cf),
    sample_colname = "sample",
    cell_id_colname = "CellId",
    gam = 20,
    BPPARAM = BPPARAM,
    load_balancing = TRUE
  )
  toc(log = TRUE, quiet = TRUE)

  rm(BPPARAM)
  gc()
}

# Export the clustering durations
writeLines(
  text = unlist(tic.log(format = TRUE)),
  con = here(out_dir, "supercell_runtime.txt")
)

tic.clearlog()

fwrite(
  x = supercell_obj$supercell_expression_matrix,
  file = here(out_dir, "supercellExpMat.csv")
)
fwrite(
  x = supercell_obj$supercell_cell_map,
  file = here(out_dir, "supercellCellMap.csv")
)
saveRDS(
  object = supercell_obj$supercell_object,
  file = here(out_dir, "supercellObj.rds")
)

# ---- 3. Louvain clustering ----
# Use only cell type markers (i.e. the one used in gating)
cell_type_markers <- paste0(
  c("CD19", "CD45", "CD10", "CD20", "CD27", "CD21", "CD38", "CD138"),
  "_asinh_cf",
  cf
)

exp_mt_toCluster <- supercell_obj$supercell_expression_matrix[, cell_type_markers, with = FALSE]

tic.clearlog()

for (rep in seq(n_rep)) {
  message(paste("Louvain clustering supercell rep", rep))

  set.seed(42)

  tic(paste0("louvain_supercell_k3_rep", rep))

  clusters <- clusterRows(
    x = exp_mt_toCluster,
    BLUSPARAM = NNGraphParam(
      k = 3,
      cluster.fun = "louvain",
      BPPARAM = MulticoreParam(RNGseed = 42)
    )
  )

  toc(log = TRUE, quiet = TRUE)

  clust_only <- data.table(louvain_cluster = clusters)

  fwrite(
    x = clust_only,
    file = here(out_dir, paste0("louvain_k3_rep", rep, "_supercell_clusters.csv"))
  )
}

# Export the clustering durations
writeLines(
  text = unlist(tic.log(format = TRUE)),
  con = paste0(here(out_dir, "louvain_supercell_runtime.txt"))
)

tic.clearlog()

# ---- 4. Run UMAP ----
supercell_mat <- supercell_obj$supercell_expression_matrix

clusters <- fread(here(out_dir, paste0("louvain_k3_rep1_supercell_clusters.csv")))

supercell_mat[, louvain_k3 := clusters]

cell_type_markers <- paste0(
    c("CD19", "CD45", "CD10", "CD20", "CD27", "CD21", "CD38", "CD138"),
    "_asinh_cf150"
)

supercell_mat <- Spectre::run.fast.umap(
    supercell_mat,
    use.cols = cell_type_markers
)

fwrite(supercell_mat, here(out_dir, "umap_sce_supercell.csv"))
