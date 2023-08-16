#' Script to run supercell on Oetjen B cell panel data
#' and to cluster the supercells using Louvain clustering.

library(data.table)
library(SuperCellCyto)
library(BiocParallel)
library(tictoc)
library(here)
library(bluster)

# ---- 1. Prepare data ----
cell_dat <- fread("output/oetjen_b_cell_panel/20230511/cell_dat_asinh.csv")


# ---- 3. Louvain clustering ----
# Use only cell type markers (i.e. the one used in gating)
cell_type_markers <- paste0(
  c("CD19", "CD45", "CD10", "CD20", "CD27", "CD21", "CD38", "CD138"),
  "_asinh_cf150"
)

exp_mt_toCluster <- cell_dat[, cell_type_markers, with = FALSE]

tic.clearlog()

tic("louvain_singlecell")

clusters <- clusterRows(
  x = exp_mt_toCluster,
  BLUSPARAM = NNGraphParam(
    k = 3,
    cluster.fun = "louvain",
    BPPARAM = MulticoreParam(RNGseed = 42)
  )
)

toc(log = TRUE, quiet = TRUE)

# Export the clustering durations
writeLines(
  text = unlist(tic.log(format = TRUE)),
  con = paste0("output/oetjen_b_cell_panel/louvain_clust_singlecell_runtime.txt")
)


