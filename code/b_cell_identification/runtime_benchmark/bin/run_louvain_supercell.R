#!/usr/bin/env Rscript

library(data.table)
library(BiocParallel)
library(tictoc)
library(bluster)

# ---- 1. Prepare data ----

args <- commandArgs(trailingOnly=TRUE)

cell_dat <- fread(args[1])
k <- as.numeric(args[2])
gam <- args[3]

# Use only cell type markers (i.e. the one used in gating)
cell_type_markers <- paste0(
    c("CD19", "CD45", "CD10", "CD20", "CD27", "CD21", "CD38", "CD138"),
    "_asinh_cf150"
)

BPPARAM <- MulticoreParam(RNGseed = 42)

exp_mt_toCluster <- cell_dat[, cell_type_markers, with = FALSE]

tic.clearlog()

tic(paste0("louvain_supercell_gamma", gam, "_k", k))

clusters <- clusterRows(
    x = exp_mt_toCluster,
    BLUSPARAM = NNGraphParam(
        k = k,
        cluster.fun = "louvain",
        BPPARAM = BPPARAM
    )
)

toc(log = TRUE, quiet = TRUE)

writeLines(
    text = unlist(tic.log(format = TRUE)),
    con = paste0(
        "duration_louvain_supercell_gamma", gam, "_k", k, ".txt"
    )
)

cluster_dt <- data.table(cluster=clusters)
fwrite(cluster_dt, paste0("cluster_louvain_supercell_gamma", gam, "_k", k, ".csv"))





