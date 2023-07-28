#!/usr/bin/env Rscript

# Run Louvain clustering on all cells

library(data.table)
library(bluster)
library(BiocParallel)
library(tictoc)

args <- commandArgs(trailingOnly=TRUE)

# ---- Parse configs ----
exp_mt <- fread(args[1])
markers <- fread(args[2])
k <- as.numeric(args[3])
out_dir <- args[4]
n_cores <- as.numeric(args[5])
n_reps <- as.numeric(args[6])

# ---- For testing ----
# setwd("/vast/scratch/users/putri.g/supercell/samusik_evaluation/nextflow_int_files/")
# exp_mt <- fread("dataset/samusik_all_asinh.csv")
# markers <- fread("dataset/samusik_all_asinh_markers_info.csv")
# n_metacluster <- 20
# grid_size <- 10

# Subset for testing
# exp_mt <- rbind(head(exp_mt, 100), tail(exp_mt, 200))

# ---- Setting up ----
set.seed(42)
markers_to_use <- paste0(
    markers[marker_class == "type"]$marker_name, "_asinh_cf5"
)
exp_mt_toCluster <- exp_mt[, markers_to_use, with=F]
# Subset for testing
# exp_mt_toCluster <- rbind(head(exp_mt_toCluster, 100), tail(exp_mt_toCluster, 200))

tic.clearlog()
# ---- Run Louvain ----
for (rep in seq(n_reps)) {

    BPPARAM <- MulticoreParam(workers = n_cores, RNGseed=42)
    # Bluster params. By default, this will be a shared nearest neighbour
    BLUSPARAM <- NNGraphParam(k=k, cluster.fun="louvain", BPPARAM=BPPARAM)

    tic(paste0("louvain_allcells_k", k, "_rep", rep))
    clusters <- clusterRows(x=exp_mt_toCluster, BLUSPARAM=BLUSPARAM)
    toc(log = TRUE, quiet = TRUE)

    rm(BPPARAM)
    gc()
}


# ---- Write out result ----

clust_only <- data.table(louvain_cluster=clusters)
fwrite(clust_only, paste0(out_dir, "/louvain_clusters_allcells.csv"))

# Export the clustering durations
writeLines(
    text = unlist(tic.log(format = TRUE)),
    con = paste0(out_dir, "/louvain_supercell_runtime.txt")
)


