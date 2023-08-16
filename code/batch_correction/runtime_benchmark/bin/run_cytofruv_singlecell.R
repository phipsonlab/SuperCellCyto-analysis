#!/usr/bin/env Rscript

library(here)
library(CytofRUV)
library(CATALYST)
library(data.table)
library(rsvd) # For cytofruv
library(tictoc)

source("/stornext/Bioinf/data/lab_phipson/givanna/github/SuperCellCyto-analysis/code/batch_correction/cytofRUV_functions.R")

args <- commandArgs(trailingOnly=TRUE)

dat <- fread(args[1])
panel_info <- fread(args[2])
meta <- as.numeric(args[3])
k <- as.numeric(args[4])
setting <- as.numeric(args[5])

# Discard the untransformed markers to save RAM
cols_to_keep <- names(dat)[c(32:68)]
dat <- dat[, cols_to_keep, with = FALSE]

markers <- names(dat)[7:37]

set.seed(42)

# ---- Run CytofRUV ----
# Create sce object

supercell_exp_mat_transposed <- as.matrix(t(dat[, markers, with = FALSE]))

rownames(supercell_exp_mat_transposed) <- markers
colnames(supercell_exp_mat_transposed) <- dat$CellID

# create a new column which concatenates the reporter and the antigen
panel_info[, reporter_marker := factor(paste(fcs_colname, antigen, sep = "_"), levels = gsub("_asinh_cf5", "", markers))]
panel_info <- DataFrame(panel_info[order(reporter_marker)])
rownames(panel_info) <- paste0(panel_info$reporter_marker, "_asinh_cf5")

sample_meta <- dat[, c("sample_id", "condition", "patient_id", "batch")]

# cluster data
celltype_markers <- as.vector(rownames(panel_info[panel_info$marker_class == 'type', ]))

sce <- SingleCellExperiment(
    list(exprs = supercell_exp_mat_transposed),
    colData = sample_meta,
    rowData = panel_info
)
tic.clearlog()

tic(paste("Run cytofruv setting", setting))
sce <- cluster_data(sce, 42, markers_to_use = celltype_markers, meta)

# saveRDS(sce, here(outdir, "sce_object.rds"))

# Use the CLL2 and the HCL1 samples as pseudo replicates.

# Extract clustered data as data.frame
supercell_clust <- data.frame(
    sample = sce$sample_id,
    cluster = cluster_ids(sce, paste0("meta", meta)),
    t(assay(sce, "exprs"))
)

norm_supercell <- run_RUVIII(
    data = supercell_clust,
    norm_clusters = list(seq(1, meta), seq(1, meta)),
    k = k,
    rep_samples = list(c("CLL2_B1", "CLL2_B2"), c("HC1_B1", "HC1_B2"))
)
toc(log = TRUE, quiet = TRUE)

# Some wrangling to store the result
norm_supercell_df <- as.data.table(norm_supercell)
norm_supercell_df$CellId <- rownames(norm_supercell)
norm_supercell_df$sample_id <- supercell_clust$sample

fwrite(norm_supercell_df, paste0("dat_cytofruv_singlecell_setting",setting,".csv"))

writeLines(
    text = unlist(tic.log(format = TRUE)),
    con = paste0("cytofruv_singlecell_setting",setting,".txt")
)