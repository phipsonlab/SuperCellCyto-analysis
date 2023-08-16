#!/usr/bin/env Rscript

library(CytofRUV)
library(CATALYST)
library(data.table)
library(rsvd) # For cytofruv
library(tictoc)

source("/stornext/Bioinf/data/lab_phipson/givanna/github/SuperCellCyto-analysis/code/batch_correction/cytofRUV_functions.R")

args <- commandArgs(trailingOnly=TRUE)

cell_dat <- fread(args[1])
panel_info <- fread(args[2])
meta <- as.numeric(args[3])
k <- as.numeric(args[4])
setting <- as.numeric(args[5])
sample_metadata <- fread(args[6])
gam <- args[7]


markers <- names(cell_dat)[c(1:31)]

set.seed(42)

# Create the marker (panel information)
# create a new column which concatenates the reporter and the antigen
panel_info[, reporter_marker := factor(paste(fcs_colname, antigen, sep = "_"), levels = gsub("_asinh_cf5", "", markers))]
panel_info <- DataFrame(panel_info[order(reporter_marker)])
rownames(panel_info) <- paste0(panel_info$reporter_marker, "_asinh_cf5")

# cluster data
celltype_markers <- as.vector(rownames(panel_info[panel_info$marker_class == 'type', ]))
# Create sce object
supercell_exp_mat_transposed <- as.matrix(t(cell_dat[, markers, with = FALSE]))
rownames(supercell_exp_mat_transposed) <- markers
colnames(supercell_exp_mat_transposed) <- cell_dat$SuperCellId

# Marry the supercell id and the metadata
supercell_metadata <- merge.data.table(
    x = cell_dat[, c("sample_id", "SuperCellId")],
    y = sample_metadata,
    by = "sample_id"
)

# Need to make sure the order is the same as the transposed supercell exp matrix.
supercell_metadata[, SuperCellId := factor(SuperCellId, levels = colnames(supercell_exp_mat_transposed))]
supercell_metadata <- DataFrame(supercell_metadata[order(SuperCellId)])
rownames(supercell_metadata) <- supercell_metadata$SuperCellId
    
sce <- SingleCellExperiment(
    list(exprs = supercell_exp_mat_transposed),
    colData = supercell_metadata,
    rowData = panel_info
)

tic.clearlog()
tic(paste("Run cytofruv supercell gam", gam, "setting", setting))
sce <- cluster_data(sce, 42, markers_to_use = celltype_markers, meta)

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

norm_supercell_df <- as.data.table(norm_supercell)
norm_supercell_df$CellId <- rownames(norm_supercell)
norm_supercell_df$sample_id <- supercell_clust$sample


fwrite(norm_supercell_df, paste0("dat_cytofruv_supercell_gam", gam, "_setting",setting,".csv"))

writeLines(
    text = unlist(tic.log(format = TRUE)),
    con = paste0("cytofruv_supercell_gam", gam, "_setting",setting,".txt")
)


