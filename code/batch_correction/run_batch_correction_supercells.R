library(here)
library(CytofRUV)
library(CATALYST)
library(data.table)
library(BiocParallel)
library(Spectre)
library(rsvd) # For cytofruv
library(cyCombine)
library(bluster)
library(scales)
library(parallel)
library(batchelor)

source(here("code", "batch_correction", "cytofRUV_functions.R"))

outdir <- here("output", "trussart_cytofruv", format(Sys.Date(), "%Y%m%d"))
dir.create(outdir)

# Read in supercell exp matrix
supercell_exp_mat <- fread(here("output", "trussart_cytofruv", "20230515_supercell_out", "supercellExpMat_rep1.csv"))

markers <- names(supercell_exp_mat)[c(1:31)]

sample_metadata <- fread(here("data", "trussart_cytofruv", "metadata", "Metadata.csv"))

set.seed(42)

# ---- Run CytofRUV ----
# Create sce object
supercell_exp_mat_transposed <- as.matrix(t(supercell_exp_mat[, markers, with = FALSE]))
rownames(supercell_exp_mat_transposed) <- markers
colnames(supercell_exp_mat_transposed) <- supercell_exp_mat$SuperCellId

# Marry the supercell id and the metadata
supercell_metadata <- merge.data.table(
    x = supercell_exp_mat[, c("sample_id", "SuperCellId")],
    y = sample_metadata,
    by = "sample_id"
)

# Need to make sure the order is the same as the transposed supercell exp matrix.
supercell_metadata[, SuperCellId := factor(SuperCellId, levels = colnames(supercell_exp_mat_transposed))]
supercell_metadata <- DataFrame(supercell_metadata[order(SuperCellId)])
rownames(supercell_metadata) <- supercell_metadata$SuperCellId

# Create the marker (panel information)
panel_info <- fread(here("data", "trussart_cytofruv", "metadata", "Panel.csv"))
# create a new column which concatenates the reporter and the antigen
panel_info[, reporter_marker := factor(paste(fcs_colname, antigen, sep = "_"), levels = gsub("_asinh_cf5", "", markers))]
panel_info <- DataFrame(panel_info[order(reporter_marker)])
rownames(panel_info) <- paste0(panel_info$reporter_marker, "_asinh_cf5")

sce <- SingleCellExperiment(
    list(exprs = supercell_exp_mat_transposed),
    colData = supercell_metadata,
    rowData = panel_info
)

# cluster data
celltype_markers <- as.vector(rownames(panel_info[panel_info$marker_class == 'type', ]))
# Not sure if 20 is a good number.
sce <- cluster_data(sce, 42, markers_to_use = celltype_markers, 20)

saveRDS(sce, here(outdir, "sce_object.rds"))

# Use the CLL2 and the HCL1 samples as pseudo replicates.

# Extract clustered data as data.frame
supercell_clust <- data.frame(
    sample = sce$sample_id,
    cluster = cluster_ids(sce, "meta20"),
    t(assay(sce, "exprs"))
)

norm_supercell <- run_RUVIII(
    data = supercell_clust,
    norm_clusters = list(seq(1, 20), seq(1, 20)),
    k = 5,
    rep_samples = list(c("CLL2_B1", "CLL2_B2"), c("HC1_B1", "HC1_B2"))
)

norm_supercell_df <- as.data.table(norm_supercell)
norm_supercell_df$SuperCellId <- rownames(norm_supercell)
norm_supercell_df$sample_id <- supercell_clust$sample

# To get the batch information
norm_supercell_df <- merge.data.table(
    x = norm_supercell_df,
    y = sample_metadata,
    by = "sample_id"
)

fwrite(norm_supercell_df, here(outdir, "supercellExpMat_postCytofRUV.csv"))

# ---- Run cyCombine ----

supercell_for_cycombine <- merge.data.table(
    supercell_exp_mat, sample_metadata,
    by = "sample_id"
)

setnames(supercell_for_cycombine, "sample_id", "sample")

cycombine_corrected <- batch_correct(
    df = supercell_for_cycombine[, c(markers, "batch", "sample", "condition"), with = FALSE],
    xdim = 8,
    ydim = 8,
    seed = 42,
    markers = markers,
    covar = "condition"
)

cycombine_corrected <- cbind(
    cycombine_corrected,
    supercell_for_cycombine[, c("SuperCellId", "patient_id")]
)

fwrite(cycombine_corrected, here(outdir, "supercellExpMat_postCycombine.csv"))

# ---- Run UMAP ----

cytofruv_corrected_supercell <- fread(here(outdir, "supercellExpMat_postCytofRUV.csv"))
cycombine_corrected_supercell <- fread(here(outdir, "supercellExpMat_postCycombine.csv"))

# Run actual UMAP
supercell_exp_mat <- run.fast.umap(dat = supercell_exp_mat, use.cols = markers)
supercell_exp_mat <- merge.data.table(supercell_exp_mat, sample_metadata, by = "sample_id")
supercell_exp_mat[, batch := factor(batch)]

cytofruv_corrected_supercell <- run.fast.umap(dat = cytofruv_corrected_supercell, use.cols = markers)
cycombine_corrected_supercell <- run.fast.umap(dat = cycombine_corrected_supercell, use.cols = markers)

cytofruv_corrected_supercell[, batch := factor(batch)]

fwrite(supercell_exp_mat, here(outdir, "supercellExpMat_rep1_umap.csv"))
fwrite(cytofruv_corrected_supercell, here(outdir, "supercellExpMat_postCytofRUV.csv"))
fwrite(cycombine_corrected_supercell, here(outdir, "supercellExpMat_postCycombine.csv"))
