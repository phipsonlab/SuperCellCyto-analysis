library(Seurat)
library(data.table)
library(harmony)
library(class)
library(stringr)
library(tictoc)

# we will use the supercell already generated from the clustering benchmarking output
# for levine_32dim data.

cytof_dt <- fread("data/clustering_benchmarking/v2_pure_R_pipeline/levine_32dim/levine_32_asinh.csv")
panel_info <- fread("data/clustering_benchmarking/v2_pure_R_pipeline/levine_32dim/levine_32_asinh_markers_info.csv")
cell_type_markers <- panel_info[marker_class == "type"]

citeseq_dat <- readRDS("data/haas_bm/97ADT_462RNA_BM_healthy_young_old.rds")

ab_names <- rownames(citeseq_dat@assays$AB@counts)
ab_names[which(str_detect(ab_names, "HLA"))]

# For each marker, see potential candidate
potential_cand <- lapply(cell_type_markers$marker_name, function(m) {
  ab_names[which(str_detect(ab_names, paste0("\\b", m, "\\b")))]
})
names(potential_cand) <- cell_type_markers$marker_name

common_markers <- unlist(potential_cand)

# create new seurat object with just the common markers
citeseq_common_ab <- CreateSeuratObject(
  GetAssayData(citeseq_dat, assay = "AB", slot="count")[common_markers, ]
)

# Extract metadata and cell type
citeseq_common_ab <- AddMetaData(
  object = citeseq_common_ab,
  metadata = rep("citeseq", ncol(citeseq_common_ab)),
  col.name = "source"
)
cell_types <- citeseq_dat$Cell_Type4
Idents(citeseq_common_ab) <- cell_types

# Normalise data using CLR
citeseq_common_ab <- NormalizeData(
  citeseq_common_ab,
  normalization.method = "CLR",
  margin = 1
)

# Reverse the common markers so anything ending with -AB is the
# name and anything not ending with -AB are the values
common_markers <- setNames(paste0(names(common_markers), "_asinh_cf5"), common_markers)

supercell_cytof <- fread("output/clustering_benchmarking/20230511/levine_32dim/supercell_runs/supercellExpMat_gamma20.csv")

# Keep only the common markers and the supercell id
supercell_cytof <- supercell_cytof[, c("SuperCellId", common_markers), with = FALSE]
setnames(supercell_cytof, "SuperCellId", "id")
setnames(supercell_cytof, common_markers, names(common_markers))


citeseq_common_ab_dt <- as.data.frame(t(GetAssayData(citeseq_common_ab)))
citeseq_common_ab_dt$id <- rownames(citeseq_common_ab_dt)


combined_dt <- rbind(citeseq_common_ab_dt, supercell_cytof)
combined_dt$id <- NULL


# Setup metadata
metadata <- data.frame(
  source = c(rep("citeseq", nrow(citeseq_common_ab_dt)), rep("cytof", nrow(supercell_cytof))),
  cell_supercell_id = c(citeseq_common_ab_dt$id, supercell_cytof$id)
)

tic.clearlog()

tic("Harmony kNN")

# Run Harmony
harmonyObj <- HarmonyMatrix(
    data_mat = combined_dt,
    meta_data = metadata,
    vars_use = 'source',
    return_object = TRUE,
    do_pca = FALSE,
)

integrated_data <- data.table(t(harmonyObj$Z_corr))
integrated_data <- cbind(integrated_data, metadata)

cytof_to_predict <- integrated_data[source == "cytof"]
cytof_to_predict[, source := NULL]
cytof_to_predict[, cell_supercell_id := NULL]

citeseq_ref <- integrated_data[source == "citeseq"]
citeseq_ref[, source := NULL]
citeseq_ref[, cell_supercell_id := NULL]

# Run KNN
knn_prediction <- knn(
  train = citeseq_ref,
  test = cytof_to_predict,
  cl = cell_types,
  k = 1
)

toc(log = TRUE, quiet = TRUE)

writeLines(
  text = unlist(tic.log(format = TRUE)),
  con = "output/label_transfer/harmony_knn_runtime.txt"
)
tic.clearlog()

# Store the annotations into a data.frame
# Then join it to the mapping of supercell id and cell id.
# So we can resolve the predicted cell type annotation of individual cell.
# Then we join it with the true cell type annotation we isolated before.
cytof_supercell_annotations <- data.table(
    SuperCellID = metadata[metadata$source == "cytof"]$cell_supercell_id,
    predicted_population = knn_prediction
)

supercell_cell_map <- fread("output/clustering_benchmarking/20230511/levine_32dim/supercell_runs/supercellCellMap_gamma20.csv")

cytof_supercell_annotations <- merge.data.table(
    x = cytof_supercell_annotations,
    y = supercell_cell_map,
    by = "SuperCellID"
)

cytof_supercell_annotations <- merge.data.table(
    x = cytof_supercell_annotations,
    y = cytof_dt[, c("CellId", "Gated_Population"), with = FALSE],
    by = "CellId"
)

fwrite(cytof_supercell_annotations, "output/label_transfer/harmony_knn.csv")


