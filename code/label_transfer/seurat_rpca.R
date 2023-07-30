library(Seurat)
library(data.table)
library(here)
library(tictoc)

# we will use the supercell already generated from the clustering benchmarking output
# for levine_32dim data.

cytof_dt <- fread(here("data",
                       "clustering_benchmarking",
                       "v2_pure_R_pipeline", "levine_32dim", "levine_32_asinh.csv"))
panel_info <- fread(
  here("data", "clustering_benchmarking", "v2_pure_R_pipeline", "levine_32dim",
       "levine_32_asinh_markers_info.csv"))
cell_type_markers <- panel_info[marker_class == "type"]

citeseq_dat <- readRDS(here("data",
                            "haas_bm",
                            "97ADT_462RNA_BM_healthy_young_old.rds"))

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

supercell_cytof <- fread(here("output", "clustering_benchmarking", "20230511",
                              "levine_32dim", "supercell_runs", "supercellExpMat_gamma20.csv"))
supercell_cytof_trans <- supercell_cytof[, common_markers, with = FALSE]
setnames(supercell_cytof_trans, common_markers, names(common_markers))

supercell_cytof_trans <- transpose(supercell_cytof_trans)
colnames(supercell_cytof_trans) <- supercell_cytof$SuperCellId
rownames(supercell_cytof_trans) <- names(common_markers)

cytof_seurat <- CreateSeuratObject(supercell_cytof_trans)
cytof_seurat <- AddMetaData(
  object = cytof_seurat,
  metadata = rep("cytof", ncol(cytof_seurat)),
  col.name = "source"
)

tic.clearlog()

tic("Seurat rPCA")

# Run rPCA
anchors <- FindTransferAnchors(
  reference = citeseq_common_ab,
  query = cytof_seurat,
  features = names(common_markers),
  normalization.method = 'LogNormalize',
  reduction = 'rpca',
  dims = 1:length(common_markers),
  approx.pca = FALSE,
  k.anchor = 20,
  npcs = length(common_markers),
)

cytof_labels <- TransferData(
    anchorset = anchors,
    refdata = Idents(citeseq_common_ab),
    weight.reduction='rpca.ref'
)

toc(log = TRUE, quiet = TRUE)

writeLines(
  text = unlist(tic.log(format = TRUE)),
  con = "output/label_transfer/seurat_rPCA_runtime.txt"
)
tic.clearlog()


Idents(cytof_seurat) <- cytof_labels$predicted.id
cytof_seurat <- AddMetaData(
  object = cytof_seurat,
  metadata = cytof_labels$prediction.score.max,
  col.name = "prediction_confidence"
)


# Store the annotations into a data.frame
# Then join it to the mapping of supercell id and cell id.
# So we can resolve the predicted cell type annotation of individual cell.
# Then we join it with the true cell type annotation we isolated before.

cytof_supercell_annotations <- data.table(
    SuperCellID = colnames(cytof_seurat),
    predicted_population = Idents(cytof_seurat),
    prediction_confidence = cytof_seurat$prediction_confidence
)

supercell_cell_map <- fread(here("output", "explore_supercell_purity_clustering", "20230511",
                                 "levine_32dim", "supercell_runs", "supercellCellMap_gamma20.csv"))

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

fwrite(cytof_supercell_annotations, here("output", "label_transfer", "seurat_rPCA.csv"))
