library(data.table)
library(ggplot2)
library(limma)
library(speckle)
library(SuperCellCyto)
library(parallel)
library(here)
library(BiocParallel)
library(Spectre)
library(pheatmap)
library(scales)
library(cyCombine)
library(HDCytoData)
library(ggrepel)
library(ggridges)
library(tictoc)

for (rep in seq(2)) {
  sce <- Krieg_Anti_PD_1_SE()

  cell_info <- data.table(as.data.frame(rowData(sce)))
  markers <- data.table(as.data.frame(colData(sce)))
  cell_dat <- data.table(as.data.frame(assay(sce)))

  cell_dat <- cbind(cell_dat, cell_info)

  # keep only the cell type and cell state markers
  markers <- markers[marker_class != "none"]
  markers_name <- markers$marker_name

  # asinh transformation with co-factor 5
  markers_name_asinh <- paste0(markers_name, "_asinh_cf5")

  monocyte_markers <- paste0(
    c(
      "CD14", "CD33", "HLA-DR", "ICAM-1", "CD64", "CD141", "CD86", "CD11c",
      "CD38", "CD274_PDL1", "CD11b"
    ),
    "_asinh_cf5"
  )

  cell_dat <- cell_dat[, c(markers_name, "group_id", "batch_id", "sample_id"), with = FALSE]

  # arc-sinh transformation with co-factor 5
  cell_dat[, (markers_name_asinh) := lapply(.SD, function(x) asinh(x / 5)), .SDcols = markers_name]

  # save ram, remove untransformed markers
  cell_dat[, c(markers_name) := NULL]

  # Change group field into factor
  cell_dat[, group_id := factor(group_id, levels = c("NR", "R"))]
  cell_dat[, cell_id := paste0("cell_", seq(nrow(cell_dat)))]

  tic("start supercell cycombine flowsom propeller")
  BPPARAM <- MulticoreParam(workers = detectCores() - 1, tasks = length(unique(cell_dat$sample_id)))

  supercell_obj <- runSuperCellCyto(
    dt = cell_dat,
    markers = markers_name_asinh,
    sample_colname = "sample_id",
    cell_id_colname = "cell_id",
    gam = 20,
    # BPPARAM = BPPARAM,
    # load_balancing = TRUE
  )

  supercell_mat <- supercell_obj$supercell_expression_matrix
  supercell_cell_map <- supercell_obj$supercell_cell_map

  sample_info <- unique(cell_info)
  supercell_mat <- merge.data.table(supercell_mat, sample_info)

  setnames(supercell_mat, "sample_id", "sample")
  setnames(supercell_mat, "batch_id", "batch")

  cycombine_corrected_sup <- batch_correct(
    df = supercell_mat[, c(markers_name_asinh, "batch", "sample", "group_id"), with = FALSE],
    xdim = 4,
    ydim = 4,
    seed = 42,
    markers = markers_name_asinh,
    covar = "group_id"
  )

  cycombine_corrected_sup <- run.flowsom(
    cycombine_corrected_sup,
    use.cols = markers_name_asinh,
    xdim = 20,
    ydim = 20,
    meta.k = 50,
    clust.seed = 42,
    meta.seed = 42
  )

  cycombine_corrected_sup$SuperCellId <- supercell_mat$SuperCellId

  expanded_supercell <- merge.data.table(
    supercell_cell_map,
    cycombine_corrected_sup[, c("SuperCellId", "FlowSOM_cluster", "FlowSOM_metacluster", "group_id", "batch")],
    by.x = "SuperCellID",
    by.y = "SuperCellId"
  )

  nsamples_min <- nrow(sample_info)
  clust_cnt <- table(expanded_supercell$FlowSOM_metacluster, expanded_supercell$Sample)
  clust_to_keep <- as.numeric(names(which(rowSums(clust_cnt > 3) >= nsamples_min)))

  expanded_supercell_sub <- expanded_supercell[FlowSOM_metacluster %in% clust_to_keep]
  supercell_mat_sub <- cycombine_corrected_sup[FlowSOM_metacluster %in% clust_to_keep]

  prop <- getTransformedProps(
    clusters = expanded_supercell_sub$FlowSOM_metacluster,
    sample = expanded_supercell_sub$Sample
  )

  sample_info[, sample_id := factor(sample_id, levels = colnames(prop$Counts))]
  sample_info <- sample_info[order(sample_id)]

  designAS <- model.matrix(~ 0 + sample_info$group_id + sample_info$batch_id)
  colnames(designAS) <- c("NR", "R", "batch29vs23")
  mycontr <- makeContrasts(NR - R, levels = designAS)

  test_res <- propeller.ttest(
    prop.list = prop, design = designAS, contrasts = mycontr,
    robust = TRUE, trend = FALSE, sort = TRUE
  )

  toc(log = TRUE, quiet = TRUE)


  setnames(cell_dat, "sample_id", "sample")
  setnames(cell_dat, "batch_id", "batch")

  tic("start cycombine flowsom propeller single cell")

  cycombine_corrected <- batch_correct(
    df = cell_dat[, c(markers_name_asinh, "batch", "sample", "group_id"), with = FALSE],
    xdim = 4,
    ydim = 4,
    seed = 42,
    markers = markers_name_asinh,
    covar = "group_id"
  )

  cycombine_corrected <- run.flowsom(
    cycombine_corrected,
    use.cols = markers_name_asinh,
    xdim = 20,
    ydim = 20,
    meta.k = 50,
    clust.seed = 42,
    meta.seed = 42
  )

  cycombine_corrected$cell_id <- cell_dat$cell_id

  nsamples_min <- nrow(sample_info)
  clust_cnt <- table(cycombine_corrected$FlowSOM_metacluster, cycombine_corrected$sample)
  clust_to_keep <- as.numeric(names(which(rowSums(clust_cnt > 3) >= nsamples_min)))

  cycombine_corrected_sub <- cycombine_corrected[FlowSOM_metacluster %in% clust_to_keep]

  prop <- getTransformedProps(
    clusters = cycombine_corrected_sub$FlowSOM_metacluster,
    sample = cycombine_corrected_sub$sample
  )

  sample_info[, sample_id := factor(sample_id, levels = colnames(prop$Counts))]
  sample_info <- sample_info[order(sample_id)]

  designAS <- model.matrix(~ 0 + sample_info$group_id + sample_info$batch_id)
  colnames(designAS) <- c("NR", "R", "batch29vs23")
  mycontr <- makeContrasts(NR - R, levels = designAS)

  test_res <- propeller.ttest(
    prop.list = prop, design = designAS, contrasts = mycontr,
    robust = TRUE, trend = FALSE, sort = TRUE
  )

  toc(log = TRUE, quiet = TRUE)
}

writeLines(
  text = unlist(tic.log(format = TRUE)),
  con = "output/da_test/supercell_vs_singlecell_duration_run1.txt"
)
















