library(data.table)
library(SuperCellCyto)
library(parallel)
library(BiocParallel)
library(HDCytoData)
library(tictoc)

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

BPPARAM <- MulticoreParam(workers = detectCores() - 1, tasks = length(unique(cell_dat$sample_id)))

tic.clearlog()

for (rep in seq(2)) {
  tic(paste("Supercell run", rep))

  supercell_obj <- runSuperCellCyto(
    dt = cell_dat,
    markers = markers_name_asinh,
    sample_colname = "sample_id",
    cell_id_colname = "cell_id",
    gam = 20,
    BPPARAM = BPPARAM,
    load_balancing = TRUE
  )

  toc(log = TRUE, quiet = TRUE)
}

writeLines(
  text = unlist(tic.log(format = TRUE)),
  con = "output/krieg_melanoma/supercell_benchmark.txt"
)





