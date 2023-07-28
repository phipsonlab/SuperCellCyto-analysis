library(data.table)
library(SuperCellCyto)
library(BiocParallel)
library(here)

outdir <- here("output", "bodenmiller_cytof", format(Sys.Date(), "%Y%m%d"))
dir.create(outdir, recursive = TRUE)

cell_dat <- fread(here("data", "bodenmiller_cytof", "cell_dat.csv"))
cell_info <- fread(here("data", "bodenmiller_cytof", "cell_info.csv"))
markers <- fread(here("data", "bodenmiller_cytof", "markers_info.csv"))

# keep only the cell type and cell state markers
markers <- markers[markers$marker_class != "none",]
markers[marker_name == "HLA-DR", marker_name := "HLA.DR"]

markers_name <- markers$marker_name
# subset data to only the markers
cell_dat <- cell_dat[, markers_name, with = FALSE]

# asinh transformation with co-factor 5
asinh_colnames <- paste0(markers_name, "_asinh_cf5")
cell_dat[, (asinh_colnames) := lapply(.SD, function(x) asinh(x / 5)), .SDcols = markers_name]

# Attach the cell id and sample id
cell_dat[, cell_id := paste0("cell_", seq(1, nrow(cell_dat)))]
cell_dat[, sample_id := cell_info$sample_id]

cell_info$cell_id <- cell_dat$cell_id

fwrite(cell_info, here(outdir, "cell_info_with_cell_id.csv"))

# ---- Run SuperCellCyto ----
BPPARAM <- MulticoreParam(workers = parallel::detectCores() - 1, tasks = length(unique(cell_dat$sample_id)))

supercell_obj <- runSuperCellCyto(
    dt = cell_dat,
    markers = asinh_colnames,
    sample_colname = "sample_id",
    cell_id_colname = "cell_id",
    gam = 20,
    BPPARAM = BPPARAM,
    load_balancing = TRUE
)

supercell_mat <- supercell_obj$supercell_expression_matrix
supercell_cell_map <- supercell_obj$supercell_cell_map

fwrite(supercell_mat, here(outdir, "supercell_gamma20_exp_mat.csv"))
fwrite(supercell_cell_map, here(outdir, "supercell_gamma20_cell_map.csv"))
saveRDS(supercell_obj$supercell_object, here(outdir, "supercell_gamma20_obj.rds"))









