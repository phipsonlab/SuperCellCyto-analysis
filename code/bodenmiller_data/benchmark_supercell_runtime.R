library(data.table)
library(SuperCellCyto)
library(BiocParallel)
library(here)
library(tictoc)

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

# ---- Run SuperCellCyto ----
BPPARAM <- MulticoreParam(workers = parallel::detectCores() - 1, tasks = length(unique(cell_dat$sample_id)))

tic.clearlog()

for (rep in seq(2)) {
  tic(paste("Supercell run", rep))

  supercell_obj <- runSuperCellCyto(
    dt = cell_dat,
    markers = asinh_colnames,
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
  con = "output/bodenmiller_cytof/supercell_benchmark.txt"
)











