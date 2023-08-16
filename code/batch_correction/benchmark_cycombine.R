library(here)
library(data.table)
library(cyCombine)
library(tictoc)

set.seed(42)

# Prepare single cell
dat <- fread("output/trussart_cytofruv/20230515_supercell_out/cell_dat_asinh.csv")
cols_to_keep <- names(dat)[c(32:68)]
dat <- dat[, cols_to_keep, with = FALSE]
markers <- names(dat)[7:37]
setnames(dat, "sample_id", "sample")

sample_metadata <- unique(dat[, c("sample", "condition", "patient_id", "batch")])

gams <- c(10,30,40)
all_dat <- lapply(gams, function(gam) {
  x = fread(paste0("output/trussart_cytofruv/20230805/supercell_mat_gamma", gam,".csv"))
  x = merge.data.table(x, sample_metadata, by.x="sample_id", by.y="sample")
  setnames(x, "sample_id", "sample")
  return(x)
})
names(all_dat) <- paste0("gam", gams)

# prepare gamma 20
x = fread("output/trussart_cytofruv/20230515_supercell_out/supercellExpMat_rep1.csv")
x = merge.data.table(x, sample_metadata, by.x="sample_id", by.y="sample")
setnames(x, "sample_id", "sample")
all_dat[["gam20"]] <- x

# add the single cell
all_dat[['singlecell']] <- dat

grid_size <- seq(6, 18, 2)
data_type <- names(all_dat)

# ---- Run cyCombine ----
tic.clearlog()
for (gr in grid_size) {

  for (data_name in data_type) {
    print(paste("Run cycombine grid", gr, "for", data_name))
    tic(paste("Run cycombine grid", gr, "for", data_name))

    cycombine_corrected <- batch_correct(
      df = all_dat[[data_name]][, c(markers, "batch", "sample", "condition"), with = FALSE],
      xdim = gr,
      ydim = gr,
      seed = 42,
      markers = markers,
      covar = "condition"
    )

    toc(log = TRUE, quiet = TRUE)
  }


}

writeLines(
  text = unlist(tic.log(format = TRUE)),
  con = "output/trussart_cytofruv/cycombine_benchmark_run2.txt"
)
