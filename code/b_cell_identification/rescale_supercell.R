library(SuperCellCyto)
library(data.table)

sc_obj <- readRDS("~/Documents/GitHub/SuperCellCyto-analysis/output/oetjen_b_cell_panel/20230511/supercellObj.rds")

dat <- fread("~/Documents/GitHub/SuperCellCyto-analysis/output/oetjen_b_cell_panel/20230511/cell_dat_asinh.csv")

gams <- c(10, 30, 40)

markers <- names(dat)[14:25]

for (gam in gams) {
  sc_recalc <- recomputeSupercells(
    dt = dat,
    sc_objects = sc_obj,
    markers = markers,
    sample_colname = "sample",
    cell_id_colname = "CellId",
    gam = gam
  )

  fwrite(sc_recalc$supercell_expression_matrix,
         paste0("~/Documents/GitHub/SuperCellCyto-analysis/output/oetjen_b_cell_panel/20230805/supercell_mat_gamma", gam, ".csv"))

  fwrite(sc_recalc$supercell_cell_map,
         paste0("~/Documents/GitHub/SuperCellCyto-analysis/output/oetjen_b_cell_panel/20230805/supercell_map_gamma", gam, ".csv"))
}

