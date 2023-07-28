# ---- Download Bodenmiller data ----
library(HDCytoData)
library(data.table)

sce <- Bodenmiller_BCR_XL_SE()

cell_info <- data.frame(rowData(sce))
markers <- data.frame(colData(sce))

cell_dat <- data.frame(assay(sce))

fwrite(cell_dat, "data/bodenmiller_cytof/cell_dat.csv")
fwrite(cell_info, "data/bodenmiller_cytof/cell_info.csv")
fwrite(markers, "data/bodenmiller_cytof/markers_info.csv")


