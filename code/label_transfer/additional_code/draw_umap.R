library(data.table)
library(Spectre)

cytof_dt <- fread("data/explore_supercell_purity_clustering/levine_32dim/levine_32_asinh.csv")

cytof_dt <- run.umap(
    cytof_dt,
    names(cytof_dt)[40:71],
    fast = TRUE
)

harmony_res <- fread(here("output", "label_transfer", "harmony_knn.csv"))
rpca_res <- fread(here("output", "label_transfer", "seurat_rPCA.csv"))

rpca_res_singlecell <- fread(here("output", "label_transfer", "seurat_rPCA_singlecell.csv"))
harmony_res_singlecell <- fread(here("output", "label_transfer", "harmony_knn_singlecell.csv"))


harmony_res$CellId <- factor(harmony_res$CellId, levels=cytof_dt$CellId)
rpca_res$CellId <- factor(rpca_res$CellId, levels=cytof_dt$CellId)
rpca_res_singlecell$CellId <- factor(rpca_res_singlecell$CellId, levels=cytof_dt$CellId)
harmony_res_singlecell$CellId <- factor(harmony_res_singlecell$CellId, levels=cytof_dt$CellId)

harmony_res <- harmony_res[order(CellId)]
rpca_res <- rpca_res[order(CellId)]
rpca_res_singlecell <- rpca_res_singlecell[order(CellId)]
harmony_res_singlecell <- harmony_res_singlecell[order(CellId)]

cytof_dt$harmony_supercell <- harmony_res$predicted_population
cytof_dt$rpca_supercell <- rpca_res$predicted_population
cytof_dt$rpca_singlecell <- rpca_res_singlecell$predicted_population
cytof_dt$harmony_singlecell <- harmony_res_singlecell$predicted_population

cytof_dt_no_unassigned <- cytof_dt[Gated_Population != 'unassigned']

cytof_dt_no_unassigned_molten <- melt(cytof_dt_no_unassigned,
     id.vars = c("CellId", "UMAP_X", "UMAP_Y"),
     measure.vars = c("harmony_supercell", "rpca_supercell", "harmony_singlecell", "rpca_singlecell"))

ggplot(cytof_dt_no_unassigned_molten, aes(x = UMAP_X, y = UMAP_Y, colour=value)) +
    geom_point(size=0.1) +
    facet_wrap(~ variable, nrow = 2, ncol = 2) +
    theme_classic() +
    theme(legend.position = "bottom") +
    guides(colour = guide_legend(override.aes = list(size=6))) +
    scale_color_viridis(option = "turbo", discrete = TRUE) +
    labs(colour = "Predicted cell type")

ggsave("~/Documents/supercell/manuscript/figures_v2_postReview/supp_figure_label_transfer_umap.png",
       dpi=600,
       width=15,
       height=15,
       bg = "white")

fwrite(cytof_dt_no_unassigned, "output/label_transfer/umap.csv")




