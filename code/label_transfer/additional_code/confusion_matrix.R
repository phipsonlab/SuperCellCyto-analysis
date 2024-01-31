library(pheatmap)
library(data.table)
library(viridis)
library(scales)
library(ggplot2)
library(RColorBrewer)
library(here)
library(grid)

# Prepare data

harmony_res <- fread(here("output", "label_transfer", "harmony_knn.csv"))
rpca_res <- fread(here("output", "label_transfer", "seurat_rPCA.csv"))

rpca_res_singlecell <- fread(here("output", "label_transfer", "seurat_rPCA_singlecell.csv"))
harmony_res_singlecell <- fread(here("output", "label_transfer", "harmony_knn_singlecell.csv"))

harmony_res <- harmony_res[Gated_Population != "unassigned"]
rpca_res <- rpca_res[Gated_Population != "unassigned"]

rpca_res_singlecell <- rpca_res_singlecell[Gated_Population != "unassigned"]
harmony_res_singlecell <- harmony_res_singlecell[Gated_Population != "unassigned"]


conf_mat_harmony_supercell <- data.table(with(harmony_res, table(predicted_population, Gated_Population)))

ggplot(conf_mat_harmony_supercell, aes(x = Gated_Population, y = predicted_population, fill= N)) +
    geom_tile() +
    geom_text(aes(label=N), color = "#81D8D0", size=3) +
    scale_fill_viridis(option="inferno") +
    labs(x = "Actual",y = "Prediction", fill="Count", title = "Harmony Combined with kNN Supercell") +
    theme_minimal() +
    theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank()
    )

ggsave("~/Documents/supercell/manuscript/figures_v2_postReview/supp_figure_harmony_knn_supercell.png",
       dpi=600,
       width=10,
       height=10,
       bg = "white")

conf_mat_harmony_singlecell <- data.table(with(harmony_res_singlecell, table(predicted_population, Gated_Population)))

ggplot(conf_mat_harmony_singlecell, aes(x = Gated_Population, y = predicted_population, fill= N)) +
    geom_tile() +
    geom_text(aes(label=N), color = "#81D8D0", size=3) +
    scale_fill_viridis(option="inferno") +
    labs(x = "Actual",y = "Prediction", fill="Count", title = "Harmony Combined with kNN Single Cell") +
    theme_minimal() +
    theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank()
    )

ggsave("~/Documents/supercell/manuscript/figures_v2_postReview/supp_figure_harmony_knn_singlecell.png",
       dpi=600,
       width=10,
       height=10,
       bg = "white")

conf_mat_rpca_supercell <- data.table(with(rpca_res, table(predicted_population, Gated_Population)))

ggplot(conf_mat_rpca_supercell, aes(x = Gated_Population, y = predicted_population, fill= N)) +
    geom_tile() +
    geom_text(aes(label=N), color = "#81D8D0", size=3) +
    scale_fill_viridis(option="inferno") +
    labs(x = "Actual",y = "Prediction", fill="Count", title = "Seurat rPCA Supercell") +
    theme_minimal() +
    theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank()
    )

ggsave("~/Documents/supercell/manuscript/figures_v2_postReview/figure_seurat_rpca_supercell.png",
       dpi=600,
       width=10,
       height=10,
       bg = "white")

conf_mat_rpca_singlecell <- data.table(with(rpca_res_singlecell, table(predicted_population, Gated_Population)))

ggplot(conf_mat_rpca_singlecell, aes(x = Gated_Population, y = predicted_population, fill= N)) +
    geom_tile() +
    geom_text(aes(label=N), color = "#81D8D0", size=3) +
    scale_fill_viridis(option="inferno") +
    labs(x = "Actual",y = "Prediction", fill="Count", title = "Seurat rPCA Single Cell") +
    theme_minimal() +
    theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank()
    )

ggsave("~/Documents/supercell/manuscript/figures_v2_postReview/supp_figure_seurat_rpca_singlecell.png",
       dpi=600,
       width=10,
       height=10,
       bg = "white")















