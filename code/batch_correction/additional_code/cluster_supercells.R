# Use clustering as a proxy to cell type annotation.
# We will use 20 as the number of metaclusters and 10x10 grid.
# We will also used seed 1234. The same as Marie's seed. Maybe not necessary..

library(data.table)
library(Spectre)
library(viridis)
library(RColorBrewer)

set.seed(42)

# supercells <- fread("output/trussart_cytofruv/20230515_supercell_out/supercellExpMat_rep1.csv")
# use this one as there is umap coordinates already..
supercells <- fread("output/trussart_cytofruv/20230524/supercellExpMat_rep1_umap.csv")

markers <- fread("data/trussart_cytofruv/metadata/Panel.csv")
markers_celltypes <- markers[marker_class == 'type']

# needed as this is what the columns in the supercell expression matrix are named like.
# Stupid me, but whatever
markers_celltypes[, channel_and_antigen := paste(fcs_colname, antigen, sep = "_")]

markers_to_cluster_on <- paste0(markers_celltypes$channel_and_antigen, "_asinh_cf5")

# --- cluster uncorrected supercells ----
supercells <- run.flowsom(
    dat = supercells,
    use.cols = markers_to_cluster_on,
    xdim = 10,
    ydim = 10,
    meta.k = 20,
    clust.seed = 42,
    meta.seed = 42
)
fwrite(supercells, "output/trussart_cytofruv/20240111/supercells_uncorrected_clustered.csv")
# supercells <- fread("output/trussart_cytofruv/20240111/supercells_uncorrected_clustered.csv")

plt <- ggplot(supercells, aes(x = UMAP_X, y = UMAP_Y, colour = factor(FlowSOM_metacluster))) +
    geom_point(size = 0.1) +
    facet_wrap(~ sample_id, ncol = 4) +
    scale_color_viridis(discrete = TRUE, option = "turbo") +
    guides(color=guide_legend(ncol=5, override.aes = list(size=5))) +
    theme_classic() +
    theme(legend.position = "bottom")

ggsave(
    filename = "~/Documents/supercell/manuscript/figures_v2_postReview/batch_correction/umap_clustered.png",
    plot = plt,
    width = 10,
    height = 15
)

plt <- ggplot(supercells, aes(x = UMAP_X, y = UMAP_Y, colour = factor(batch))) +
    geom_point(size = 0.1) +
    facet_wrap(~ sample_id, ncol = 4) +
    scale_color_brewer(palette = "Set1") +
    guides(color=guide_legend(ncol=5, override.aes = list(size=5))) +
    theme_classic() +
    theme(legend.position = "bottom")

ggsave(
    filename = "~/Documents/supercell/manuscript/figures_v2_postReview/batch_correction/umap_batch.png",
    plot = plt,
    width = 10,
    height = 15
)

supercells <- supercells[order(batch)]
plt <- ggplot(supercells, aes(x = UMAP_X, y = UMAP_Y, colour = factor(batch))) +
    geom_point(size = 0.01) +
    scale_color_manual(values=c("1" = "yellow", "2" = "dodgerblue")) +
    guides(color=guide_legend(ncol=5, override.aes = list(size=5))) +
    theme_classic() +
    theme(legend.position = "bottom") +
    labs(color = "batch", title = "Uncorrected")
ggsave(
    filename = "~/Documents/supercell/manuscript/figures_v2_postReview/batch_correction/umap_uncorrected_batch_v2.png",
    plot = plt,
    width = 10,
    height = 10
)


# --- cluster cytofruv output ----
cytofruv_supercells <- fread("output/trussart_cytofruv/20230524/supercellExpMat_postCytofRUV.csv")
cytofruv_supercells <- run.flowsom(
    dat = cytofruv_supercells,
    use.cols = markers_to_cluster_on,
    xdim = 10,
    ydim = 10,
    meta.k = 20,
    clust.seed = 42,
    meta.seed = 42
)
fwrite(cytofruv_supercells, "output/trussart_cytofruv/20240111/supercells_cytofruv_clustered.csv")
# cytofruv_supercells <- fread("output/trussart_cytofruv/20240111/supercells_cytofruv_clustered.csv")

plt <- ggplot(cytofruv_supercells, aes(x = UMAP_X, y = UMAP_Y, colour = factor(FlowSOM_metacluster))) +
    geom_point(size = 0.1) +
    facet_wrap(~ sample_id, ncol = 4) +
    scale_color_viridis(discrete = TRUE, option = "turbo") +
    guides(color=guide_legend(ncol=5, override.aes = list(size=5))) +
    theme_classic() +
    theme(legend.position = "bottom")

ggsave(
    filename = "~/Documents/supercell/manuscript/figures_v2_postReview/batch_correction/umap_cytofruv_clustered.png",
    plot = plt,
    width = 10,
    height = 15
)

plt <- ggplot(cytofruv_supercells, aes(x = UMAP_X, y = UMAP_Y, colour = factor(batch))) +
    geom_point(size = 0.1) +
    facet_wrap(~ sample_id, ncol = 4) +
    scale_color_brewer(palette = "Set1") +
    guides(color=guide_legend(ncol=5, override.aes = list(size=5))) +
    theme_classic() +
    theme(legend.position = "bottom")

ggsave(
    filename = "~/Documents/supercell/manuscript/figures_v2_postReview/batch_correction/umap_cytofruv_batch.png",
    plot = plt,
    width = 10,
    height = 15
)

cytofruv_supercells <- cytofruv_supercells[order(batch)]
plt <- ggplot(cytofruv_supercells, aes(x = UMAP_X, y = UMAP_Y, colour = factor(batch))) +
    geom_point(size = 0.01) +
    scale_color_manual(values=c("1" = "yellow", "2" = "dodgerblue")) +
    guides(color=guide_legend(ncol=5, override.aes = list(size=5))) +
    theme_classic() +
    theme(legend.position = "bottom") +
    labs(color = "batch", title = "CytofRUV")
ggsave(
    filename = "~/Documents/supercell/manuscript/figures_v2_postReview/batch_correction/umap_cytofruv_batch_v2.png",
    plot = plt,
    width = 10,
    height = 10
)

# --- cluster cyCombine output ----
cycombine_supercells <- fread("output/trussart_cytofruv/20230704/supercellExpMat_postCycombine.csv")
cycombine_supercells <- run.flowsom(
    dat = cycombine_supercells,
    use.cols = markers_to_cluster_on,
    xdim = 10,
    ydim = 10,
    meta.k = 20,
    clust.seed = 42,
    meta.seed = 42
)
fwrite(cycombine_supercells, "output/trussart_cytofruv/20240111/supercells_cycombine_clustered.csv")
# cycombine_supercells <- fread("output/trussart_cytofruv/20240111/supercells_cycombine_clustered.csv")

plt <- ggplot(cycombine_supercells, aes(x = UMAP_X, y = UMAP_Y, colour = factor(FlowSOM_metacluster))) +
    geom_point(size = 0.1) +
    facet_wrap(~ sample, ncol = 4) +
    scale_color_viridis(discrete = TRUE, option = "turbo") +
    guides(color=guide_legend(ncol=5, override.aes = list(size=5))) +
    theme_classic() +
    theme(legend.position = "bottom")

ggsave(
    filename = "~/Documents/supercell/manuscript/figures_v2_postReview/batch_correction/umap_cycombine_clustered.png",
    plot = plt,
    width = 10,
    height = 15
)

plt <- ggplot(cycombine_supercells, aes(x = UMAP_X, y = UMAP_Y, colour = factor(batch))) +
    geom_point(size = 0.1) +
    facet_wrap(~ sample, ncol = 4) +
    scale_color_brewer(palette = "Set1") +
    guides(color=guide_legend(ncol=5, override.aes = list(size=5))) +
    theme_classic() +
    theme(legend.position = "bottom")

ggsave(
    filename = "~/Documents/supercell/manuscript/figures_v2_postReview/batch_correction/umap_cycombine_batch.png",
    plot = plt,
    width = 10,
    height = 15
)

cycombine_supercells <- cycombine_supercells[order(batch)]
plt <- ggplot(cycombine_supercells, aes(x = UMAP_X, y = UMAP_Y, colour = factor(batch))) +
    geom_point(size = 0.01) +
    scale_color_manual(values=c("1" = "yellow", "2" = "dodgerblue")) +
    guides(color=guide_legend(ncol=5, override.aes = list(size=5))) +
    theme_classic() +
    theme(legend.position = "bottom") +
    labs(color = "batch", title = "cyCombine")
ggsave(
    filename = "~/Documents/supercell/manuscript/figures_v2_postReview/batch_correction/umap_cycombine_batch_v2.png",
    plot = plt,
    width = 10,
    height = 10
)












