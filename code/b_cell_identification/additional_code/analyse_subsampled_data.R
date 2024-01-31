library(data.table)
library(Spectre)
library(BiocParallel)
library(bluster)
library(ggplot2)
library(scales)
library(viridis)
library(pheatmap)

# subsample the cells to 415711 cells, same number of supercells we get.

dat <- fread("output/oetjen_b_cell_panel/20230511/cell_dat_asinh.csv")

proportion <- 415711 / nrow(dat)
subsampling_targets <- data.table(round(table(dat$sample) * proportion))
setnames(subsampling_targets, "V1", "sample")

set.seed(42)
# stratified subsamping based on the samples
dat_sub <- do.subsample(
    dat = dat,
    targets = subsampling_targets,
    divide.by = "sample",
    seed = 42
)

fwrite(dat_sub, "output/oetjen_b_cell_panel/20240112/dat_subsampled.csv")

dat_sub <- fread("output/oetjen_b_cell_panel/20240112/dat_subsampled.csv")

cell_type_markers <- paste0(
    c("CD19", "CD45", "CD10", "CD20", "CD27", "CD21", "CD38", "CD138"),
    "_asinh_cf",
    150
)

exp_mt_toCluster <- dat_sub[, cell_type_markers, with = FALSE]

clusters <- clusterRows(
    x = exp_mt_toCluster,
    BLUSPARAM = NNGraphParam(
        k = 3,
        cluster.fun = "louvain",
        BPPARAM = MulticoreParam(RNGseed = 42)
    )
)

dat_sub[, louvain_k3 := paste0("cluster_", clusters)]

fwrite(dat_sub, "output/oetjen_b_cell_panel/20240112/dat_subsampled_clustered.csv")

median_exp <- dat_sub[, lapply(.SD, median), by = louvain_k3, .SDcols = cell_type_markers]
median_exp_mat <- as.matrix(median_exp[, cell_type_markers, with = FALSE])
row.names(median_exp_mat) <- median_exp$louvain_k3
colnames(median_exp_mat) <- gsub("_asinh_cf150", "", colnames(median_exp_mat))

pheatmap(
    median_exp_mat,
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    display_numbers = TRUE,
    number_format = "%.2f",
    scale = "column",
    main = "Median Expression of Cell Type Markers Across Clusters"
)

annotation <- rep("CD45- and/or CD19-", 35)
# 2 = Transitional B cells
# 19 = CD20- CD10-
# 20 = Naïve Mature B cells
# 21 = Exhausted/Tissue-like Memory
# 22 = Resting Memory
# 31-33 = CD20- CD10-
# 34 = Pre-B cell
# 35 = Immature B cells
annotation[34] <- "Pre B cells"
annotation[35] <- "Immature B cells"
annotation[2] <- "Transitional B cells"
annotation[c(31:33)] <- "CD20- CD10-"
annotation[20] <- "Naive Mature B cells"
annotation[21] <- "Exhausted Tissue like Memory"
annotation[22] <- "Resting Memory"
# Not found: "Plasmablast", "Plasma cells", "Activated Mature"
annotation_df <- data.table(cluster = paste0("cluster_", seq(1, 35)), cell_type = annotation)

# The order is for ordering the heatmap so the same cell types are together
heatmap_annot_col <- annotation_df[order(cell_type)]
heatmap_annot_col <- heatmap_annot_col[cell_type != "CD45- and/or CD19-"]
heatmap_annot_col <- data.frame(heatmap_annot_col, row.names = "cluster")
annot_cell_types <- unique(heatmap_annot_col$cell_type)

annot_colours <- viridis(n = length(annot_cell_types), option = "turbo")
names(annot_colours) <- annot_cell_types

# Re-order heatmap to match annotation so the same cell types are together
median_exp_mat_noUnassigned <- median_exp_mat[rownames(heatmap_annot_col), ]

pheatmap(
    t(median_exp_mat_noUnassigned),
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    annotation_col = heatmap_annot_col,
    # scale = "column",
    main = "Median Expression of Cell Type Markers Across Clusters",
    annotation_colors = list(cell_type = annot_colours),
    angle_col = 90
)

# run umap
dat_sub <- run.umap(
    dat = dat_sub,
    use.cols = cell_type_markers,
    fast = TRUE
)

dat_sub[, louvain_k3 := paste0("cluster_", louvain_k3)]

# add the annotation in
dat_sub <- merge.data.table(
    dat_sub,
    annotation_df,
    by.x = "louvain_k3",
    by.y = "cluster"
)

fwrite(dat_sub, "output/oetjen_b_cell_panel/20240112/dat_subsampled_annotated.csv")
# dat_sub <- fread("output/oetjen_b_cell_panel/20240112/dat_subsampled_annotated.csv")

make.multi.plot(
    dat = dat_sub,
    x.axis = "UMAP_X",
    y.axis = "UMAP_Y",
    plot.by = cell_type_markers
)

# so the umap is pretty
cell_types <- sort(unique(dat_sub$cell_type))
cell_types <- cell_types[cell_types != "CD45- and/or CD19-"]
colour_scheme <- viridis(n = length(cell_types), option = "turbo")
names(colour_scheme) <- cell_types
colour_scheme["CD45- and/or CD19-"] <- "#D3D3D3"

# So unknown is at the bottom.
dat_sub[, cell_type := factor(cell_type, levels = names(colour_scheme))]
dat_sub <- dat_sub[order(-cell_type)]

ggplot(dat_sub, aes(x=UMAP_X, y=UMAP_Y)) +
    geom_point(aes(color=cell_type), size=0.1) +
    theme_classic() +
    theme(legend.title = element_blank()) +
    scale_colour_manual(values = colour_scheme) +
    ggtitle("Annotation") +
    guides(colour = guide_legend(override.aes = list(size=5)))

ggsave("output/oetjen_b_cell_panel/20240112/umap_annotation.png")
