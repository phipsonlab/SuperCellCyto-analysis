# check CD34 expression of the progenitors
library(Seurat)

citeseq_dat <- readRDS("data/haas_bm/97ADT_462RNA_BM_healthy_young_old.rds")

table(citeseq_dat$Cell_Type4)

as.matrix(unique(citeseq_dat$Cell_Type4))

Idents(citeseq_dat) <- citeseq_dat$Cell_Type4

# Start with the progenitors
progenitor_subsets <- c('Late erythroid progenitor',
                 'Early erythroid progenitor',
                 'Lymphoid-primed multipotent progenitors',
                 'Plasmacytoid dendritic cell progenitors',
                 'Megakaryocyte progenitors',
                 'Erythro-myeloid progenitors',
                 'Eosinophil-basophil-mast cell progenitors',
                 'NK cell progenitors',
                 'HSCs & MPPs')
progenitor_dat <- subset(x = citeseq_dat, idents = progenitor_subsets)

FeaturePlot(progenitor_dat, features = c("CD34"), label = TRUE)
ggsave("~/Documents/supercell/manuscript/figures_v2_postReview/cd34_progenitor_umap_plot.png", width = 10, height = 10)

RidgePlot(progenitor_dat, features = c("CD34"), ncol = 1)
ggsave("~/Documents/supercell/manuscript/figures_v2_postReview/cd34_progenitor_ridge_plot.png", width = 10, height = 10)

# Late erythroid progenitor are mostly CD34-
# Plasmacytoid dendritic cell progenitors are also mostly CD34-
# 'Early erythroid progenitor' and 'Eosinophil-basophil-mast cell progenitors' are
# mix bag CD34- and CD34+, but mostly CD34-. So remove.

# Remove them and draw the ridgeplot again
progenitor_subsets <- c('Lymphoid-primed multipotent progenitors',
                        'Megakaryocyte progenitors',
                        'Erythro-myeloid progenitors',
                        'NK cell progenitors',
                        'HSCs & MPPs')
progenitor_dat <- subset(x = citeseq_dat, idents = progenitor_subsets)

RidgePlot(progenitor_dat, features = c("CD34"), ncol = 1)

# ggsave("~/Documents/supercell/manuscript/figures_v2_postReview/cd34_progenitor_ridge_plot_final.png", width = 10, height = 10)

FeaturePlot(progenitor_dat, features = c("CD34"), label = TRUE)
# ggsave("~/Documents/supercell/manuscript/figures_v2_postReview/cd34_progenitor_umap_plot_final.png", width = 10, height = 10)

# B cells
b_cells <- c('CD11c+ memory B cells',
             'Mature naive B cells',
             'Class switched memory B cells',
             'Nonswitched memory B cells',
             'Immature B cells',
             'Pre-pro-B cells',
             'Small pre-B cell',
             'pro-B cells')

bcell_dat <- subset(x = citeseq_dat, idents = b_cells)

RidgePlot(bcell_dat, features = c("CD19", "CD19-AB", "IgG-AB", "IgD-AB"), ncol = 1)
FeaturePlot(bcell_dat, features = c("CD19", "CD19-AB", "IgG-AB", "IgD-AB"), label = TRUE)
DimPlot(bcell_dat)

# Immature B cells kind of float between mature and pre or pro B cells, best leave it out.
# Also Pre-pro Bcells are best left out. It is not clear what it is doing.

# Reconcile (or map) the labels in the CITEseq data with the labels in the cytometry data
# For drawing giant UMAP with reconciled cell types in the CITEseq data.
reconcile_label <- function(label) {
    if (label == 'CD56brightCD16- NK cells'){
        return('CD16-_NK_cells')
    }
    else if (label == 'CD56dimCD16+ NK cells'){
        return('CD16+_NK_cells')
    }
    else if (label %in% c('CD4+ cytotoxic T cells', 'CD4+ memory T cells',
                          'Naive CD4+T cells')){
        return('CD4_T_cells')
    }
    else if (label %in% c("CD8+ central memory T cells",
                          "CD8+ effector memory T cells",
                          "CD8+ naive T cells",
                          "CD8+CD103+ tissue resident memory T cells")){
        return('CD8_T_cells')
    }
    else if (label %in% c('CD11c+ memory B cells',
                          'Mature naive B cells',
                          'Class switched memory B cells',
                          'Nonswitched memory B cells')){
        return('Mature_B_cells')
    }
    else if (label %in% c('Classical Monocytes', 'Non-classical monocytes')){
        return('Monocytes')
    }
    else if (label == 'Plasmacytoid dendritic cells'){
        return('pDCs')
    }
    else if (label == 'Plasma cells'){
        return('Plasma_B_cells')
    }
    else if (label == 'Small pre-B cell'){
        return('Pre_B_cells')
    }
    else if (label == 'pro-B cells'){
        return('Pro_B_cells')
    }
    else if (label %in% c('Lymphoid-primed multipotent progenitors', 'Megakaryocyte progenitors',
                       'Erythro-myeloid progenitors', 'NK cell progenitors', 'HSCs & MPPs')){
        return('CD34+_HSCs_and_HSPCs')
    }
    else
        return(NA)
}

reconciled_cell_types <- sapply(citeseq_dat$Cell_Type4, reconcile_label)
citeseq_dat$reconciled_cell_type <- reconciled_cell_types

# UMAP for original annotation
Idents(citeseq_dat) <- citeseq_dat$Cell_Type4
DimPlot(citeseq_dat) + theme(legend.position = "right") +
    guides(color=guide_legend(ncol=1, override.aes = list(size=5))) +
    ggtitle("CITEseq data original annotations")
ggsave("~/Documents/supercell/manuscript/figures_v2_postReview/label_transfer/supp_figure_umap_citeseq_origAnnot.png",
       width = 12, height = 10)

DimPlot(citeseq_dat, label=TRUE) + NoLegend()
ggsave("~/Documents/supercell/manuscript/figures_v2_postReview/supp_figure_umap_citeseq_origAnnot_noLegend.png",
       width = 10, height = 10)


# UMAP for reconciled annotation
Idents(citeseq_dat) <- citeseq_dat$reconciled_cell_type
DimPlot(citeseq_dat) + theme(legend.position = "right") +
    guides(color=guide_legend(ncol=1, override.aes = list(size=5))) +
    ggtitle("CITEseq data mapped annotations")
ggsave("~/Documents/supercell/manuscript/figures_v2_postReview/label_transfer/supp_figure_umap_citeseq_reconciledAnnot.png",
       width = 12, height = 10)

DimPlot(citeseq_dat, label=TRUE) + NoLegend()
ggsave("~/Documents/supercell/manuscript/figures_v2_postReview/supp_figure_umap_citeseq_reconciledAnnot_noLegend.png",
       width = 10, height = 10)

FeaturePlot(citeseq_dat, features = c("CD34", "CD34-AB"))
ggsave("~/Documents/supercell/manuscript/figures_v2_postReview/label_transfer/supp_figure_umap_citeseq_cd34.png",
       width = 10, height = 5)

FeaturePlot(citeseq_dat, features = c("CD19", "CD19-AB", "IgD-AB", "IgG-AB"), ncol=2)
ggsave("~/Documents/supercell/manuscript/figures_v2_postReview/label_transfer/supp_figure_umap_citeseq_cd19.png",
       width = 10, height = 10)







