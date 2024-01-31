library(data.table)
library(ggplot2)
library(scales)
library(RColorBrewer)


weighted_accuracies <- fread("output/label_transfer/weighted_accuracy.csv")
weighted_accuracies[, method := gsub("_reconciled", "", method)]
weighted_accuracies[, method := gsub("harmony", "Harmony_kNN", method)]
weighted_accuracies[, method := gsub("rpca", "rPCA", method)]

weighted_accuracies[, algorithm := ifelse(grepl("rpca", method), "rPCA", "Harmony_kNN")]
weighted_accuracies[, resolution := ifelse(grepl("supercell", method), "Supercell", "Single_cell")]

ggplot(weighted_accuracies, aes(x = method, fill=method, y=weighted_accuracy)) +
    geom_bar(stat="identity", position=position_dodge()) +
    scale_y_continuous(breaks=pretty_breaks(10), limits = c(0,1)) +
    scale_fill_brewer(palette = "Paired") +
    theme_classic() +
    theme(legend.position = "bottom") +
    labs(y = "Weighted accuracy", x = "Method", fill = "Method",
         title = "Weighted accuracy of label transfer") +
    guides(fill=guide_legend(ncol=2))

ggsave("~/Documents/supercell/manuscript/figures_v2_postReview/figure_weighted_accuracy.png",
       width = 7,
       height = 7)

accuracy_per_celltype <- fread("output/label_transfer/accuracy_per_celltype.csv")
setnames(accuracy_per_celltype, "V1", "cell_type")

accuracy_per_celltype <- melt(accuracy_per_celltype, id.vars = "cell_type")
accuracy_per_celltype[, method := gsub("_reconciled", "", variable)]
accuracy_per_celltype[, method := gsub("harmony", "Harmony_kNN", method)]
accuracy_per_celltype[, method := gsub("rpca", "rPCA", method)]
accuracy_per_celltype

ggplot(accuracy_per_celltype, aes(x = cell_type, fill=method, y=value)) +
    geom_bar(stat="identity", position=position_dodge()) +
    scale_y_continuous(breaks=pretty_breaks(10), limits = c(0,1)) +
    coord_flip() +
    theme_classic() +
    labs(y = "Accuracy", x = "Cell type", fill = "Method",
         title = "Accuracy of label transfer per cell type") +
    theme(
        # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "bottom"
        ) +
    scale_fill_brewer(palette = "Paired") +
    guides(fill=guide_legend(ncol=2))
ggsave("~/Documents/supercell/manuscript/figures_v2_postReview/figure_accuracy_per_celltype.png",
       height=7, width=10)



