library(data.table)
library(ggplot2)
library(scales)
library(RColorBrewer)

metrics <- fread("output/trussart_cytofruv/20240111/bio_conservation_metrics.csv")
metrics[, Metric := ifelse(V1=="silhouette", "ASW_label", V1)]
metrics[, V1 := NULL]

metrics_melt <- melt(metrics, id.vars = "Metric", variable.name = "Method", value.name = "Score")
metrics_melt$Metric <- factor(metrics_melt$Metric, levels = c("ASW_label", "NMI", "ARI"))
metrics_melt$Method <- factor(metrics_melt$Method, levels = c("cycombine", "cytofruv", "uncorrected"))

ggplot(metrics_melt, aes(x = Method, y = Metric, fill = Score)) +
    geom_tile(color = "black") +
    geom_text(aes(label = round(Score, 3)), color="white") +
    scale_fill_distiller(palette = "Blues", direction = +1, limits = c(0, 1)) +
    theme_minimal()

ggsave("~/Documents/supercell/manuscript/figures_v2_postReview/batch_correction/fig_metrics.png",
       dpi = 600,
       bg="white")


