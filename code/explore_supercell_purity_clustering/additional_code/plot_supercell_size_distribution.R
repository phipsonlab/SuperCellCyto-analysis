# Bloody reviewer..
# want the supercell size distribution for samusik_all and levine_32dim

library(data.table)
library(scales)
library(ggplot2)

gamma <- seq(5,50,5)
supercell_map <- lapply(c("levine_32dim", "samusik_all"), function(dataset) {
    x <- lapply(gamma, function(gam) {
      x <- fread(paste0("output/explore_supercell_purity_clustering/20230511/",
                   dataset, "/supercell_runs/supercellCellMap_gamma",
                   gam, ".csv"))
      x$gamma <- gam
      return(x)
    })
    x <- rbindlist(x)
    x$dataset <- dataset
    return(x)
})
supercell_map <- rbindlist(supercell_map)

# Count number of cells per supercell and dataset
supercell_sizes <- supercell_map[, .(n_cells = .N), by = .(SuperCellID, dataset, gamma)]
# Log the size as the plot is super skewed
supercell_sizes[, log_n_cells := log10(n_cells)]

# Plot
ggplot(supercell_sizes, aes(x=factor(gamma), y=log_n_cells)) +
    geom_violin() +
    stat_summary(fun = mean, geom = "point", size = 1, color = "red") +
    facet_wrap(~ dataset, scales = "free_y") +
    scale_y_continuous(labels = function(x) format(10^x, scientific = FALSE)) +
    theme_bw() +
    labs(x = "Gamma", y = "Number of Cells", title = "Distribution of supercell sizes")
ggsave("~/Documents/supercell/manuscript/figures_v2_postReview/supp_figure_supercell_sizes.png",
       dpi=600,
       width=10,
       height=5)


