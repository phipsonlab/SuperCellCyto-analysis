library(data.table)
library(ggplot2)
library(scales)
library(viridis)

#---- samusik_all data ----

# get the annotation
annotation <- fread("data/explore_supercell_purity_clustering/samusik_all/samusik_all_asinh.csv")
annotation <- annotation[Gated_Population != 'unassigned']
annotation <- annotation[, c("CellId", "Gated_Population")]

purity_scores <- fread("output/explore_supercell_purity_clustering/20230511/samusik_all/evaluation/supercell_purities.csv")

# get purity less than 1.
# purity_filtered <- purity_scores[purity < 1]
purity_filtered <- purity_scores

# get the map
supercell_cell_map <- lapply(seq(5, 50, 5), function(g) {
    x <- fread(paste0(
        "output/explore_supercell_purity_clustering/20230511/samusik_all/supercell_runs/supercellCellMap_gamma",
        g, ".csv"))
    x$gamma <- g
    return(x)
})
supercell_cell_map <- rbindlist(supercell_cell_map)

supercell_cell_map_fil <- merge.data.table(
    supercell_cell_map,
    purity_filtered,
    by = c('SuperCellID', 'gamma')
)

supercell_cell_map_merged <- merge.data.table(
    supercell_cell_map_fil,
    annotation,
    by = "CellId"
)

# count how many different cell type per supercell
n_cell_type_per_supercell <- supercell_cell_map_merged[, .(n_cell_type = uniqueN(Gated_Population)),  by = c("SuperCellID", "gamma")]
n_cell_type_per_supercell <- merge.data.table(
    n_cell_type_per_supercell,
    purity_filtered,
    by = c("SuperCellID", "gamma")
)
# wrangle so count how many supercell with 2,3,4,etc. cell types
n_supercell_per_n_cell_type_samusik <- n_cell_type_per_supercell[, .(cnt = .N), by = c("gamma", "n_cell_type")]
n_supercell_per_n_cell_type_samusik$dataset <- "samusik_all"

#---- levine_32 data ----

# get the annotation
annotation <- fread("data/explore_supercell_purity_clustering/levine_32dim/levine_32_asinh.csv")
annotation <- annotation[Gated_Population != 'unassigned']
annotation <- annotation[, c("CellId", "Gated_Population")]

purity_scores <- fread("output/explore_supercell_purity_clustering/20230511/levine_32dim/evaluation/supercell_purities.csv")

# get purity less than 1.
# purity_filtered <- purity_scores[purity < 1]
purity_filtered <- purity_scores

# get the map
supercell_cell_map <- lapply(seq(5, 50, 5), function(g) {
    x <- fread(paste0(
        "output/explore_supercell_purity_clustering/20230511/levine_32dim/supercell_runs/supercellCellMap_gamma",
        g, ".csv"))
    x$gamma <- g
    return(x)
})
supercell_cell_map <- rbindlist(supercell_cell_map)

supercell_cell_map_fil <- merge.data.table(
    supercell_cell_map,
    purity_filtered,
    by = c('SuperCellID', 'gamma')
)

supercell_cell_map_merged <- merge.data.table(
    supercell_cell_map_fil,
    annotation,
    by = "CellId"
)

# count how many different cell type per supercell
n_cell_type_per_supercell <- supercell_cell_map_merged[, .(n_cell_type = uniqueN(Gated_Population)),  by = c("SuperCellID", "gamma")]
n_cell_type_per_supercell <- merge.data.table(
    n_cell_type_per_supercell,
    purity_filtered,
    by = c("SuperCellID", "gamma")
)

# wrangle so count how many supercell with 2,3,4,etc. cell types
n_supercell_per_n_cell_type_levine <- n_cell_type_per_supercell[, .(cnt = .N), by = c("gamma", "n_cell_type")]
n_supercell_per_n_cell_type_levine$dataset <- "levine_32dim"

#---- Add everything and plot ----
n_supercell_per_n_cell_type <- rbind(n_supercell_per_n_cell_type_samusik, n_supercell_per_n_cell_type_levine)
n_supercell_per_n_cell_type[, total := sum(cnt), by=.(gamma, dataset)]
n_supercell_per_n_cell_type[, proportion := cnt/total]
n_supercell_per_n_cell_type[, perc := round(100 * proportion, 2)]
n_supercell_per_n_cell_type[, perc_count := paste0(cnt, " (", perc, "%)")]

# Collapse everything greater than 10
n_supercell_per_n_cell_type[, n_cell_type := ifelse(n_cell_type > 9, ">=10", n_cell_type)]
n_supercell_per_n_cell_type[, n_cell_type := factor(n_cell_type, levels = c(seq(1, 9), ">=10"))]

ggplot(n_supercell_per_n_cell_type, aes(x = factor(gamma), fill = n_cell_type, y = perc)) +
    geom_histogram(stat="identity") +
    facet_wrap(~ dataset, scales = "free_y") +
    labs(
        title = "Number of cell types captured in supercells",
        y = "Percentage of supercells",
        fill = "Number of cell types",
        x = "Gamma"
    ) +
    scale_y_continuous(breaks = pretty_breaks(10)) +
    scale_fill_brewer(palette = "Paired") +
    theme_bw()

ggsave("~/Documents/supercell/manuscript/figures_v2_postReview/supp_figure_n_cells_per_supercell.png",
       dpi=600,
       width=10,
       height=5)







