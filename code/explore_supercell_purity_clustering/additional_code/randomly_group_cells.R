# Again reviewer.. wanting a comparison with random subsampling.
# I don't know how to calculate purity by randomly subsampling cells, but
# I can randomly group cells into supercells and compute purity.
#
# So for each gamma, count how many supercells you get, then randomly group
# the cells into those number of supercells then compute purity

set.seed(42)

library(data.table)
library(SuperCell)
library(ggplot2)

gamma <- seq(5, 50, 5)

annotations <- c(
    "data/explore_supercell_purity_clustering/levine_32dim/levine_32_asinh.csv",
    "data/explore_supercell_purity_clustering/samusik_all/samusik_all_asinh.csv"
)
datasets <- c("levine_32dim", "samusik_all")

random_assign_purity <- lapply(seq(1, length(annotations)), function(i) {
    dataset <- datasets[i]
    # get the annotation
    annotation <- fread(annotations[i])

    # get the map
    supercell_cell_map <- lapply(gamma, function(g) {
        x <- fread(paste0(
            "output/explore_supercell_purity_clustering/20230511/", dataset,
            "/supercell_runs/supercellCellMap_gamma",
            g, ".csv"))
        x$gamma <- g
        return(x)
    })
    supercell_cell_map <- rbindlist(supercell_cell_map)

    # Count number of supercells
    n_supercells <- supercell_cell_map[, .(n_supercell = uniqueN(SuperCellID)), by = gamma]

    n_cells <- nrow(annotation)

    random_alloc <- lapply(gamma, function(gam) {
        set.seed(42)
        n_supercell <- n_supercells[gamma == gam]$n_supercell
        ids <- data.table(
            CellId = annotation$CellId,
            random_id = sample(1:n_supercell, size = n_cells, replace = TRUE),
            gamma = gam
        )
        return(ids)
    })
    random_alloc <- rbindlist(random_alloc)


    # attach the annotation
    random_alloc_annot <- merge.data.table(
        random_alloc,
        annotation[, c("CellId", "Gated_Population")],
        by = "CellId"
    )

    fwrite(random_alloc_annot,
           paste0("output/explore_supercell_purity_clustering/20240103/random_assignment_",
                  dataset ,".csv"))
    # remove unassigned
    random_alloc_annot <- random_alloc_annot[Gated_Population != 'unassigned']

    # purity for random allocation
    random_alloc_purity <- lapply(gamma, function(gam) {
        dt <- random_alloc_annot[gamma == gam]
        purity <- supercell_purity(
            clusters = dt$Gated_Population,
            supercell_membership = dt$random_id
        )
        purity_dt <- data.table(
            random_id = names(purity),
            score = purity,
            gamma = gam
        )
        return(purity_dt)
    })
    random_alloc_purity <- rbindlist(random_alloc_purity)
    random_alloc_purity$dataset <- dataset

    return(random_alloc_purity)
})
random_assign_purity <- rbindlist(random_assign_purity)


ggplot(random_assign_purity, aes(x = factor(gamma), y = score)) +
    geom_violin() +
    stat_summary(fun = mean, geom = "point", size = 1, color = "red") +
    facet_wrap(~ dataset) +
    scale_y_continuous(breaks = pretty_breaks(n=10)) +
    labs(y = "Purity", x = "Gamma value",
         title = "Distribution of Purity in Randomly Assigned Groups",
         subtitle = "Number of groups is determined by the number of supercells generated for each gamma value")+
    theme_bw()

ggsave("~/Documents/supercell/manuscript/figures_v2_postReview/supp_figure_random_assignment.png",
       dpi=600,
       width=10,
       height=6)

fwrite(random_assign_purity, "output/explore_supercell_purity_clustering/20240103/random_assignment_purity.csv")

