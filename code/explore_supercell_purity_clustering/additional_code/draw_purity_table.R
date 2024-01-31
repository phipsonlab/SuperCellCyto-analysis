library(data.table)
library(ggplot2)
library(scales)

purity_scores <- lapply(c("samusik_all", "levine_32dim"), function(dataset) {
    x <- fread(paste0("output/explore_supercell_purity_clustering/20230511/", dataset, "/evaluation/supercell_purities.csv"))
    x$dataset <- dataset
    return(x)
})
purity_scores <- rbindlist(purity_scores)

# categorise purity score
purity_category <- function(score) {
  if (score == 1) {
    return("1")
  } else if (score >= 0.9) {
    return("0.9-1")
  } else if (score >= 0.5) {
    return("0.5-0.9")
  } else {
    return("Below_0.5")
  }
}

# Apply the function to categorise the purity score
purity_scores[, purity_category := sapply(purity, purity_category)]
purity_scores[, purity_category := factor(purity_category, levels = c("Below_0.5", "0.5-0.9", "0.9-1", "1"))]

purity_scores_cnt_per_cat <- purity_scores[, .(count = .N), by = .(gamma, dataset, purity_category)]

# Proportion per purity category
purity_scores_cnt_per_cat[, total := sum(count), by = .(gamma, dataset)]
purity_scores_cnt_per_cat[, proportion := count / total]
purity_scores_cnt_per_cat[, perc := round(100 * proportion, 2)]

# Ordering for easy viewing
purity_scores_cnt_per_cat <- purity_scores_cnt_per_cat[order(dataset, gamma, purity_category)]

# Pretty table?
purity_scores_cnt_per_cat[, perc_count := paste0(count, " (", perc, "%)")]

samusik_purity <- dcast(purity_scores_cnt_per_cat[dataset == 'samusik_all'],
                        gamma ~ purity_category,
                        value.var = "perc_count")
samusik_purity$dataset <- "samusik_all"
# fwrite(samusik_purity, "output/explore_supercell_purity_clustering/20240103/samusik_purity_category.csv")
levine_purity <- dcast(purity_scores_cnt_per_cat[dataset == 'levine_32dim'],
                        gamma ~ purity_category,
                        value.var = "perc_count")
levine_purity$dataset <- "levine_32dim"
# fwrite(levine_purity, "output/explore_supercell_purity_clustering/20240103/levine_32dim_purity_category.csv")
all_purity <- rbind(samusik_purity, levine_purity)
fwrite(all_purity, "output/explore_supercell_purity_clustering/20240103/purity_category.csv")

# how many get purity > 0.9?




