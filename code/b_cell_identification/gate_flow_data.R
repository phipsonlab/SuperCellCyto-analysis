# Manually gate the flow cytometry data to include only live single cells.
# The gating strategy is as such:
# 1. FSC-H and FSC-A to isolate only the single events. (Also check SSC-H vs SSC-A).
# 2. FSC-A and SSC-A to remove debris.
# 3. Live/Dead and SSC-A to isolate live cells.

library(CytoExploreR)
library(tidyverse)

setwd("data/oetjen_bm_dataset/FR-FCM-ZYQ9_flow/B_cell_panel/")

gs <- cyto_setup(gatingTemplate = paste0(Sys.Date(), "_single_cell_gates.csv"))
gs <- cyto_compensate(gs)

# Just for gating live single cells, we only need to transform the Live/Dead
# marker for gating out the dead cells.
gs <- cyto_transform(gs, type = "logicle", channels = "V545-A")

cyto_gate_draw(gs,
  parent = "root",
  alias = "Single events",
  channels = c("FSC-A", "FSC-H")
)

cyto_gate_draw(gs,
  parent = "Single events",
  alias = "Cells",
  channels = c("FSC-A", "SSC-A")
)

cyto_gate_draw(gs,
  parent = "Cells",
  alias = "Live single cells",
  channels = c("LIVE/DEAD", "SSC-A")
)

# Check that other channels have not been transformed?
cyto_plot(gs, parent = "Live single cells", channels = c("LIVE/DEAD", "SSC-A"))

# 18,648,656 cells
all_cells_list <- cyto_extract(gs, raw = TRUE, parent = "root")
all_cells <- lapply(names(all_cells_list), function(samp_name) {
  dat_df <- data.frame(all_cells_list[[samp_name]])
  dat_df$sample <- str_split_i(samp_name, "_", 3)
  return(dat_df)
}) |> bind_rows()
dim(all_cells)

# 8,314,260 cells
single_live_cells_list <- cyto_extract(gs, raw = TRUE, parent = "Live single cells")
single_live_cells <- lapply(names(single_live_cells_list), function(samp_name) {
  dat_df <- data.frame(single_live_cells_list[[samp_name]])
  dat_df$sample <- str_split_i(samp_name, "_", 3)
  return(dat_df)
}) |> bind_rows()
dim(single_live_cells)

# Rename the channel to markers
markers <- read_csv("040423-Experiment-Markers.csv")

markers_for_renaming <- markers |>
  filter(!is.na(marker) & marker != "LIVE/DEAD") |>
  mutate(
    ch_name = gsub("[^[:alnum:] ]", ".", channel),
    marker_name = gsub("[^[:alnum:] ]", "_", ifelse(is.na(marker), channel, marker))
  ) |>
  pull(ch_name, name = marker_name)

single_live_cells <- single_live_cells |>
  select(all_of(c("sample", markers_for_renaming)))

all_cells <- all_cells |>
  select(all_of(c("sample", markers_for_renaming)))

dim(single_live_cells)
dim(all_cells)

# For saving intermediate data
# TODO: manually move the data out to output/oetjen_b_cell_panel folder.
write_csv(single_live_cells, paste0(Sys.Date(), "_single_live_cells.csv"))
write_csv(all_cells, paste0(Sys.Date(), "_all_cells.csv"))


# Visualise the single live cells gate for manuscript.
cyto_plot_save("single_live_cells.png", width = 12, height = 5)
cyto_plot_gating_scheme(gs, header = "Single live cells")
cyto_plot_complete()

# Test the single live cells to see if they are intact
dim(read_csv(paste0(Sys.Date(), "_single_live_cells.csv")))
