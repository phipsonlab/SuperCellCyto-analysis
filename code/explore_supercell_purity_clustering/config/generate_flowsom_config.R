# Generate CSV file containing FlowSOM config for clustering.

metaclusters <- seq(15, 30, by = 5)
grid_size <- seq(10, 14, by = 1)

config <- expand.grid(metacluster = metaclusters, grid_size = grid_size)

write.csv(
  x = config,
  file = "code/levine_benchmarking_data/nextflow_pipeline/config/flowsom_config.csv",
  row.names = FALSE,
  quote = FALSE
)
