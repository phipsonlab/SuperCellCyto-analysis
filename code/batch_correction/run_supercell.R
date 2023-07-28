library(Spectre)
library(here)
library(data.table)
library(stringr)
library(SuperCellCyto)
library(tictoc)
library(BiocParallel)
# library(parallel) # only needed to detectCores

supercell_out <- here(outdir, paste0(format(Sys.Date(), "%Y%d%m"), "_supercell_out"))
dir.create(supercell_out)

dat <- rbind(
    rbindlist(read.files(
        file.loc = here("output", "trussart_cytofruv", "Norm_Raw_RUV1b_RUV3b_simult", "RUV1b"),
        file.type = ".fcs"
    )),
    rbindlist(read.files(
        file.loc = here("output", "trussart_cytofruv", "Norm_Raw_RUV1b_RUV3b_simult", "RUV3b"),
        file.type = ".fcs"
    ))
)

# Remove unassigned
dat <- dat[! dat$FileName %in% c('Unassigned', 'Run3_Unassigned')]

dat$FileName <- paste0(dat$FileName, ".fcs")

metadata <- fread(here("data", "trussart_cytofruv", "metadata", "Metadata.csv"))

# Somehow, for some FCS file, we don't have the batch and sample information
# i.e. C3-D5 for both Run1 and Run3.
# We will abandon these files.
dat <- merge.data.table(
    x = dat,
    y = metadata,
    by.x = "FileName",
    by.y = "file_name"
)

# Keep the markers we have the metadata for
panel_metadata <- fread(here("data", "trussart_cytofruv", "metadata", "Panel.csv"))
panel_metadata[, concat := paste0(fcs_colname, "_", antigen)]

# Rename some weird symbols
names(dat) <- gsub("-", "_", names(dat))

markers_to_keep <- sapply(names(dat), function(col_name) {
    split_name <- str_split_1(col_name, "_")
    if (length(split_name) < 3 | split_name[3] == "Barcode")
        return(NA)
    return(paste(split_name[c(1, 3:length(split_name))], collapse = "_"))
})
markers_to_keep <- markers_to_keep[complete.cases(markers_to_keep)]
length(markers_to_keep) # must be at least 31, how many markers we have in the panel

setnames(dat, names(markers_to_keep), markers_to_keep)

cols_to_keep <- c(intersect(markers_to_keep, panel_metadata$concat),
                     "FileName", "sample_id", "condition", "patient_id", "batch"
)

# Remove columns we don't need
dat <- dat[, cols_to_keep, with = FALSE]

dat$CellID <- paste0("Cell_", seq(nrow(dat)))

markers <- cols_to_keep[1:31]

dat <- do.asinh(
    dat = dat,
    use.cols = markers,
    cofactor = 5,
    append.cf = TRUE
)

fwrite(dat, here(supercell_out, "cell_dat_asinh.csv"))

# Discard the untransformed markers to save RAM
cols_to_keep <- names(dat)[c(32:68)]
dat <- dat[, cols_to_keep, with = FALSE]

rm(metadata)
rm(panel_metadata)
gc()

n_rep <- 2

tic.clearlog()

markers <- paste0(markers, "_asinh_cf5")

for (rep in seq(n_rep)) {
    tic(paste0("supercell_rep", rep))

    BPPARAM <- MulticoreParam(workers = parallel::detectCores(), tasks = length(unique(dat$sample)))

    supercell_res <- runSuperCellCyto(
        dt = dat,
        markers = markers,
        sample_colname = "sample_id",
        cell_id_colname = "CellID",
        BPPARAM=BPPARAM,
        load_balancing = TRUE
    )

    toc(log = TRUE, quiet = TRUE)

    rm(BPPARAM)
    gc()

    fwrite(
        x = supercell_res$supercell_expression_matrix,
        file = here(supercell_out, paste0("supercellExpMat_rep", rep, ".csv"))
    )
    fwrite(
        x = supercell_res$supercell_cell_map,
        file = here(supercell_out, paste0("supercellCellMap_rep", rep, ".csv"))
    )
    saveRDS(
        object = supercell_res$supercell_object,
        file = here(supercell_out, paste0("supercellObj_rep", rep, ".rds"))
    )
}
# Export the supercell durations
writeLines(
    text = unlist(tic.log(format = TRUE)),
    con = here(supercell_out, "supercell_runtime.txt")
)

tic.clearlog()









