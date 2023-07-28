BiocManager::install("CATALYST", version = "3.10", force = TRUE)

# title: "CytofRUV_preprocessing_dataset"
# author: "Marie Trussart"
# date: "31/03/2020"

library(readxl)
library(flowCore)
library(matrixStats)
library(ggplot2)
library(ggridges)
library(reshape2)
library(dplyr)
library(limma)
library(ggrepel)
library(FlowSOM)
library(cowplot)
library(lme4)
library(multcomp)
library(grid)
library(gridExtra)
library(CATALYST)

# Define work directory that contains the raw data
wd <- "/vast/projects/phipson_supercell/trussart_cytofruv/"

# Load the file containing all the fcs files
file_list <- read.table(paste(wd, "Supp/all_files.txt", sep = ""), stringsAsFactors = F)
full_path_files <- as.vector(paste(wd, "fcs_files/", file_list$V1, sep = ""))

# Load all the files from the RUV1b run:
ff1 <- concatFCS(full_path_files[1:12])

# Load all the files from the RUV3B run:
ff2 <- concatFCS(full_path_files[13:24])

common_markers <- intersect(colnames(ff1), colnames(ff2))
expr_ff2 <- exprs(ff2)
expr_ff2 <- expr_ff2[, common_markers]
expr_ff1 <- exprs(ff1)
expr_ff1 <- expr_ff1[, common_markers]
expr_ffall <- rbind(expr_ff1, expr_ff2)
ffall <- ff1
exprs(ffall) <- expr_ffall
es <- flowCore::exprs(ffall)
es_t <- asinh(es / 5)
chs <- flowCore::colnames(ffall)

# ---- Custom functions ----
# The script modified some functions that were in catalyst. 
# Source the functions only after running the code above as it depends on 
# the variables that are only assigned using the code above.
source(paste0(wd, "prepare_data_customFunc.R"))


# ---- Run bead-based normalisation ----
### Normalise the data using Finck normalisation on both runs together
out_dir <- paste(wd, "Norm_Raw_RUV1b_RUV3b_simult/", sep = "")
dir.create(out_dir)

norm_ffall <- normCytof2(
    x = ffall, 
    y = "dvs", 
    k = 80, 
    remove = F, 
    out_path = out_dir
)

norm_ffall_without_beads <- normCytof2(
    x = ffall, 
    y = "dvs", 
    k = 80, 
    out_path = NULL, 
    remove = T
)

to_be_removed <- norm_ffall_without_beads[[3]]
exp_toberem <- exprs(to_be_removed)

### Debarcode Run1
norm_ff1 <- norm_ffall[1:dim(ff1)[1], ]
exp_norm_ff1 <- exprs(norm_ff1)
exp_norm_ff1_without_beads <- anti_join(as.data.frame(exp_norm_ff1), as.data.frame(exp_toberem))
exprs(norm_ff1) <- as.matrix(exp_norm_ff1_without_beads)
ind <- dim(ff1)[1] + 1
data(sample_key)
re0 <- assignPrelim(x = norm_ff1, y = sample_key, verbose = TRUE)
re0
# estimate separation cutoffs
re <- estCutoffs(x = re0)
# use global separation cutoff
applyCutoffs(x = re, sep_cutoffs = 0.35)
dir.create(paste0(wd, "Norm_Raw_RUV1b_RUV3b_simult/RUV1b/"), recursive = TRUE)
outFCS(x = re, y = norm_ff1, out_path = paste(wd, "Norm_Raw_RUV1b_RUV3b_simult/RUV1b/", sep = ""))


### Debarcode Run2
norm_ff2 <- norm_ffall[ind:dim(norm_ffall)[1], ]
exp_norm_ff2 <- exprs(norm_ff2)
exp_norm_ff2_without_beads <- anti_join(as.data.frame(exp_norm_ff2), as.data.frame(exp_toberem))
exprs(norm_ff2) <- as.matrix(exp_norm_ff2_without_beads)
data(sample_key)
re0 <- assignPrelim(x = norm_ff2, y = sample_key, verbose = TRUE)
re0
# estimate separation cutoffs
re <- estCutoffs(x = re0)
# use global separation cutoff
applyCutoffs(x = re, sep_cutoffs = 0.35)
dir.create(paste0(wd, "Norm_Raw_RUV1b_RUV3b_simult/RUV3b/"), recursive = TRUE)
outFCS(x = re, y = norm_ff2, out_path = paste(wd, "Norm_Raw_RUV1b_RUV3b_simult/RUV3b/", sep = ""))

# Rename RUV3b files with Run3_ prefix
setwd(wd)
setwd("Norm_Raw_RUV1b_RUV3b_simult/RUV3b")

for (f in list.files()) {
    file.rename(from = f, to = paste0("Run3_", f))
}


