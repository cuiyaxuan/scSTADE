############################################################
## STEP 2: Compute FVG from saved DEG results
############################################################


source("functions_FVG.R")
library(tidyverse)

rounded_number <- 3

fre <- list.files(
  path = ".",
  pattern = "\\.csv$",
  full.names = TRUE
) %>%
  map_dfr(read.csv) %>%        #  rbind
  count(x) %>%                 # table(vec1$x)
  filter(n >= rounded_number) %>%
  pull(x)

write.csv(fre, "df1.csv", row.names = FALSE)

# ==============================================================================
# File: main.R
# Description: Main execution script for FVG calculation
# Usage: Make sure functions_FVG.R is in the same directory before running
# ==============================================================================

library(Seurat)
library(dplyr)
library(hdf5r)
library(foreach)
library(parallel)
library(doParallel)

# ==============================================================================
# 1. Set paths and load data
# ==============================================================================
h5_path <- "/home/cuiyaxuan/spatialLIBD/151673/151673_filtered_feature_bc_matrix.h5"
spatial_path <- "/home/cuiyaxuan/spatialLIBD/151673/spatial/tissue_positions_list.csv"
gene_list_path <- "df1.csv"

cat("Loading data...\n")
hc1 <- Read10X_h5(h5_path)
tissue_local <- read.csv(spatial_path, row.names = 1, header = FALSE)

# ==============================================================================
# 2. Seurat preprocessing
# ==============================================================================
pbmc <- CreateSeuratObject(counts = hc1, project = "HC_1", min.cells = 10)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 30000)

# ==============================================================================
# 3. Data alignment
# ==============================================================================
common_cells <- intersect(colnames(pbmc), rownames(tissue_local))
cat(sprintf("Common cells found: %d\n", length(common_cells)))

pbmc <- subset(pbmc, cells = common_cells)
tissue_local <- tissue_local[common_cells, ]

# Load target gene list
df1 <- read.csv(gene_list_path, header = TRUE)
target_genes <- df1[, 1]

# Keep genes present in expression matrix
valid_genes <- intersect(target_genes, rownames(pbmc))

# Extract expression matrix (cells × genes)
mat_sparse <- GetAssayData(pbmc, layer = "data")[valid_genes, ]
dat1 <- as.matrix(t(mat_sparse))

# Extract spatial coordinates (first two columns)
x_y_list <- tissue_local[, 1:2]

if (!all(rownames(dat1) == rownames(x_y_list))) {
  stop("Data alignment error!")
}

# ==============================================================================
# 4. Parallel computation
# ==============================================================================
n_cores <- 24
cat(sprintf("Starting parallel processing with %d cores...\n", n_cores))

cls <- makeCluster(n_cores)
registerDoParallel(cls)

# Export required objects to worker nodes
clusterExport(cls, c("cripar_opt", "getmode", "x_y_list"))

crinum <- foreach(
  gene_vec = iter(dat1, by = "col"),
  .combine = "c",
  .packages = "stats"
) %dopar% {
  cripar_opt(gene_vec, x_y_list, Ccri = 10, highval = 500, lowval = 50)
}

stopCluster(cls)

# ==============================================================================
# 5. Save results
# ==============================================================================
result_df <- data.frame(gene = valid_genes, crinum = crinum)

# Merge with original gene list (keep original order)
final_res <- merge(df1, result_df,
                   by.x = names(df1)[1],
                   by.y = "gene",
                   all.x = TRUE)

fvg_sort <- final_res[order(-final_res$crinum), ]
write.csv(fvg_sort, "fvg.csv", row.names = FALSE)

cat("Done! Results saved to 'fvg.csv'.\n")




