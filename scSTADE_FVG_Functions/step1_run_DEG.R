# =========================
# 1. Load libraries
# =========================
source("functions_DEG.R")
library(DEGman)
library(Seurat)
library(dplyr)
library(hdf5r)
library(philentropy)
library(foreach)
library(doParallel)
library(doSNOW)



# =========================
# 3. Load data
# =========================
hc1 <- Read10X_h5(
  "/home/cuiyaxuan/spatialLIBD/151673/151673_filtered_feature_bc_matrix.h5"
)

label <- read.csv(
  "/home/cuiyaxuan/metric_change/revise_R2/est_151673/conlabel.csv",
  header = TRUE,
  row.names = 1
)

# =========================
# 4. Preprocess (once)
# =========================
matlist <- prep_matlist(hc1, label)
k <- length(matlist)    # number of clusters

# =========================
# 5. Parallel + progress bar
# =========================
ncore <- min(k, parallel::detectCores() - 1)
cl <- makeCluster(ncore)
registerDoSNOW(cl)

pb <- txtProgressBar(min = 0, max = k, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

files <- foreach(
  n = 1:k,
  .packages = c("DEGman"),
  .options.snow = opts
) %dopar% {
  run_one_cluster(matlist, n)
}

close(pb)
stopCluster(cl)

# =========================
# 6. Done
# =========================
print(files)



