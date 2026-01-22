library(DEGman)
library(Seurat)
library(dplyr)
library(hdf5r)
library(philentropy)
library(foreach)
library(doParallel)
library(doSNOW)

prep_matlist <- function(hc1, label) {
  pbmc <- CreateSeuratObject(counts = hc1, project = "HC_1", min.cells = 10)
  pbmc <- NormalizeData(pbmc,
                        normalization.method = "LogNormalize",
                        scale.factor = 100000)

  mat <- as.matrix(pbmc@assays$RNA@data)

  label <- as.numeric(as.vector(label[,1]))
  mat1 <- rbind(label, mat)
  mat1 <- t(mat1)
  mat1 <- mat1[order(mat1[,1]), ]

  k <- max(label)
  matlist <- vector("list", k)

  for (i in 1:k) {
    xx <- t(mat1[mat1[,1] == i, ])
    matlist[[i]] <- xx[-1, ]
  }

  return(matlist)
}

run_one_cluster <- function(matlist, n) {
  k <- length(matlist)

  amat <- matlist[[n]]
  for (i in 1:k) {
    if (i == n) next
    amat <- cbind(amat, matlist[[i]])
  }

  result <- DEGman(
    amat,
    dim(matlist[[n]])[2],
    dim(amat)[2] - dim(matlist[[n]])[2]
  )[, 1]

  file_name <- paste0("output_", n, ".csv")
  write.csv(result, file_name, row.names = FALSE)

  return(file_name)
}


# =========================
# 2. Functions
# =========================

# ---- 2.1 预处理函数（只运行一次） ----
prep_matlist <- function(hc1, label) {
  pbmc <- CreateSeuratObject(counts = hc1,
                             project = "HC_1",
                             min.cells = 10)

  pbmc <- NormalizeData(pbmc,
                        normalization.method = "LogNormalize",
                        scale.factor = 100000)

  mat <- as.matrix(pbmc@assays$RNA@data)

  label <- as.numeric(as.vector(label[, 1]))
  mat1 <- rbind(label, mat)
  mat1 <- t(mat1)
  mat1 <- mat1[order(mat1[, 1]), ]

  k <- max(label)
  matlist <- vector("list", k)

  for (i in 1:k) {
    xx <- t(mat1[mat1[, 1] == i, ])
    matlist[[i]] <- xx[-1, ]
  }

  return(matlist)
}

# ---- 2.2 单 cluster DEG 计算函数 ----
run_one_cluster <- function(matlist, n) {
  k <- length(matlist)

  amat <- matlist[[n]]
  for (i in 1:k) {
    if (i == n) next
    amat <- cbind(amat, matlist[[i]])
  }

  result <- DEGman(
    amat,
    dim(matlist[[n]])[2],
    dim(amat)[2] - dim(matlist[[n]])[2]
  )[, 1]

  file_name <- paste0("output_", n, ".csv")
  write.csv(result, file_name, row.names = FALSE)

  return(file_name)
}