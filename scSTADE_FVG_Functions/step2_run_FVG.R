############################################################
## STEP 2: Compute FVG from saved DEG results
############################################################


source("functions_FVG.R")

# ==============================================================================
# File: main.R
# Description: 主执行脚本，包含数据IO和并行控制
# Usage: 运行此脚本前请确保 function.R 在同一目录
# ==============================================================================

library("Seurat")
library("dplyr")
library("hdf5r")
library("foreach")
library("parallel")
library("doParallel")

# 1. 加载功能函数
# # ------------------------------------------------------------------------------


# 2. 设置路径和读取数据
# ------------------------------------------------------------------------------
# 请根据实际环境修改以下路径
h5_path <- "/home/cuiyaxuan/spatialLIBD/151673/151673_filtered_feature_bc_matrix.h5"
spatial_path <- "/home/cuiyaxuan/spatialLIBD/151673/spatial/tissue_positions_list.csv"
gene_list_path <- "df1.csv"

cat("Loading data...\n")
hc1 <- Read10X_h5(h5_path)
tissue_local <- read.csv(spatial_path, row.names = 1, header = FALSE)

# 3. Seurat 预处理
# ------------------------------------------------------------------------------
pbmc <- CreateSeuratObject(counts = hc1, project = "HC_1", min.cells = 10)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 30000)

# 4. 数据对齐 (Data Alignment) - 内存优化版
# ------------------------------------------------------------------------------
# 找到公共细胞
common_cells <- intersect(colnames(pbmc), rownames(tissue_local))
cat(sprintf("Common cells found: %d\n", length(common_cells)))

# 过滤数据
pbmc <- subset(pbmc, cells = common_cells)
tissue_local <- tissue_local[common_cells, ]

# 读取目标基因
df1 <- read.csv(gene_list_path, header = TRUE)
target_genes <- df1[, 1] # 假设第一列是基因名

# 确保只计算存在的基因
valid_genes <- intersect(target_genes, rownames(pbmc))

# 提取表达矩阵 (行=细胞, 列=基因)
# 使用 GetAssayData 提取并转置，仅提取需要的基因
mat_sparse <- GetAssayData(pbmc, layer = "data")[valid_genes, ]
dat1 <- as.matrix(t(mat_sparse)) 

# 提取坐标 (假设 tissue_local 前两列为坐标，请根据实际情况调整列索引)
# 通常 tissue_local 的列可能是: row, col, imagerow, imagecol 等
# 这里保留原逻辑 dat[,3:4] 对应的位置，即 tissue_local 的前两列
x_y_list <- tissue_local[, 1:2]

if(!all(rownames(dat1) == rownames(x_y_list))) stop("Data alignment error!")

# 5. 并行计算 (Parallel Execution)
# ------------------------------------------------------------------------------
n_cores <- 24
cat(sprintf("Starting parallel processing with %d cores...\n", n_cores))

cls <- makeCluster(n_cores)
registerDoParallel(cls)

# 重要：将主环境中的函数和变量导出到所有从节点
clusterExport(cls, c("cripar_opt", "getmode", "x_y_list"))

# 并行循环
crinum <- foreach(gene_vec = iter(dat1, by = "col"), 
                  .combine = 'c', 
                  .packages = c("stats")) %dopar% {
                    # 调用 function.R 中的函数
                    cripar_opt(gene_vec, x_y_list, Ccri=10, highval=500, lowval=50)
                  }

stopCluster(cls)

# 6. 结果输出
# ------------------------------------------------------------------------------
# 整合结果
result_df <- data.frame(gene = valid_genes, crinum = crinum)

# 与原始列表合并 (保留原始顺序，未计算的为NA)
final_res <- merge(df1, result_df, by.x = names(df1)[1], by.y = "gene", all.x = TRUE)

# 排序并保存
fvg_sort <- final_res[order(-final_res$crinum), ]
write.csv(fvg_sort, "fvg.csv", row.names = FALSE)

cat("Done! Results saved to 'fvg.csv'.\n")


