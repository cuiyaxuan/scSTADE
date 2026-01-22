# ==============================================================================
# File: function.R
# Description: 包含核心计算函数，供主程序调用
# ==============================================================================

# 1. 计算众数的辅助函数
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# 2. 核心计算函数 (cripar_opt)
# 优化说明：
# - 移除内部 library(philentropy)，改用 stats::dist
# - 输入改为向量 (gene_vec) 以减少内存传输
cripar_opt <- function(gene_vec, x_y_list, Ccri=50, highval=500, lowval=50){
  
  # 确保输入是数值型
  geneval <- as.numeric(gene_vec)
  
  # 1. 过滤表达量 > 0 的点
  valid_idx <- which(geneval > 0)
  
  # 如果有效点太少，直接返回0
  if(length(valid_idx) == 0) return(0)
  
  # 构建临时矩阵 A: [x, y, value]
  A <- cbind(x_y_list[valid_idx, ], val = geneval[valid_idx])
  b <- round(A[,3], 1)
  len_b <- length(b)
  
  # 初始化结果
  cri <- 0
  subset_A <- NULL
  
  # 2. 根据数量进行逻辑分支
  if(len_b <= Ccri){
    cri <- 0
  } else {
    if(len_b < highval && len_b > lowval){
      # 在阈值范围内，保留所有点
      subset_A <- A
    } else {
      # 超过高阈值，进行模式(Mode)过滤
      result <- getmode(b)
      
      if(floor(result) == ceiling(result)){
        # 整数情况 (区间: result-0.5 到 result+0.5)
        keep <- (A[,3] > (floor(result) - 0.5)) & (A[,3] < (ceiling(result) + 0.5))
      } else {
        # 非整数情况 (区间: floor 到 ceiling)
        keep <- (A[,3] > floor(result)) & (A[,3] < ceiling(result))
      }
      subset_A <- A[keep, , drop=FALSE]
    }
    
    # 3. 计算欧氏距离并评分
    if(!is.null(subset_A) && nrow(subset_A) > 0){
      # 使用 R 基础包 dist 计算欧氏距离，效率高且无需额外包
      cri_dat <- subset_A[, 1:2]
      dist_mat <- as.matrix(dist(cri_dat, method = "euclidean"))
      
      # 统计距离 <= 2 的点对数量
      # dist_mat 对称且对角线为0，需减去对角线数量并除以2
      count_le2 <- sum(dist_mat <= 2)
      n_rows <- nrow(subset_A)
      
      cri <- (count_le2 - n_rows) / 2
    }
  }
  
  return(cri)
}








