
estimate_spatial<-function(hc1=hc1){
  library("Seurat")
  library("dplyr")
  library("hdf5r")
  pbmc=CreateSeuratObject(counts = hc1, project = "HC_1")
  pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 5000)
  all.genes <- rownames(pbmc)
  mat<-as.matrix(pbmc[["RNA"]]@data)
  a <- VariableFeatures(pbmc)
  str(a)
  mat=mat[rownames(mat) %in% a,]
  
  dim(mat)
  max(mat)
  mat=log2(mat+1)
  mat=t(mat)
  
  png("./estimate_plot.png", width = 800, height = 600)
  library(ClusterR)
  opt = Optimal_Clusters_KMeans(mat, max_clusters = 20, plot_clusters = T,
                                
                                criterion = 'distortion_fK', fK_threshold = 0.85,
                                
                                initializer = 'optimal_init', tol_optimal_init = 0.2)
  dev.off()
  
}





