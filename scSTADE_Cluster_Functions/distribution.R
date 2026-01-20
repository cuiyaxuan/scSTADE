library(DEGman)
library("Seurat")
library("dplyr")
library("hdf5r")
library(philentropy)




distri<-function(hc1,label,n){
  current_directory <- getwd()

  # 获取当前目录下的所有文件
  all_files <- list.files(current_directory)

  # 筛选出所有CSV文件
  csv_files <- all_files[grep("\\.csv$", all_files, ignore.case = TRUE)]

  # 删除所有CSV文件
  for (csv_file in csv_files) {
    file_path <- file.path(current_directory, csv_file)
    unlink(file_path)
  }
  pbmc=CreateSeuratObject(counts = hc1, project = "HC_1", min.cells = 10)
  pbmc=NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 100000)

  mat<-pbmc@assays$RNA@data
  mat=as.matrix(mat)
  dim(mat)
  dim(label)
  label=t(label)
  mat1=rbind(label,mat)
  mat1=t(mat1)
  mat1<-mat1[order(mat1[,1]),]
  matlist=list()
  k=max(label)
    

  for (i in 1:k) {
    xx=t(mat1[mat1[,1]==i,])
    matlist[[i]]=xx[-1,]
  }



  
    library(foreach)
    library(doParallel)
    cl <- makeCluster((k-1))
    registerDoParallel(cl)
    foreach(i=1:k) %dopar% {
      if (i != n) {
        library(DEGman)
        amat=cbind(matlist[[n]],matlist[[i]])
        result <- DEGman(amat,dim(matlist[[n]])[2],dim(matlist[[i]])[2])[,1]
        #result=DEGman(amat,dim(matlist[[n]])[2],dim(matlist[[i]])[2])[,1]
        write.csv(result,paste(i,"_","result.csv", sep = ""),row.names=FALSE,col.names = TRUE)
      }
    }
    stopCluster(cl)
  
  
}









