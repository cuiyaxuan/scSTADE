library(DEGman)
library("Seurat")
library("dplyr")
library("hdf5r")
library(philentropy)




distri<-function(hc1,label,n){
  current_directory <- getwd()


  all_files <- list.files(current_directory)


  csv_files <- all_files[grep("\\.csv$", all_files, ignore.case = TRUE)]


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

  amat=matlist[[n]]
  for (i in 1:k) {
    if (i == n) {
      next
    }
    amat=cbind(amat,matlist[[i]])
  }

  dim(amat)

  result=DEGman(amat,dim(matlist[[n]])[2],(dim(amat)[2]-dim(matlist[[n]])[2]))[,1]

  file_name <- paste0("output_", n, ".csv")
  write.csv(result,file_name,row.names=FALSE,col.names = TRUE)
  
  
}









