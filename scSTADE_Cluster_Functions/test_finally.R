library("Seurat")
library("dplyr")
library("hdf5r")



###################################################################################

library(foreach)
library(parallel)
library(doParallel)

cripar<-function(q,dat1,x_y_list, Ccri=50, highval=500,lowval=50){
  # Create the function.
  getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }
  # Create the vector with numbers.
  geneval=as.numeric(dat1[,q])
  gene1=cbind(x_y_list,geneval)
  A=gene1[gene1[,3]>0,]
  b=round(A[,3],1)
  # Calculate the mode using the user function.
  if(length(b)<=Ccri){
    cri=0
  }
  else if(length(b)<highval && length(b)>lowval){
    library(philentropy)
    cri_dat=A[,1:2]
    xx=distance(cri_dat, method="euclidean")
    sum(xx<=2)
    dim(A)[1]
    cri=(sum(xx<=2)-dim(A)[1])/2
  }
  else{
    result <- getmode(b)
    if(floor(result)==ceiling(result)){
      A=A[A[,3]>(floor(result)-0.5),]
      A=A[A[,3]<(ceiling(result)+0.5),]
    }
    else{
      A=A[A[,3]>floor(result),]
      A=A[A[,3]<ceiling(result),]
    }
    library(philentropy)
    cri_dat=A[,1:2]
    xx=distance(cri_dat, method="euclidean")
    sum(xx<=2)
    dim(A)[1]
    cri=(sum(xx<=2)-dim(A)[1])/2
  }
  return(cri)
}

