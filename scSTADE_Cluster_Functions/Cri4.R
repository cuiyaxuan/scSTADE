library(SingleCellExperiment)
library(SC3)
library(ggplot2)
library("Seurat")
library("ggplot2")
library("dplyr")
library("hdf5r")
library(foreach)
library(parallel)
library(doParallel)




select_feature<-function(hc1,initial,rate){
  feature_dat=vector(mode = 'list', length = 3)
  for (i in 1:3) {
    pbmc=CreateSeuratObject(counts = hc1, project = "HC_1")
    pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = initial)
    all.genes <- rownames(pbmc)
    mat<-as.matrix(pbmc[["RNA"]]@data)
    a <- VariableFeatures(pbmc)
    mat=mat[rownames(mat) %in% a,]
    feature_dat[[i]] = mat
    initial=initial+rate
  }
  return(feature_dat)
}




pearson_metric<-function(X,feature,k){
  library(SingleCellExperiment)
  library(SC3)
  library(foreach)
  library(parallel)
  library(doParallel)
  pearson_pre <- function(mat,k){
    library(SingleCellExperiment)
    library(SC3)
    sce <- SingleCellExperiment(
      assays = list(
        counts = as.matrix(mat),
        logcounts = log2(as.matrix(mat) + 1)
      )
    )
    rowData(sce)$feature_symbol <- rownames(sce)
    sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]
    sce <- sc3_prepare(sce)
    sce <- sc3_calc_dists(sce)
    metadata(sce)$sc3$distances[[1]]=metadata(sce)$sc3$distances[[2]]
    metadata(sce)$sc3$distances[[3]]=metadata(sce)$sc3$distances[[2]]
    sce <- sc3_calc_transfs(sce)
    sce <- sc3_kmeans(sce, ks = k)
    sce <- sc3_calc_consens(sce)
    return(sce@metadata$sc3$consensus[[1]]$consensus)
  }
  
  if(X==1){
    pearson_dat1=pearson_pre(feature[[1]],k)
    write.table(pearson_dat1,file="data1.txt",sep=" ",quote=TRUE)
  }
  if(X==2){
    pearson_dat2=pearson_pre(feature[[2]],k)
    write.table(pearson_dat2,file="data2.txt",sep=" ",quote=TRUE)
  }
  if(X==3){
    pearson_dat3=pearson_pre(feature[[3]],k)
    write.table(pearson_dat3,file="data3.txt",sep=" ",quote=TRUE)
  }
  # if(X==4){
  #   pearson_dat4=pearson_pre(feature[[4]],k)
  #   write.table(pearson_dat4,file="data4.txt",sep=" ",quote=TRUE)
  # }
  # if(X==5){
  #   pearson_dat5=pearson_pre(feature[[5]],k)
  #   write.table(pearson_dat5,file="data5.txt",sep=" ",quote=TRUE)
  # }
}



#View(feature[[1]])

construct_adj_matrix<-function(feature,tissue_local){
  mat=t(feature)
  aa <- rownames(mat)
  tissue_local=tissue_local[rownames(tissue_local) %in% aa,]
  DF1 <- mutate(tissue_local, id = rownames(tissue_local))
  class(tissue_local)
  class(mat)
  mat=as.data.frame(mat)
  DF2 <- mutate(mat, id = rownames(mat))
  dat=merge(DF1,DF2,by="id")
  dat1=dat[,1:4]
  hyperG=matrix(0,dim(mat)[1],dim(mat)[1])
  for (i in 1:dim(mat)[1]) {
    x = dat1[i,3]
    y = dat1[i,4]
    for (j in 1:dim(mat)[1]) {
      x1 = dat1[j,3]
      y1 = dat1[j,4]
      radius=(x-x1)^2+(y-y1)^2
      if(radius<=16)
        hyperG[i,j]=1
    }
  }
  return(hyperG)
}






library(foreach)
library(parallel)
library(doParallel)

spectral_nei<-function(X,K){
  W1=read.table("./data1.txt",header = T,quote = "",sep=' ')
  W2=read.table("./data2.txt",header = T,quote = "",sep=' ')
  W3=read.table("./data3.txt",header = T,quote = "",sep=' ')
  hyperG=read.table("./adj_matrix.txt",header = T,quote = "",sep=' ')
  file.remove("./data1.txt")
  file.remove("./data2.txt")
  file.remove("./data3.txt")
  file.remove("./adj_matrix.txt")
  iter<-function(hyperG,W,K){
    hyperG=hyperG*0.9
    hyperG[hyperG==0]=0.1
    spectralClustering <- function(affinity, K, type=3) {
      
      ###This function implements the famous spectral clustering algorithms. There are three variants. The default one is the third type. 
      ###THe inputs are as follows:
      
      #affinity: the similarity matrix;
      #K: the number of clusters
      # type: indicators of variants of spectral clustering 
      
      d = rowSums(affinity)
      d[d == 0] = .Machine$double.eps
      D = diag(d)
      L = D - affinity
      if (type == 1) {
        NL = L
      } else if (type == 2) {
        Di = diag(1 / d)
        NL = Di %*% L
      } else if(type == 3) {
        Di = diag(1 / sqrt(d))
        NL = Di %*% L %*% Di
      }
      eig = eigen(NL)
      res = sort(abs(eig$values),index.return = TRUE)
      U = eig$vectors[,res$ix[1:K]]
      normalize <- function(x) x / sqrt(sum(x^2))
      if (type == 3) {
        U = t(apply(U,1,normalize))
      }
      eigDiscrete = .discretisation(U)
      eigDiscrete = eigDiscrete$discrete
      labels = apply(eigDiscrete,1,which.max)
      
      
      
      return(labels)
    }
    
    
    .discretisation <- function(eigenVectors) {
      
      normalize <- function(x) x / sqrt(sum(x^2))
      eigenVectors = t(apply(eigenVectors,1,normalize))
      
      n = nrow(eigenVectors)
      k = ncol(eigenVectors)
      
      R = matrix(0,k,k)
      R[,1] = t(eigenVectors[round(n/2),])
      
      mini <- function(x) {
        i = which(x == min(x))
        return(i[1])
      }
      
      c = matrix(0,n,1)
      for (j in 2:k) {
        c = c + abs(eigenVectors %*% matrix(R[,j-1],k,1))
        i = mini(c)
        R[,j] = t(eigenVectors[i,])
      }
      
      lastObjectiveValue = 0
      for (i in 1:20) {
        eigenDiscrete = .discretisationEigenVectorData(eigenVectors %*% R)
        
        svde = svd(t(eigenDiscrete) %*% eigenVectors)
        U = svde[['u']]
        V = svde[['v']]
        S = svde[['d']]
        
        NcutValue = 2 * (n-sum(S))
        if(abs(NcutValue - lastObjectiveValue) < .Machine$double.eps) 
          break
        
        lastObjectiveValue = NcutValue
        R = V %*% t(U)
        
      }
      
      return(list(discrete=eigenDiscrete,continuous =eigenVectors))
    }
    
    .discretisationEigenVectorData <- function(eigenVector) {
      
      Y = matrix(0,nrow(eigenVector),ncol(eigenVector))
      maxi <- function(x) {
        i = which(x == max(x))
        return(i[1])
      }
      j = apply(eigenVector,1,maxi)
      Y[cbind(1:nrow(eigenVector),j)] = 1
      
      return(Y)
      
    }
    
    .dominateset <- function(xx,KK=20) {
      ###This function outputs the top KK neighbors.	
      
      zero <- function(x) {
        s = sort(x, index.return=TRUE)
        x[s$ix[1:(length(x)-KK)]] = 0
        return(x)
      }
      normalize <- function(X) X / rowSums(X)
      A = matrix(0,nrow(xx),ncol(xx));
      for(i in 1:nrow(xx)){
        A[i,] = zero(xx[i,]);
        
      }
      
      
      return(normalize(A))
    }
    # Calculate the mutual information between vectors x and y.
    .mutualInformation <- function(x, y) {
      classx <- unique(x)
      classy <- unique(y)
      nx <- length(x)
      ncx <- length(classx)
      ncy <- length(classy)
      
      probxy <- matrix(NA, ncx, ncy)
      for (i in 1:ncx) {
        for (j in 1:ncy) {
          probxy[i, j] <- sum((x == classx[i]) & (y == classy[j])) / nx
        }
      }
      
      probx <- matrix(rowSums(probxy), ncx, ncy)
      proby <- matrix(colSums(probxy), ncx, ncy, byrow=TRUE)
      result <- sum(probxy * log(probxy / (probx * proby), 2), na.rm=TRUE)
      return(result)
    }
    
    # Calculate the entropy of vector x.
    .entropy <- function(x) {
      class <- unique(x)
      nx <- length(x)
      nc <- length(class)
      
      prob <- rep.int(NA, nc)
      for (i in 1:nc) {
        prob[i] <- sum(x == class[i])/nx
      }
      
      result <- -sum(prob * log(prob, 2))
      return(result)
    }
    
    .repmat = function(X,m,n){
      ##R equivalent of repmat (matlab)
      if (is.null(dim(X))) {
        mx = length(X)
        nx = 1
      } else {
        mx = dim(X)[1]
        nx = dim(X)[2]
      }
      matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
    }
    
    normalize <- function(X) X / rowSums(X)
    W=normalize(W)
    W = (W+t(W))/2
    hyperG=normalize(hyperG)
    W=as.matrix(hyperG) %*% as.matrix(W) %*% as.matrix(t(hyperG))
    W = (W+t(W))/2
    label<-spectralClustering(W, K)
    return(label)
  }
  
  
  
  if(X==1){
    label1=iter(hyperG,W1,K)
    write.table(label1,file="label1.txt",sep=" ",quote=TRUE)
  }
  if(X==2){
    label2=iter(hyperG,W2,K)
    write.table(label2,file="label2.txt",sep=" ",quote=TRUE)
  }
  if(X==3){
    label3=iter(hyperG,W3,K)
    write.table(label3,file="label3.txt",sep=" ",quote=TRUE)
  }
  
}












