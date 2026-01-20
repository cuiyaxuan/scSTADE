
conlabel<-function(hc1,k,true_label,compare=F){
    pbmc=CreateSeuratObject(counts = hc1, project = "HC_1")
    pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
    all.genes <- rownames(pbmc)
    mat<-as.matrix(pbmc[["RNA"]]@data)
    a <- VariableFeatures(pbmc)
    mat=mat[rownames(mat) %in% a,]
    # dim(mat)

    group1=read.csv("./label_4000.csv",row.names = 1)
    group2=read.csv("./label_4500.csv",row.names = 1)
    group3=read.csv("./label_5000.csv",row.names = 1)
    group4=read.table("./label1.txt",row.names = 1)
    group5=read.table("./label2.txt",row.names = 1)
    group6=read.table("./label3.txt",row.names = 1)

    file.remove("./label_4000.csv")
    file.remove("./label_4500.csv")
    file.remove("./label_5000.csv")
    file.remove("./label1.txt")
    file.remove("./label2.txt")
    file.remove("./label3.txt")

    group1=t(group1)
    group2=t(group2)
    group3=t(group3)
    group4=t(group4)
    group5=t(group5)
    group6=t(group6)
    mm=result(mat,group1,group2,group3,group4,group5,group6)

    label<-spectralClustering(mm, K=k)
    pre_label=label
    pre_label[1] 
    pre_label=as.data.frame(pre_label)
    rownames(pre_label)=colnames(mat)
    write.csv(pre_label,"label.csv")


    if(compare==TRUE){
        a <- rownames(pre_label)
        aa <- rownames(true_label)
        pre_label=pre_label[rownames(pre_label) %in% aa,]

        library("mclust")
        true_label=as.array(true_label[,1])
        ari=adjustedRandIndex(pre_label, true_label)
        paste("ARI:",ari)
        write.csv(ari,"ARI.csv")
    }
    
}

