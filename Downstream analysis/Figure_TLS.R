library(Seurat)
library(ggplot2)
library(patchwork)
library(clusterProfiler)
library(dplyr)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(CellChat)

#setwd("/Users/cuiyang/Desktop/R/spatial")
plotColor <- c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
               "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB",
               "#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4",
               "#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")

human_breast <- Load10X_Spatial('/Users/cyx/spatialLIBD/3.human_breast_cancer/') # your path


plot1 <- VlnPlot(human_breast, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(human_breast, features = "nCount_Spatial") + theme(legend.position = "right")
# wrap_plots(plot1, plot2)

human_breast <- SCTransform(human_breast, assay = "Spatial", verbose = FALSE)

human_breast <- human_breast %>% 
  RunPCA(assay = "SCT", verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  RunUMAP(reduction = "pca", dims = 1:30)

metadata <- read.table('/Users/cyx/spatialLIBD/3.human_breast_cancer/label.csv', sep = ',', header = TRUE, row.names = 1)

human_breast <- AddMetaData(human_breast, metadata)

#-----------FigA----------------
domain20_col <- plotColor
names(domain20_col) <- names(table(human_breast$domain30))
p2 <- SpatialDimPlot(human_breast, group.by = 'domain30', label = TRUE, cols = domain20_col) + NoLegend()
p2


#-----------FigB----------------
# SpatialFeaturePlot(human_breast, features = c('HLA-DPA1', 'HLA-DMA'), pt.size.factor = 1.6, ncol = 2, crop = TRUE)
# SpatialFeaturePlot(human_breast, features = c('PTPRC', 'CD3D'), pt.size.factor = 1.6, ncol = 2, crop = TRUE)
# SpatialFeaturePlot(human_breast, features = c('CD79A', 'CD68'), pt.size.factor = 1.6, ncol = 2, crop = TRUE)
# SpatialFeaturePlot(human_breast, features = c('PTPRC', 'CD3D', 'CD79A', 'CD68'), pt.size.factor = 1.6, ncol = 2, crop = TRUE)

SpatialFeaturePlot(human_breast, features = c('PTPRC'))
SpatialFeaturePlot(human_breast, features = c('CD3D'))
SpatialFeaturePlot(human_breast, features = c('CD79A'))
SpatialFeaturePlot(human_breast, features = c('CD68'))
SpatialFeaturePlot(human_breast, features = c('ITGAX'))
SpatialFeaturePlot(human_breast, features = c('HLA-DPA1'))
SpatialFeaturePlot(human_breast, features = c('HLA-DMA'))
#-----------FigC----------------
cluster13_markers <- FindMarkers(human_breast, group.by = 'domain30',ident.1 = '13')
cluster13_markers <- cluster13_markers[cluster13_markers$p_val_adj < 0.05,]
gene.df <- bitr(rownames(cluster13_markers),fromType="SYMBOL",toType="ENTREZID", OrgDb = org.Hs.eg.db)#Convert to ENTREZID    
gene <- gene.df$ENTREZID
ego_CC <- enrichGO(gene = gene,
                   OrgDb=org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "CC",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = TRUE)

ego_BP <- enrichGO(gene = gene,
                   OrgDb=org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "BP",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = TRUE)

ego_MF <- enrichGO(gene = gene,
                   OrgDb=org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "MF",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = TRUE)
ego_result_BP <- as.data.frame(ego_BP)
ego_result_CC <- as.data.frame(ego_CC)
ego_result_MF <- as.data.frame(ego_MF)
display_number = c(10, 10, 10)#show top 10 enrichment result
ego_result_BP <- as.data.frame(ego_BP)[1:display_number[1], ]
ego_result_CC <- as.data.frame(ego_CC)[1:display_number[2], ]
ego_result_MF <- as.data.frame(ego_MF)[1:display_number[3], ]

ego_result_BP$pvalue <- -log10(ego_result_BP$pvalue)#-log10 transfer
ego_result_CC$pvalue <- -log10(ego_result_CC$pvalue)
ego_result_MF$pvalue <- -log10(ego_result_MF$pvalue)

go_enrich_df <- data.frame(
  ID=c(ego_result_BP$ID, ego_result_CC$ID, ego_result_MF$ID),                         
  Description=c(ego_result_BP$Description,ego_result_CC$Description,ego_result_MF$Description),
  GeneNumber=c(ego_result_BP$pvalue, ego_result_CC$pvalue, ego_result_MF$pvalue),
  type=factor(c(rep("biological process", display_number[1]), 
                rep("cellular component", display_number[2]),
                rep("molecular function", display_number[3])), 
              levels=c("biological process", "cellular component","molecular function" )))

go_enrich_df$type_order=factor(rev(as.integer(rownames(go_enrich_df))),labels=rev(go_enrich_df$Description))
COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")

ggplot(data=go_enrich_df, aes(x=type_order,y=GeneNumber, fill=type)) + 
  geom_bar(stat="identity", width=0.8) + 
  scale_fill_manual(values = COLS) + 
  coord_flip() + 
  xlab("GO term") + 
  ylab("-log10(P_value)") + 
  labs(title = "The Most Enriched GO Terms")+
  theme_bw()

#-----------FigD----------------
metadata <- human_breast@meta.data
domain30_13 <- subset(metadata, metadata$domain30 == '13')
table(domain30_13$ground_truth)
domain30_13_ids <- rownames(domain30_13)
metadata[domain30_13_ids,]$ground_truth <- 'TLS' 
human_breast@meta.data <- metadata
names(domain20_col) <- names(table(human_breast$ground_truth))
SpatialDimPlot(human_breast, group.by = 'ground_truth', label = TRUE, cols = domain20_col) + NoLegend()
#cellchat
visium.brain <- human_breast
# Prepare input data for CelChat analysis
data.input = GetAssayData(visium.brain, slot = "data", assay = "SCT") # normalized data matrix
meta <- data.frame(labels = metadata$ground_truth, row.names = rownames(metadata))
unique(meta$labels) # check the cell labels
# Spatial locations of spots from full (NOT high/low) resolution images are required
spatial.locs = GetTissueCoordinates(visium.brain, scale = NULL, cols = c("imagerow", "imagecol")) 
# Scale factors and spot diameters of the full resolution images 
scale.factors = jsonlite::fromJSON(txt = file.path("/Users/cyx/spatialLIBD/3.human_breast_cancer/spatial", 'scalefactors_json.json'))
scale.factors = list(spot.diameter = 65, spot = scale.factors$spot_diameter_fullres, # these two information are required
                     fiducial = scale.factors$fiducial_diameter_fullres, hires = scale.factors$tissue_hires_scalef, lowres = scale.factors$tissue_lowres_scalef # these three information are not required
)

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels",
                           datatype = "spatial", coordinates = spatial.locs, scale.factors = scale.factors)
cellchat


CellChatDB <- CellChatDB.human # use CellChatDB.human if running on human data

CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchat@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multicore", workers = 4) # do parallel

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, 
                              distance.use = TRUE, interaction.length = 200, scale.distance = 0.01)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = rowSums(cellchat@net$count), weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = rowSums(cellchat@net$weight), weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
#-----------FigE----------------
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", height = 30)
netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", height = 30) 








