library(Seurat)
library(monocle3)
Trans <- readRDS("~/SMC12.rds")
data <- GetAssayData(combinedc6n, assay = 'RNA', slot = 'counts')
cell_metadata <- combinedc6n@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds=preprocess_cds(cds,norm_method='log',method='PCA')#PCA
cds<-reduce_dimension(cds)#降维
cds.embed <- cds@int_colData$reducedDims$UMAP #选择用Seurat的降维图像
int.embed <- Embeddings(combinedc6n, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
cds <- cluster_cells(cds) #聚类
cds<-learn_graph(cds)
cds <- order_cells(cds) #手动选择root
plot_cells(cds, color_cells_by = "pseudotime") #查看轨迹
plot_cells(cds, color_cells_by = "pseudotime",cell_size=2,show_trajectory_graph = FALSE) #查看轨迹
