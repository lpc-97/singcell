library(dplyr)
library(Seurat)
library(cowplot)
library(patchwork)
library(monocle)
scRNA <- readRDS("scRNA.rds")
data<-as(as.matrix(scRNA@assays$RNA@counts),"sparseMatrix")
fData<-data.frame(gene_short_name=row.names(data),row.names=row.names(data))
pd<-new('AnnotatedDataFrame',data=scRNA@meta.data)
fd<-new('AnnotatedDataFrame',data=fData)
monocle_cods<-newCellDataSet(data,phenoData=pd,featureData=fd,lowerDetectionLimit=0.5,expressionFamily=negbinomial.size())#转化为monocle对象
HSMM<-monocle_cods
HSMM<-estimateSizeFactors(HSMM)
HSMM<-estimateDispersions(HSMM) #归一化
HSMM<-detectGenes(HSMM,min_expr = 0.1) #过滤基因
expressed_genes<-row.names(subset(fData(HSMM),num_cells_expressed>=10))
diff_test_res<-differentialGeneTest(HSMM[expressed_genes,],fullModelFormulaStr = "~seurat_clusters") #差异基因
ordering_genes<-row.names(subset(diff_test_res,qval<0.01)) #筛选出显著的差异基因 这里的q可以调节 q越小 定义轨迹的marker基因差异性越显著
HSMM<-setOrderingFilter(HSMM,ordering_genes)
HSMM<-reduceDimension(HSMM,max_components = 2,method='DDRTree') #降维
HSMM<-orderCells(HSMM)
plot_cell_trajectory(HSMM,color_by = "seurat_clusters")
plot_cell_trajectory(HSMM, color_by = "sample") + facet_wrap(~celltype_human, nrow = 1)
