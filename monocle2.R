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
monocle_cods<-newCellDataSet(data,phenoData=pd,featureData=fd,lowerDetectionLimit=0.5,expressionFamily=negbinomial.size())#ת��Ϊmonocle����
HSMM<-monocle_cods
HSMM<-estimateSizeFactors(HSMM)
HSMM<-estimateDispersions(HSMM) #��һ��
HSMM<-detectGenes(HSMM,min_expr = 0.1) #���˻���
expressed_genes<-row.names(subset(fData(HSMM),num_cells_expressed>=10))
diff_test_res<-differentialGeneTest(HSMM[expressed_genes,],fullModelFormulaStr = "~seurat_clusters") #�������
ordering_genes<-row.names(subset(diff_test_res,qval<0.01)) #ɸѡ�������Ĳ������ �����q���Ե��� qԽС ����켣��marker���������Խ����
HSMM<-setOrderingFilter(HSMM,ordering_genes)
HSMM<-reduceDimension(HSMM,max_components = 2,method='DDRTree') #��ά
HSMM<-orderCells(HSMM)
plot_cell_trajectory(HSMM,color_by = "seurat_clusters")
plot_cell_trajectory(HSMM, color_by = "sample") + facet_wrap(~celltype_human, nrow = 1)