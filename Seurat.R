library(Seurat)
library(cowplot)
library(patchwork)
library(dplyr)
data1<-Read10X(data_dir1 <- "~/data/" ) #��ȡ����
data1<-CreateSeuratObject(counts =data1,project = "data1",min.cells = 3,min.features = 200) #����ת�������ҹ���
data1[["percent.mt"]]<-PercentageFeatureSet(data1,pattern = "MT-") #�鿴������
HB.genes_total<-c("HBA1","HBA2","HBB","HBD","HBE1","HBG2","HBM","HBQ1","HB2")
HB_m<-match(HB.genes_total,rownames(data1@assays$RNA))
HB.genes<-rownames(data1@assays$RNA)[HB_m]
HB.genes<-HB.genes[!is.na(HB.genes)]
data1[["percent.HB"]]<-PercentageFeatureSet(data1,features = HB.genes) #�����ϸ���ٷֱ�
VlnPlot(data1,features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.HB")) #С����ͼ�鿴����
data1 <- subset(data1, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
data1<-NormalizeData(data1)
data1 <- FindVariableFeatures(data1, selection.method = "vst", nfeatures = 2000) #�ɱ����
#data<-FindIntegrationAnchors(object.list = list(data1, data2), dims = 1:20) #��������ϲ�
#data<-IntegrateData(anchorset=data,dims=1:20)
#DefaultAssay(data)<-"integrated" #�ϲ���������
#scale����
all.genes<-rownames(data)
data<-ScaleData(data,features = all.genes,vars.to.regress = "percent.mt")

data<-RunPCA(data,npcs=30,verbose=FALSE) #PCA
data <- FindNeighbors(data, dims = 1:10)
data <- FindClusters(data, resolution = 0.5) #��ά����
data<-RunUMAP(data,dims = 1:10)
data.markers <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

new.cluster.ids<-c("Naive CD4 T","Memory CD4 T","CD14+Mono","B","CD8 T","FCGR3A+Mono","NK","DC","Platelet")
names(new.cluster.ids)<-levels(data)
data <- RenameIdents(data, new.cluster.ids)
DimPlot(data, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

#ϸ�����ڷ���
?cc.genes
length(c(cc.genes$s.genes,cc.genes$g2m.genes))
head(c(cc.genes$s.genes,cc.genes$g2m.genes))
CaseMatch(c(cc.genes$s.genes,cc.genes$g2m.genes),VariableFeatures(data))
g2m_genes=cc.genes$g2m.genes
g2m_genes=CaseMatch(search = g2m_genes,match = rownames(data))
s_genes=cc.genes$s.genes
s_genes=CaseMatch(search = s_genes,match = rownames(data))
data<-CellCycleScoring(object =data,g2m.features = g2m_genes,s.features=s_genes)
head(data)
DimPlot(data,reduction = "umap",group.by = "Phase")

singleRע��ϸ����Ⱥ
library(Seurat) 
library(SingleR)
library(ggplot2)
library(reshape2)
library(celldex)
##�������ݼ�
hpca.se=HumanPrimaryCellAtlasData() 
Blue.se=BlueprintEncodeData() 
Immune.se=DatabaseImmuneCellExpressionData()
Nover.se=NovershternHematopoieticData()
MonacoIm.se=MonacoImmuneData()
ImmGen.se=ImmGenData() #(��)
Mouse.se=MouseRNAseqData() #(��)
meta=data@meta.data 
data_for_SingleR <- GetAssayData(data, slot="data") ##��ȡ��׼������
data.hesc <- SingleR(test = data_for_SingleR, ref = hpca.se, labels = hpca.se$label.main) # ʹ��HumanPrimaryCellAtlasData�ο����ݼ�
table(data.hesc[[i]]$labels,meta$seurat_clusters) 
plotScoreHeatmap(data.hesc)
data@meta.data$labels <-data.hesc$labels
DimPlot(data, group.by = c("seurat_clusters", "labels"),reduction = "umap")

#ϸ����Ⱥ����Է���
av <-AverageExpression(pbmc,group.by = "Cluster",assays = "RNA")
av=av[[1]]
cg=names(tail(sort(apply(av, 1, sd)),1000))#ɸѡ��1000��������з���
View(av[cg,])
View(cor(av[cg,],method = 'spearman'))
pheatmap::pheatmap(cor(av[cg,],method = 'spearman'))#������ͼ