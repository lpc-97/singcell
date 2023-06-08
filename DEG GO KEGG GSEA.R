library(Seurat)
library(clusterProfiler)
library(tidyverse)
library(patchwork)
library(monocle)
library(org.Hs.eg.db) #其他物种看其他的就可以（这个是人的）

dge.cluster <- FindMarkers(data,ident.1 = "group1",ident.2 = "group2",group.by="Group") #这里可以选择分组
sig_dge.celltype <- subset(dge.cluster, p_val_adj<0.01&abs(avg_log2FC)>0.25)
write.csv(sig_dge.celltype,file="DEG.csv") #R语言abs函数的意思是绝对值

#GO分析
ego_ALL <- enrichGO(gene          = row.names(sig_dge.celltype),OrgDb         = 'org.Hs.eg.db',keyType       = 'SYMBOL',ont           = "ALL",pAdjustMethod = "BH",pvalueCutoff  = 0.01,qvalueCutoff  = 0.05)
ego_ALL <- data.frame(ego_ALL)
write.csv(ego_ALL,'enrichGO.csv')

#KEGG分析
genelist = bitr(row.names(sig_dge.celltype), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
genelist <- pull(genelist,ENTREZID)
ekegg <- enrichKEGG(gene = genelist, organism = 'hsa',pvalueCutoff = 0.01,qvalueCutoff = 0.01,minGSSize = 1)
write.csv(as.data.frame(ekegg@result), file="kegg.csv") #KEGG结果保存
genelist = bitr(genelist, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db") #可以利用这个对以上基因进行转化
write.csv(genelist,'KEGG_genelist.csv')

#如果是CSV文件
#GO分析
filepath<-file.choose()
df2<-read.csv(filepath,header = TRUE,stringsAsFactors = F)
ego_ALL <- enrichGO(gene          =df2$genelist,OrgDb         = 'org.Hs.eg.db',keyType       = 'SYMBOL',ont           = "ALL",pAdjustMethod = "BH",pvalueCutoff  = 0.01,qvalueCutoff  = 0.05)
ego_ALL <- data.frame(ego_ALL)
write.csv(ego_ALL,'enrichGO.csv')

#KEGG分析
genelist = bitr(df2$genelist, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
genelist <- pull(genelist,ENTREZID)
ekegg <- enrichKEGG(gene = genelist, organism = 'hsa',pvalueCutoff = 0.01,qvalueCutoff = 0.01,minGSSize = 1)
write.csv(as.data.frame(ekegg@result), file="kegg.csv") #KEGG结果保存
genelist = bitr(genelist, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db") #可以利用这个对以上基因进行转化
write.csv(genelist,'KEGG_genelist.csv')

#GSEA分析
library(Seurat)
library(clusterProfiler)
library(tidyverse)
library(patchwork)
library(org.Hs.eg.db) #其他物种看其他的就可以（这个是人的）
library(msigdbr)
library(dplyr)
library(tibble)
library(fgsea)
library(ggplot2)
library(enrichplot)
geneset <- read.gmt("c2.cp.kegg.v7.5.1.entrez.gmt") #选择想要使用的数据集

markers <- FindMarkers(combinedcc6n3,ident.1 = "Artery",ident.2 = c("Tip",'Vein','AVM','BBB','EndMT'),min.pct = 0.25, logfc.threshold = 0) #寻找想要的marker

gs <-bitr(rownames(markers), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db") #进行gene ID转化
markers1<-cbind(markers[gs[,1],],gs) 
geneList = markers1$avg_log2FC #这里的是基因名和avg_log2FC两列值
names(geneList) = markers1$ENTREZID
geneList = sort(geneList,decreasing = T) #进行排序
egmt <- GSEA(geneList, TERM2GENE=geneset,verbose=F)
egmt1<- setReadable(egmt,OrgDb=org.Hs.eg.db, keyType = "ENTREZID")
y=data.frame(egmt1)
write.csv(y,file = "BBB.csv") #保存富集的通路

for(i in seq_along(egmt@result$ID)){
  p <- gseaplot2(egmt, geneSetID = i, title = egmt@result$ID[i])
  filename <- paste0('GSEA', egmt@result$ID[i], '.png')
  ggsave(filename = filename, p, width = 4, height = 4)} #作图并保存前面所有分析出

gseaplot2(egmt,"KEGG_RIBOSOME",title = "KEGG_RIBOSOME") #对特定的通路进行分析