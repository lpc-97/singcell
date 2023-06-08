library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
library(NMF)
options(stringsAsFactors = FALSE)

#1.单独分析
#创建CellChat对象
data.input  <- EC_skin@assays$RNA@data
identity = data.frame(group =EC_skin$Cluster   , row.names = names(EC_skin$Cluster)) #创建细胞的cell label
cellchat <- createCellChat(data.input)
cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels")
groupSize <- as.numeric(table(cellchat@idents))
groupSize #查看每组的细胞数目
CellChatDB <- CellChatDB.human #导入配体受体数据库
#str(CellChatDB) #查看数据库信息
#包含interaction complex cofactor geneInfo四个dataframe
#showDatabaseCategory(CellChatDB) #查看数据库构成

#CellChatDB.use<-subsetDB(CellChatDB,search="Secreted Signaling") #选择其中的Secreted Signaling来进行后续分析
#cellchat@DB<-CellChatDB.use
cellchat@DB<-CellChatDB
#预处理
cellchat<-subsetData(cellchat)
future::plan("multiprocess",workers=4)
cellchat<-identifyOverExpressedGenes(cellchat)
cellchat<-identifyOverExpressedInteractions(cellchat)
cellchat<-projectData(cellchat,PPI.human) #找到配体受体关系后，projectData将配体受体对的表达值投射到PPI上，来对@data.signaling中的表达值进行校正
#推断配体受体水平细胞通讯网络 推断信号通路水平上的通信概率
#根据表达值推测细胞互作的概率
cellchat<-computeCommunProb(cellchat,raw.use=FALSE)
cellchat<-filterCommunication(cellchat)
df.net<-subsetCommunication(cellchat)
write.csv(df.net,"net_lr.csv") #配体受体水平细胞通讯网络

#通过计算链路的数量或汇总通信概率计算细胞之间的聚合通信网络
cellchat<-computeCommunProbPathway(cellchat)
df.netp<-subsetCommunication(cellchat,slot.name='netP')
write.csv(df.netp,'net_pathway.csv')  #信号通路水平的细胞通讯网络

#统计细胞与细胞之间的通信数目
cellchat<-aggregateNet(cellchat)
groupSize<-as.numeric(table(cellchat@idents))
par(mfrow=c(1,1),xpd=TRUE)
netVisual_circle(cellchat@net$count,vertex.weight=groupSize,weight.scale=T,label.edge=F,title.name='Number of interactions')
netVisual_circle(cellchat@net$weight,vertex.weight=groupSize,weight.scale=T,label.edge=F,title.name='Interaction weights/strength')

#保存所有信号通路
cellchat@netP$pathways #查看信号通路
pathways.show.all<-cellchat@netP$pathways
levels(cellchat@idents)
dir.create("all_pathways_com_circle")
setwd("all_pathways_com_circle")
for ( i in 1:length(pathways.show.all)) {
  netVisual(cellchat,signaling=pathways.show.all[i],out.format=c('pdf'),
            vertex.receiver=vertex.reciver)#绘制网络图
  gg<-netAnalysis_contribution(cellchat,signaling=pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"),
         plot=gg, width=5, height=2.5, units='in',dpi=300)
}
setwd('../')

levels(cellchat@idents)
p=netVisual_bubble(cellchat,sources.use=c(3,4),targets.use=c(1,2,5,6),remove.isolate=FALSE)
p #气泡图
ggsave('bubble.pdf',p,width=8,height=12)
#绘制特定通路气泡图
netVisual_bubble(cellchat,sources.use=c(3,4),targets.use=c(1,2,5,6),signaling=c('TGFb'),remove.isolate=FALSE)
plotGeneExpression(cellchat,signaling="TGFb")#计算特定信号通路表达量

pathways.show<-c('TGFb') #查看亚群对于特定通路中的关系、
#计算细胞群的网络中心性指标 识别发送者 接收者 调解者 影响者
cellchat<-netAnalysis_computeCentrality(cellchat,slot.name='netP')
netAnalysis_signalingRole_network(cellchat,signaling=pathways.show,width=8,height=5,font.size=8)
netAnalysis_signalingRole_scatter(cellchat)#在 2D 空间中可视化占主导地位的发送器（源）和接收器（目标）
netAnalysis_signalingRole_scatter(cellchat,signaling = "FN1")#在 2D 空间中可视化占主导地位的发送器（源）和接收器（目标）选择特定的信号通路

#识别细胞的信号流模式
ht1<-netAnalysis_signalingRole_heatmap(cellchat,pattern='outgoing',font.size=5)
ht2<-netAnalysis_signalingRole_heatmap(cellchat,pattern='incoming',font.size=5)
ht1+ht2

#非负矩阵分解
#1 outgoing
selectK(cellchat,pattern='outgoing') #运行慢
nPatterns=3 #写出非负矩阵曲线中第一个出现下降的点
dev.off()
cellchat<-identifyCommunicationPatterns(cellchat,pattern='outgoing',k=nPatterns,width=5,height=8,font.size=6)#热图展示pattern
netAnalysis_river(cellchat,pattern='outgoing')#冲积图展示
netAnalysis_dot(cellchat,pattern='outgoing')#气泡图展示

#2 incoming
selectK(cellchat,pattern='incoming') #运行慢
dev.off()
nPatterns=3 #写出非负矩阵曲线中第一个出现下降的点
cellchat<-identifyCommunicationPatterns(cellchat,pattern='incoming',k=nPatterns,width=5,height=8,font.size=6)#热图展示pattern
netAnalysis_river(cellchat,pattern='incoming')#冲积图展示
netAnalysis_dot(cellchat,pattern='incoming')#气泡图展示


#信号网络的多重和分类学习分析

##信号网络聚类
# 按功能相似性聚类
cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
#可视化
p = netVisual_embedding(cellchat, type = "functional")
p = netVisual_embeddingZoomIn(cellchat, type = "functional")
# 按结构相似性聚类
cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")
#可视化
p = netVisual_embedding(cellchat, type = "structural")
p = netVisual_embeddingZoomIn(cellchat, type = "structural")


#2.对比分析

#配对分析必须保证细胞类型是一样的，才可以进行配对。如果 两个样本的细胞类型不一样又想进行配对分析时，可以用subset把两个样本的细胞类型取成一致的
cco.til <- createCellChat(SMC_CM@assays$RNA@data, meta = SMC_CM@meta.data, group.by = "Cluster")
cco.pbmc <- createCellChat(SMC_Ctrl@assays$RNA@data, meta = SMC_Ctrl@meta.data, group.by = "Cluster")


dir.create("./Compare")
setwd("./Compare")

#对cco.pbmc进行分析
cellchat <- cco.pbmc
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 5)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")
cco.pbmc <- cellchat

#对cco.til进行分析
cellchat <- cco.til
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
#cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
#cellchat <- filterCommunication(cellchat, min.cells = 5)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
#cellchat <- computeNetSimilarity(cellchat, type = "functional")
#cellchat <- netEmbedding(cellchat, type = "functional")
#cellchat <- netClustering(cellchat, type = "functional")
#cellchat <- computeNetSimilarity(cellchat, type = "structural")
#cellchat <- netEmbedding(cellchat, type = "structural")
#cellchat <- netClustering(cellchat, type = "structural")
cco.til <- cellchat

#合并cellchat对象
cco.list<-list(CM=cco.til,Ctrl=cco.pbmc)
cellchat<-mergeCellChat(cco.list,add.names=names(cco.list),cell.prefix=TRUE)

#可视化
#所以细胞群总体观：通讯数目与强度对比
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "count")
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
p <- gg1 + gg2
ggsave("Overview_number_strength.pdf", p, width = 6, height = 4)

#网络图
par(mfrow = c(1,2))
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
#红色是case相对于control上调 蓝色是下调

#数量与强度差异热图
par(mfrow=c(1,1))
h1<-netVisual_heatmap(cellchat)
h2<-netVisual_heatmap(cellchat,measure='weight')
h1+h2

#细胞互作数目对比网络图
par(mfrow = c(1,2))
weight.max <- getMaxWeight(cco.list, attribute = c("idents","count"))
for (i in 1:length(cco.list)) {
  netVisual_circle(cco.list[[i]]@net$count, weight.scale = T, label.edge= F,
                   edge.weight.max = weight.max[2], edge.width.max = 12,
                   title.name = paste0("Number of interactions - ", names(cco.list)[i])) }

#指定细胞互作数量对比网络图
par(mfrow = c(1,2)) s.cell <- c("CD4+ T cells", "CD8+ T cells", "Monocytes")
count1 <- cco.list[[1]]@net$count[s.cell, s.cell]
count2 <- cco.list[[2]]@net$count[s.cell, s.cell]
weight.max <- max(max(count1), max(count2))
netVisual_circle(count1, weight.scale = T, label.edge= T, edge.weight.max=weight.max,edge.width.max=12, title.name = paste0("Number of interactions-", names(cco.list)[1]))
netVisual_circle(count1, weight.scale = T, label.edge= T, edge.weight.max=weight.max,edge.width.max=12, title.name = paste0("Number of interactions-", names(cco.list)[2]))



#保守和特异性信号通路的识别与可视化
## 通路信号强度对比分析
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
p <- gg1 + gg2
ggsave("Compare_pathway_strengh.pdf", p, width = 10, height = 6)

#流行学习识别差异信号通路
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional") #netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5) #netVisual_embeddingPairwiseZoomIn(cellchat, type = "functional", nCol = 2)
cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural") #netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5) #netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2)
p <- rankSimilarity(cellchat, type = "structural") + ggtitle("Structural similarity of pathway")
ggsave("Pathway_Similarity.pdf", p, width = 8, height = 5)

#细胞信号模式对比
library(ComplexHeatmap)
pathway.union <- union(cco.list[[1]]@netP$pathways, cco.list[[2]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(cco.list[[1]], pattern = "all", signaling = pathway.union,title = names(cco.list)[1], width = 8, height = 10)
ht2 = netAnalysis_signalingRole_heatmap(cco.list[[2]], pattern = "all", signaling = pathway.union, title = names(cco.list)[2], width = 8, height = 10)

draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

#输出信号模式对比
pathway.union <- union(cco.list[[1]]@netP$pathways, cco.list[[2]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(cco.list[[1]], pattern = "outgoing", signaling = pathway.union, title = names(cco.list)[1], width = 6, height = 7)
ht2 = netAnalysis_signalingRole_heatmap(cco.list[[2]], pattern = "outgoing", signaling = pathway.union,  title = names(cco.list)[2], width = 6, height = 7)

draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))


#输出信号模式对比
pathway.union <- union(cco.list[[1]]@netP$pathways, cco.list[[2]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(cco.list[[1]], pattern = "incoming", signaling = pathway.union, title = names(cco.list)[1], width = 6, height = 7)
ht2 = netAnalysis_signalingRole_heatmap(cco.list[[2]], pattern = "incoming", signaling = pathway.union,  title = names(cco.list)[2], width = 6, height = 7)

draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

#特定信号通路的对比 只能选择两组中共有的信号
cellchat@netP$CM$pathways#查看CM中的通路
cellchat@netP$Ctrl$pathways#查看CM中的通路
pathways.show <- c("TGFb")
weight.max <- getMaxWeight(cco.list, slot.name = c("netP"), attribute = pathways.show)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(cco.list)) {
  netVisual_aggregate(cco.list[[i]], signaling = pathways.show, layout = "circle",
                      edge.weight.max = weight.max[1], edge.width.max = 10,
                      signaling.name = paste(pathways.show, names(cco.list)[i])) }

#热图
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(cco.list)) {
  ht[[i]] <- netVisual_heatmap(cco.list[[i]], signaling = pathways.show, color.heatmap = "Reds", title.name = paste(pathways.show, "signaling ",names(cco.list)[i]))}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

#和弦图
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(cco.list)) {
  netVisual_aggregate(cco.list[[i]], signaling = pathways.show, layout = "chord", pt.title = 3,title.space=0.05,vertex.label.cex = 0.6, signaling.name = paste(pathways.show, names(cco.list)[i])) }

#配体受体对比分析 气泡图
levels(cellchat@idents$joint)
p <- netVisual_bubble(cellchat, sources.use = c(1,2,3), targets.use = c(1,2,3), comparison =c(1,2),angle.x=45)
ggsave("Compare_LR_bubble.pdf", p, width = 12, height = 8)
netVisual_bubble(cellchat,sources.use = c(1,2,3), targets.use = c(1,2,3), comparison =c(1,2),signaling=c('FN1'),remove.isolate=FALSE) #特定信号通路气泡图展示 但只能用于共同通路

#气泡图展示上调或下调的配体受体对
p1 <- netVisual_bubble(cellchat, sources.use =  c(1,2,3), targets.use =  c(1,2,3), comparison =c(1,2), max.dataset = 2, title.name = "Increased signaling in Ctrl", angle.x = 45,remove.isolate=T)
p2 <- netVisual_bubble(cellchat,  sources.use =  c(1,2,3), targets.use =  c(1,2,3), comparison =c(1,2), max.dataset = 1, title.name = "Increased signaling in CM", angle.x = 45,remove.isolate=T)
p1+p2

#和弦图
par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(cco.list)) {
  netVisual_chord_gene(cco.list[[i]], sources.use =  c(1,2,3), targets.use =  c(1,2,3), signaling="TGFb", lab.cex = 0.6, legend.pos.x = 10, legend.pos.y = 20, title.name = paste0("Signaling from Treg - ", names(cco.list)[i])) }


netAnalysis_signalingChanges_scatter(cellchat,idents.use='Fibroblast_like')
