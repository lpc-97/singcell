library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
library(NMF)
options(stringsAsFactors = FALSE)

#1.��������
#����CellChat����
data.input  <- EC_skin@assays$RNA@data
identity = data.frame(group =EC_skin$Cluster   , row.names = names(EC_skin$Cluster)) #����ϸ����cell label
cellchat <- createCellChat(data.input)
cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels")
groupSize <- as.numeric(table(cellchat@idents))
groupSize #�鿴ÿ���ϸ����Ŀ
CellChatDB <- CellChatDB.human #���������������ݿ�
#str(CellChatDB) #�鿴���ݿ���Ϣ
#����interaction complex cofactor geneInfo�ĸ�dataframe
#showDatabaseCategory(CellChatDB) #�鿴���ݿ⹹��

#CellChatDB.use<-subsetDB(CellChatDB,search="Secreted Signaling") #ѡ�����е�Secreted Signaling�����к�������
#cellchat@DB<-CellChatDB.use
cellchat@DB<-CellChatDB
#Ԥ����
cellchat<-subsetData(cellchat)
future::plan("multiprocess",workers=4)
cellchat<-identifyOverExpressedGenes(cellchat)
cellchat<-identifyOverExpressedInteractions(cellchat)
cellchat<-projectData(cellchat,PPI.human) #�ҵ����������ϵ��projectData����������Եı���ֵͶ�䵽PPI�ϣ�����@data.signaling�еı���ֵ����У��
#�ƶ���������ˮƽϸ��ͨѶ���� �ƶ��ź�ͨ·ˮƽ�ϵ�ͨ�Ÿ���
#���ݱ���ֵ�Ʋ�ϸ�������ĸ���
cellchat<-computeCommunProb(cellchat,raw.use=FALSE)
cellchat<-filterCommunication(cellchat)
df.net<-subsetCommunication(cellchat)
write.csv(df.net,"net_lr.csv") #��������ˮƽϸ��ͨѶ����

#ͨ��������·�����������ͨ�Ÿ��ʼ���ϸ��֮��ľۺ�ͨ������
cellchat<-computeCommunProbPathway(cellchat)
df.netp<-subsetCommunication(cellchat,slot.name='netP')
write.csv(df.netp,'net_pathway.csv')  #�ź�ͨ·ˮƽ��ϸ��ͨѶ����

#ͳ��ϸ����ϸ��֮���ͨ����Ŀ
cellchat<-aggregateNet(cellchat)
groupSize<-as.numeric(table(cellchat@idents))
par(mfrow=c(1,1),xpd=TRUE)
netVisual_circle(cellchat@net$count,vertex.weight=groupSize,weight.scale=T,label.edge=F,title.name='Number of interactions')
netVisual_circle(cellchat@net$weight,vertex.weight=groupSize,weight.scale=T,label.edge=F,title.name='Interaction weights/strength')

#���������ź�ͨ·
cellchat@netP$pathways #�鿴�ź�ͨ·
pathways.show.all<-cellchat@netP$pathways
levels(cellchat@idents)
dir.create("all_pathways_com_circle")
setwd("all_pathways_com_circle")
for ( i in 1:length(pathways.show.all)) {
  netVisual(cellchat,signaling=pathways.show.all[i],out.format=c('pdf'),
            vertex.receiver=vertex.reciver)#��������ͼ
  gg<-netAnalysis_contribution(cellchat,signaling=pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"),
         plot=gg, width=5, height=2.5, units='in',dpi=300)
}
setwd('../')

levels(cellchat@idents)
p=netVisual_bubble(cellchat,sources.use=c(3,4),targets.use=c(1,2,5,6),remove.isolate=FALSE)
p #����ͼ
ggsave('bubble.pdf',p,width=8,height=12)
#�����ض�ͨ·����ͼ
netVisual_bubble(cellchat,sources.use=c(3,4),targets.use=c(1,2,5,6),signaling=c('TGFb'),remove.isolate=FALSE)
plotGeneExpression(cellchat,signaling="TGFb")#�����ض��ź�ͨ·������

pathways.show<-c('TGFb') #�鿴��Ⱥ�����ض�ͨ·�еĹ�ϵ��
#����ϸ��Ⱥ������������ָ�� ʶ������ ������ ������ Ӱ����
cellchat<-netAnalysis_computeCentrality(cellchat,slot.name='netP')
netAnalysis_signalingRole_network(cellchat,signaling=pathways.show,width=8,height=5,font.size=8)
netAnalysis_signalingRole_scatter(cellchat)#�� 2D �ռ��п��ӻ�ռ������λ�ķ�������Դ���ͽ�������Ŀ�꣩
netAnalysis_signalingRole_scatter(cellchat,signaling = "FN1")#�� 2D �ռ��п��ӻ�ռ������λ�ķ�������Դ���ͽ�������Ŀ�꣩ѡ���ض����ź�ͨ·

#ʶ��ϸ�����ź���ģʽ
ht1<-netAnalysis_signalingRole_heatmap(cellchat,pattern='outgoing',font.size=5)
ht2<-netAnalysis_signalingRole_heatmap(cellchat,pattern='incoming',font.size=5)
ht1+ht2

#�Ǹ�����ֽ�
#1 outgoing
selectK(cellchat,pattern='outgoing') #������
nPatterns=3 #д���Ǹ����������е�һ�������½��ĵ�
dev.off()
cellchat<-identifyCommunicationPatterns(cellchat,pattern='outgoing',k=nPatterns,width=5,height=8,font.size=6)#��ͼչʾpattern
netAnalysis_river(cellchat,pattern='outgoing')#���ͼչʾ
netAnalysis_dot(cellchat,pattern='outgoing')#����ͼչʾ

#2 incoming
selectK(cellchat,pattern='incoming') #������
dev.off()
nPatterns=3 #д���Ǹ����������е�һ�������½��ĵ�
cellchat<-identifyCommunicationPatterns(cellchat,pattern='incoming',k=nPatterns,width=5,height=8,font.size=6)#��ͼչʾpattern
netAnalysis_river(cellchat,pattern='incoming')#���ͼչʾ
netAnalysis_dot(cellchat,pattern='incoming')#����ͼչʾ


#�ź�����Ķ��غͷ���ѧϰ����

##�ź��������
# �����������Ծ���
cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
#���ӻ�
p = netVisual_embedding(cellchat, type = "functional")
p = netVisual_embeddingZoomIn(cellchat, type = "functional")
# ���ṹ�����Ծ���
cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")
#���ӻ�
p = netVisual_embedding(cellchat, type = "structural")
p = netVisual_embeddingZoomIn(cellchat, type = "structural")


#2.�Աȷ���

#��Է������뱣֤ϸ��������һ���ģ��ſ��Խ�����ԡ���� ����������ϸ�����Ͳ�һ�����������Է���ʱ��������subset������������ϸ������ȡ��һ�µ�
cco.til <- createCellChat(SMC_CM@assays$RNA@data, meta = SMC_CM@meta.data, group.by = "Cluster")
cco.pbmc <- createCellChat(SMC_Ctrl@assays$RNA@data, meta = SMC_Ctrl@meta.data, group.by = "Cluster")


dir.create("./Compare")
setwd("./Compare")

#��cco.pbmc���з���
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

#��cco.til���з���
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

#�ϲ�cellchat����
cco.list<-list(CM=cco.til,Ctrl=cco.pbmc)
cellchat<-mergeCellChat(cco.list,add.names=names(cco.list),cell.prefix=TRUE)

#���ӻ�
#����ϸ��Ⱥ����ۣ�ͨѶ��Ŀ��ǿ�ȶԱ�
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "count")
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
p <- gg1 + gg2
ggsave("Overview_number_strength.pdf", p, width = 6, height = 4)

#����ͼ
par(mfrow = c(1,2))
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
#��ɫ��case�����control�ϵ� ��ɫ���µ�

#������ǿ�Ȳ�����ͼ
par(mfrow=c(1,1))
h1<-netVisual_heatmap(cellchat)
h2<-netVisual_heatmap(cellchat,measure='weight')
h1+h2

#ϸ��������Ŀ�Ա�����ͼ
par(mfrow = c(1,2))
weight.max <- getMaxWeight(cco.list, attribute = c("idents","count"))
for (i in 1:length(cco.list)) {
  netVisual_circle(cco.list[[i]]@net$count, weight.scale = T, label.edge= F,
                   edge.weight.max = weight.max[2], edge.width.max = 12,
                   title.name = paste0("Number of interactions - ", names(cco.list)[i])) }

#ָ��ϸ�����������Ա�����ͼ
par(mfrow = c(1,2)) s.cell <- c("CD4+ T cells", "CD8+ T cells", "Monocytes")
count1 <- cco.list[[1]]@net$count[s.cell, s.cell]
count2 <- cco.list[[2]]@net$count[s.cell, s.cell]
weight.max <- max(max(count1), max(count2))
netVisual_circle(count1, weight.scale = T, label.edge= T, edge.weight.max=weight.max,edge.width.max=12, title.name = paste0("Number of interactions-", names(cco.list)[1]))
netVisual_circle(count1, weight.scale = T, label.edge= T, edge.weight.max=weight.max,edge.width.max=12, title.name = paste0("Number of interactions-", names(cco.list)[2]))



#���غ��������ź�ͨ·��ʶ������ӻ�
## ͨ·�ź�ǿ�ȶԱȷ���
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
p <- gg1 + gg2
ggsave("Compare_pathway_strengh.pdf", p, width = 10, height = 6)

#����ѧϰʶ������ź�ͨ·
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional") #netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5) #netVisual_embeddingPairwiseZoomIn(cellchat, type = "functional", nCol = 2)
cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural") #netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5) #netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2)
p <- rankSimilarity(cellchat, type = "structural") + ggtitle("Structural similarity of pathway")
ggsave("Pathway_Similarity.pdf", p, width = 8, height = 5)

#ϸ���ź�ģʽ�Ա�
library(ComplexHeatmap)
pathway.union <- union(cco.list[[1]]@netP$pathways, cco.list[[2]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(cco.list[[1]], pattern = "all", signaling = pathway.union,title = names(cco.list)[1], width = 8, height = 10)
ht2 = netAnalysis_signalingRole_heatmap(cco.list[[2]], pattern = "all", signaling = pathway.union, title = names(cco.list)[2], width = 8, height = 10)

draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

#����ź�ģʽ�Ա�
pathway.union <- union(cco.list[[1]]@netP$pathways, cco.list[[2]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(cco.list[[1]], pattern = "outgoing", signaling = pathway.union, title = names(cco.list)[1], width = 6, height = 7)
ht2 = netAnalysis_signalingRole_heatmap(cco.list[[2]], pattern = "outgoing", signaling = pathway.union,  title = names(cco.list)[2], width = 6, height = 7)

draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))


#����ź�ģʽ�Ա�
pathway.union <- union(cco.list[[1]]@netP$pathways, cco.list[[2]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(cco.list[[1]], pattern = "incoming", signaling = pathway.union, title = names(cco.list)[1], width = 6, height = 7)
ht2 = netAnalysis_signalingRole_heatmap(cco.list[[2]], pattern = "incoming", signaling = pathway.union,  title = names(cco.list)[2], width = 6, height = 7)

draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

#�ض��ź�ͨ·�ĶԱ� ֻ��ѡ�������й��е��ź�
cellchat@netP$CM$pathways#�鿴CM�е�ͨ·
cellchat@netP$Ctrl$pathways#�鿴CM�е�ͨ·
pathways.show <- c("TGFb")
weight.max <- getMaxWeight(cco.list, slot.name = c("netP"), attribute = pathways.show)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(cco.list)) {
  netVisual_aggregate(cco.list[[i]], signaling = pathways.show, layout = "circle",
                      edge.weight.max = weight.max[1], edge.width.max = 10,
                      signaling.name = paste(pathways.show, names(cco.list)[i])) }

#��ͼ
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(cco.list)) {
  ht[[i]] <- netVisual_heatmap(cco.list[[i]], signaling = pathways.show, color.heatmap = "Reds", title.name = paste(pathways.show, "signaling ",names(cco.list)[i]))}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

#����ͼ
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(cco.list)) {
  netVisual_aggregate(cco.list[[i]], signaling = pathways.show, layout = "chord", pt.title = 3,title.space=0.05,vertex.label.cex = 0.6, signaling.name = paste(pathways.show, names(cco.list)[i])) }

#��������Աȷ��� ����ͼ
levels(cellchat@idents$joint)
p <- netVisual_bubble(cellchat, sources.use = c(1,2,3), targets.use = c(1,2,3), comparison =c(1,2),angle.x=45)
ggsave("Compare_LR_bubble.pdf", p, width = 12, height = 8)
netVisual_bubble(cellchat,sources.use = c(1,2,3), targets.use = c(1,2,3), comparison =c(1,2),signaling=c('FN1'),remove.isolate=FALSE) #�ض��ź�ͨ·����ͼչʾ ��ֻ�����ڹ�ͬͨ·

#����ͼչʾ�ϵ����µ������������
p1 <- netVisual_bubble(cellchat, sources.use =  c(1,2,3), targets.use =  c(1,2,3), comparison =c(1,2), max.dataset = 2, title.name = "Increased signaling in Ctrl", angle.x = 45,remove.isolate=T)
p2 <- netVisual_bubble(cellchat,  sources.use =  c(1,2,3), targets.use =  c(1,2,3), comparison =c(1,2), max.dataset = 1, title.name = "Increased signaling in CM", angle.x = 45,remove.isolate=T)
p1+p2

#����ͼ
par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(cco.list)) {
  netVisual_chord_gene(cco.list[[i]], sources.use =  c(1,2,3), targets.use =  c(1,2,3), signaling="TGFb", lab.cex = 0.6, legend.pos.x = 10, legend.pos.y = 20, title.name = paste0("Signaling from Treg - ", names(cco.list)[i])) }


netAnalysis_signalingChanges_scatter(cellchat,idents.use='Fibroblast_like')