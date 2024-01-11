#####cellChat细胞通讯############
dir.create('results/05.cellchat')
sce_anno=data.frame(cell=rownames(tumor_sc@meta.data),
                                          AS_group=tumor_sc@meta.data$AS_group)

library(CellChat)
sce1=sce
colnames(sce1@meta.data)
table(sce1$cell_type2)
table(sce1$cell_type)
sce1$AS_group=NA
sce1$AS_group[sce_anno[,'cell']] = ifelse(tumor_sc$AS[sce_anno[,'cell']]<0,'Low','High')
table(sce1$AS_group)

sce1$group=ifelse(sce1$cell_type2=='cancer cell',paste0(sce1$AS_group,' neoplastic'),sce1$cell_type)
table(sce1$group)

#教程链接 dooooob https://www.bilibili.com/read/cv12949609 
# #创建cellchat对象
# DefaultAssay(sce1) <- 'RNA'
# cellchat <- createCellChat(object = sce1, meta = sce1@meta.data, group.by = "group")
# cellchat@DB <- subsetDB(CellChatDB = CellChatDB.human, search = "Secreted Signaling")
# cellchat <- subsetData(cellchat)
# # cellchat@data.signaling
# 
# cellchat <- identifyOverExpressedGenes(cellchat)
# cellchat <- identifyOverExpressedInteractions(cellchat)
# cellchat <- projectData(cellchat, PPI.human)
# #分析胞间通讯和互作
# 
# #首先计算得到每一个互作的communication probability（具体原理参见发表原文），并通过permutation test得到对应互作的p值。这些结果存在于cellchat@net
# cellchat <- computeCommunProb(cellchat)
# cellchat <- computeCommunProbPathway(cellchat)
# cellchat <- aggregateNet(cellchat)
# cellchat <- netAnalysis_computeCentrality(cellchat)
#之后可以关注communication pattern（outgoing和incoming），得到cellchat@netP$pattern
#先寻找合适的k，即pattern数目，以outgoing为例
# library(NMF)
# selectK(cellchat, pattern = "outgoing")
# selectK(cellchat, pattern = "incoming")
# save(cellchat,file='results/06.cellchat/cellchat.RData')
load('results/06.cellchat/cellchat.RData')
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = 3)
#这里根据上面的结果设置pattern为3，因为两种score都在3的时候下降很多 

#####细胞通讯可视化######
groupSize <- as.numeric(table(cellchat@idents))

# par(mfrow=c(1,2))
# netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, 
#                  title.name = "Number of interactions")
# netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, 
#                  title.name = "Interaction weights/strength") 
# dev.off()
# 
# 看各个细胞类群的Signaling pathway情况
# 
# 1. 通过heatmap看outgoing和incoming的communication probability
# netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
# netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
# dev.off()


#查看特定通路在群体间的联系，单个信号通路细胞互作可视化
cellchat@netP$pathways  #16
pdf('results/06.cellchat/cellchat.pdf',height = 12,width = 16)
par(mfrow=c(3,5))
for (i in 1:length(cellchat@netP$pathways )) {
  netVisual_aggregate(cellchat, signaling = cellchat@netP$pathways[i], layout = "circle") 
}
dev.off()

#计算每个细胞群的网络中心指标，识别没累细胞在信号通路中的角色
for (i in 1:length(cellchat@netP$pathways)) {
  pdf(paste0('results/06.cellchat/',cellchat@netP$pathways[i],'_signaling.pdf'),height = 5,width = 7)
  netAnalysis_signalingRole_network(cellchat, signaling = cellchat@netP$pathways[i], 
                                    width = 8, height = 2.5, font.size = 10)
  dev.off()
}
#自动可视化该信号通路内富集的重要基因
plotGeneExpression(cellchat, signaling = cellchat@netP$pathways[1])



P=5
netVisual_aggregate(cellchat, signaling = cellchat@netP$pathways[P], layout = "circle",) 
pdf('results/06.cellchat/SPP1_pairLR_net_2.pdf',height = 5,width = 6)
netAnalysis_signalingRole_network(cellchat, signaling = cellchat@netP$pathways[P], 
                                  width = 8, height = 2.5, font.size = 10)
dev.off()
netAnalysis_contribution(cellchat, signaling = cellchat@netP$pathways[P]) 
pairLR <- extractEnrichedLR(cellchat, signaling = cellchat@netP$pathways[P] , geneLR.return = FALSE)
pairLR
LR.show <- pairLR[2,]
pdf('results/06.cellchat/SPP1_pairLR_net_1.pdf',height = 6,width = 6)
netVisual_individual(cellchat, signaling = cellchat@netP$pathways[P], pairLR.use = LR.show, layout = "circle")
dev.off()



P=1
netVisual_aggregate(cellchat, signaling = cellchat@netP$pathways[P], layout = "circle",) 
pdf('results/06.cellchat/MK_pairLR_net_2.pdf',height = 5,width = 6)
netAnalysis_signalingRole_network(cellchat, signaling = cellchat@netP$pathways[P], 
                                  width = 8, height = 2.5, font.size = 10)
dev.off()
netAnalysis_contribution(cellchat, signaling = cellchat@netP$pathways[P]) 
pairLR <- extractEnrichedLR(cellchat, signaling = cellchat@netP$pathways[P] , geneLR.return = FALSE)
pairLR
LR.show <- pairLR[5,]
pdf('results/06.cellchat/MK_pairLR_net_1.pdf',height = 6,width = 6)
netVisual_individual(cellchat, signaling = cellchat@netP$pathways[P], pairLR.use = LR.show, layout = "circle")
dev.off()

P=10
netVisual_aggregate(cellchat, signaling = cellchat@netP$pathways[P], layout = "circle") 
pdf('results/06.cellchat/GDF_pairLR_net_2.pdf',height = 5,width = 6)
netAnalysis_signalingRole_network(cellchat, signaling = cellchat@netP$pathways[P], 
                                  width = 8, height = 2.5, font.size = 10)
dev.off()
netAnalysis_contribution(cellchat, signaling = cellchat@netP$pathways[P]) 
pairLR <- extractEnrichedLR(cellchat, signaling = cellchat@netP$pathways[P] , geneLR.return = FALSE)
pairLR
LR.show <- pairLR[1,]
pdf('results/06.cellchat/GDF_pairLR_net_1.pdf',height = 6,width = 12)
#par(mfrow=c(1,2))
netVisual_individual(cellchat, signaling = cellchat@netP$pathways[P], 
                     pairLR.use = LR.show, layout = "circle")
dev.off()


