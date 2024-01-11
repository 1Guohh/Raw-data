############4.1 单细胞基础流程+细胞注释########
library(Seurat)
library(dplyr)
library(ggplot2)
library(magrittr)
library(gtools)
library(stringr)
library(Matrix)
library(tidyverse)
library(patchwork)
library(data.table)
library(RColorBrewer)
dir.create('results/04.scRNA')
dir_name=list.files('origin_datas/GEO/GSE203360_RAW/',pattern = '.txt')
datalist=list()
for (i in 1:length(dir_name)){
  #i=1
  files = paste0("origin_datas/GEO/GSE203360_RAW/",dir_name[i])
  counts=fread(file = files,data.table = T,sep = '\t',check.names = F)
  counts=data.frame(counts)
  rownames(counts)=counts[,1]
  counts=counts[,-1]
  rownames(counts) <- gsub("_","-", rownames(counts))
  Samples1=stringr::str_split_fixed(dir_name[i],'_',4)[,1]
  Patient=stringr::str_split_fixed(dir_name[i],'_',4)[,2]
  Type=stringr::str_split_fixed(dir_name[i],'_',4)[,3]
  colnames(counts)=paste0(Samples1,'_',colnames(counts))
  datalist[[i]]<- CreateSeuratObject(counts=counts,project = Samples1,min.cells = 3, min.features = 200) 
  datalist[[i]] <- AddMetaData(datalist[[i]] , Samples1,col.name = "Samples")
  datalist[[i]] <- AddMetaData(datalist[[i]] , Patient,col.name = "Patient")
  datalist[[i]] <- AddMetaData(datalist[[i]] , Type,col.name = "Type")
}
names(datalist)=stringr::str_split_fixed(dir_name,'_',4)[,1]


for (i in 1:length(datalist)){
  sce <- datalist[[i]]
  sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")# 计算线粒体占比
  sce[["percent.Ribo"]] <- PercentageFeatureSet(sce, pattern = "^RP[SL]")# 计算rRNA占比
  datalist[[i]] <- sce
  rm(sce)
}
sce <- merge(datalist[[1]],y=datalist[2:length(datalist)])


#细胞数统计
raw_meta=sce@meta.data
raw_count <- table(raw_meta$Samples)
raw_count
sum(raw_count)#19383
pearplot_befor<-VlnPlot(sce,group.by ='Samples', 
                        features = c("nFeature_RNA", "nCount_RNA","percent.mt"),
                        pt.size = 0, 
                        ncol = 3)
pearplot_befor


#过滤
datalist <- lapply(X = datalist, FUN = function(x) {
  x<-subset(x,subset =  nFeature_RNA < 4500 & 
              nFeature_RNA > 200 &
              percent.mt<25)
})


#合并数据
sce <- merge(datalist[[1]],y=datalist[2:length(datalist)])
clean_meta=sce@meta.data
clean_count <- table(clean_meta$Samples)
clean_count
sum(clean_count)#18134
pearplot_after <- VlnPlot(sce,group.by ='Samples', 
                          features = c("nFeature_RNA", "nCount_RNA","percent.mt"),
                          pt.size = 0, 
                          ncol = 3)
pearplot_after

summary_cells <- as.data.frame(cbind(raw_count,clean_count))
summary_cells
writeMatrix(summary_cells,'results/04.scRNA/summary_cells.txt')

#过滤后细胞数的统计
head(summary_cells)
summary_cells$Samples=rownames(summary_cells)
summary_cells1=summary_cells
summary_cells1$fit=summary_cells1$raw_count-summary_cells$clean_count
summary_cells1=reshape2::melt(summary_cells1[,c("clean_count","Samples","fit")])
summary_cells1$variable=ifelse(summary_cells1$variable=='fit','Filtering number','Number after filtering')
summary_cells1$variable=factor(summary_cells1$variable,levels = c('Filtering number','Number after filtering'))
subset_barplot<-ggbarplot(summary_cells1, x = "Samples", y="value", color="black", fill="variable",
                          legend="right", 
                          legend.title="", main="Before and after filtering",
                          font.main = c(14,"bold", "black"), font.x = c(12, "bold"), 
                          font.y=c(12,"bold")) + 
  theme_bw() +
  rotate_x_text() + scale_fill_manual(values=c("#66C2A5","#FC8D62" ))+
  labs(x = "", y = "Cell Number") + 
  theme(axis.text.x=element_text(angle = 0),
        axis.ticks.x=element_blank(),
        axis.title = element_text(face = "bold"), 
        plot.title = element_text(face = "bold"), 
        legend.title = element_text(face = "bold")) +coord_flip()
subset_barplot

save(datalist,file = 'results/04.scRNA/datalist.RData')

#降维
load('results/04.scRNA/datalist.RData')
sce <- merge(datalist[[1]],y=datalist[2:length(datalist)])
#数据标准化
sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
#使用均值与方差之间的关系，来挑选高变基因，便于进行特征选择筛选生物学相关基因，默认返回前2000个高变基因进入下游分析
sce <- FindVariableFeatures(sce,
                            selection.method = "vst",
                            nfeatures = 2000,
                            mean.cutoff=c(0.0125,3),
                            dispersion.cutoff =c(1.5,Inf))
#使用[ScaleData()进行数据归一化。归一化的功能：使得每一个基因在所有cell中的表达均值为0，方差为1
sce <- ScaleData(sce, features =  rownames(sce))
#对归一化后的数据进行PCA分析
sce <- RunPCA(sce, features = VariableFeatures(sce))

# dimplot1 <- DimPlot(sce, reduction = "pca",group.by = 'Samples')
# elbowplot1 <- ElbowPlot(sce, ndims=50, reduction="pca")
# sc_pca <- dimplot1+elbowplot1
# pdf('results/01.scRNA/PCA.pdf',height = 8,width = 15)
# sc_pca
# dev.off()

Dims <- 30
sce <- RunTSNE(sce, dims=1:Dims, reduction="pca")

before_batch=DimPlot(sce,group.by='Samples',
                     reduction="tsne",
                     label = "T",
                     pt.size = 0.2,
                     label.size = 0)+
  ggtitle('Beforer batch')


#去批次
# Normalizing the data
load('results/04.scRNA/datalist.RData')
for (i in 1:length(datalist)){
  datalist[[i]]<-NormalizeData(datalist[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
  datalist[[i]]<-FindVariableFeatures(datalist[[i]],
                                      selection.method = "vst",
                                      nfeatures = 2000,
                                      mean.cutoff=c(0.0125,3),
                                      dispersion.cutoff =c(1.5,Inf))
}

#使用CCA的方法进行剔除批次效应
datalist <- FindIntegrationAnchors(object.list = datalist, dims = 1:40,
                                   reduction = c("cca", "rpca")[1])
sce <- IntegrateData(anchorset = datalist, dims = 1:40)
#ScaleData
DefaultAssay(sce) <- "integrated"

#sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
sce <- FindVariableFeatures(sce,
                            selection.method = "vst",
                            nfeatures = 2000,
                            mean.cutoff=c(0.0125,3),
                            dispersion.cutoff =c(1.5,Inf))
sce <- ScaleData(sce, features =  rownames(sce))
sce <- RunPCA(sce, features = VariableFeatures(sce))

dimplot <- DimPlot(sce, reduction = "pca",group.by = 'Samples')
elbowplot <- ElbowPlot(sce, ndims=50, reduction="pca")
sc_pca <- dimplot+elbowplot
pdf('results/04.scRNA//FigS2.pdf',height = 6,width = 15)
sc_pca
dev.off()

pdf('results/04.scRNA//FigS2_PCA.pdf',height = 6,width = 6)
elbowplot
dev.off()


###tsne 降维
Dims <- 40
sce <- RunTSNE(sce,
               dims=1:Dims,
               reduction="pca")
after_batch=DimPlot(sce,group.by='Samples',
                    reduction="tsne",
                    label = "T",
                    pt.size = 0.2,
                    label.size = 0)+
  ggtitle('After batch')+
  theme(legend.position = 'right')


# #聚类
library(clustree)
sce <- FindNeighbors(sce, dims = 1:30)
sce <- FindClusters(
  object = sce,
  #resolution = c(seq(.2,1,.2))
  resolution = 0.5
)
colnames(sce@meta.data)
length(table(sce@meta.data$seurat_clusters))
#17

# pdf('results/04.scRNA/clust.snn_res.pdf',he=15,wi=15)
# clustree(sce@meta.data, prefix = "integrated_snn_res.")
# dev.off()

######细胞注释#######
#Monocytes/macrophages  0,2,4,6,7,8,11
# VlnPlot(sce,features = c('FCGR3A','CD14','FCGR1A','CD68','CD163','MSR1'),pt.size = 0,group.by = 'seurat_clusters')
VlnPlot(sce,features = c('FCGR1A','CD68'),pt.size = 0,group.by = 'seurat_clusters') 
#T_cells	10
# VlnPlot(sce,features = c('CD3D','CD3G','CD2','TRAC','TRBC1','TRBC2'),pt.size = 0,group.by = 'seurat_clusters')
VlnPlot(sce,features = c('CD3D','CD2','TRAC'),pt.size = 0,group.by = 'seurat_clusters')

#Myeloid dendritic cell 	3,16

# VlnPlot(sce,features = c('CD1A','CD1C','CD207','CCL17','CCL22'),pt.size = 0,group.by = 'seurat_clusters')
VlnPlot(sce,features = c('CD1A','CCL22'),pt.size = 0,group.by = 'seurat_clusters')

#B_cells   9
# VlnPlot(sce,features = c('CD79A','IGHG1','IGHM','IGHG3'),pt.size = 0,group.by = 'seurat_clusters')
VlnPlot(sce,features = c('CD79A'),pt.size = 0,group.by = 'seurat_clusters')

#Mast_cells	12
# VlnPlot(sce,features = c('TPSAB1','TPSB2','CPA3','MS4A2'),pt.size = 0,group.by = 'seurat_clusters')
VlnPlot(sce,features = c('TPSAB1'),pt.size = 0,group.by = 'seurat_clusters')

##########肺上皮细胞
# #Ciliated_cells	纤毛细胞  
# VlnPlot(sce,features = c('CAPS'),pt.size = 0,group.by = 'seurat_clusters')

#Alveolar epithelial cell   15
VlnPlot(sce,features = c('CAV1'),pt.size = 0,group.by = 'seurat_clusters')

#cancer cell 1,5,13,14
VlnPlot(sce,features = c('EPCAM'),pt.size = 0,group.by = 'seurat_clusters')

# #CAV1 VWF CLDN5 CD34 MYLK RAMP2
# VlnPlot(sce,features = c('PRF1','RPF1','STAT1','TAGAP','THEMIS','TIGIT','TNFRSF17','TRBC1','TRBC2'),
#         pt.size = 0,group.by = 'seurat_clusters')

feat_gene<-c('FCGR1A','CD68',
             'CD3D','CD2','TRAC',
             'CD1A','CCL22','CD79A','TPSAB1',
             'CAV1','EPCAM')           
length(feat_gene)
#11

#marker基因的注释
library(reshape2)
vln.df=as.data.frame(sce[["RNA"]]@data[feat_gene,])
vln.df$gene=rownames(vln.df)
vln.df=melt(vln.df,id="gene")
colnames(vln.df)[c(2,3)]=c("CB","exp")

anno=sce@meta.data[,"seurat_clusters",drop=F]
anno$CB=rownames(sce@meta.data)
vln.df=inner_join(vln.df,anno,by="CB")
vln.df$gene=factor(vln.df$gene,levels = feat_gene) #为了控制画图的基因顺序
table(vln.df$seurat_clusters)

pdf('results/04.scRNA/FigS2.pdf',height = 15,width = 15)
vln.df%>%ggplot(aes(gene,exp))+geom_violin(aes(fill=seurat_clusters),scale = "width")+
  facet_grid(vln.df$seurat_clusters~.,scales = "free_y")+
  #scale_fill_brewer(palette = c("Set1"),direction = 1)+
  scale_x_discrete("")+scale_y_continuous("")+
  theme_bw()+
  theme(
    axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    legend.position = "none")
dev.off()

pdf('results/04.scRNA/FigS3.pdf',he=15,width =16)
FeaturePlot(sce,
            features = feat_gene,
            pt.size = 0.1,reduction = 'tsne',ncol = 4)
dev.off()



cell_anno=readMatrix('results/04.scRNA/cell_anno.txt',row = F,header = T)

head(cell_anno)
table(cell_anno$cell_type)

sce$cell_type<-sce$seurat_clusters
for (i in 1:nrow(cell_anno)){
  sce$cell_type=gsub(paste0('^',cell_anno$seurat_clusters[i],'$'),as.character(cell_anno$cell_type[i]),sce$cell_type)
}
table(sce$cell_type)
#length(feat_gene)

sce$cell_type2<-sce$seurat_clusters
for (i in 1:nrow(cell_anno)){
  sce$cell_type2=gsub(paste0('^',cell_anno$seurat_clusters[i],'$'),as.character(cell_anno$cell_type2[i]),sce$cell_type2)
}
table(sce$cell_type2)
colnames(sce@meta.data)

Idents(sce)='cell_type'
sc_marker_dotplot <- DotPlot(object = sce, features = feat_gene,
                             cols=c("blue", "red"),scale = T)+
  RotatedAxis()+ ggtitle("Marker Genes")+
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab('')+ylab('')
sc_marker_dotplot



fig4a=DimPlot(sce,group.by = 'seurat_clusters',
              reduction="tsne",
              cols = c(brewer.pal(9, "Set1"),brewer.pal(8, "Accent"),brewer.pal(9, "Paired")),
              label = "T", 
              pt.size = 0.2,
              label.size = 5)
fig4a
fig4b=DimPlot(sce,group.by = 'cell_type',
              reduction="tsne",
              cols =pal_d3()(8),
              label = "F", 
              pt.size = 0.1,
              label.size = 3.5)
fig4b

fig4c=DimPlot(sce,group.by = 'cell_type2',
              reduction="tsne",
              cols =pal_d3()(8),
              label = "F", 
              pt.size = 0.1,
              label.size = 3.5)
fig4c



#寻找差异基因时的差异倍数
Logfc = 0.5
#差异基因时最小的表达比例
Minpct = 0.25
DefaultAssay(sce) <- "RNA"
Idents(sce)<-'cell_type'

sce.markers <- FindAllMarkers(object = sce,logfc.threshold = Logfc, min.pct = Minpct,only.pos = T)
head(sce.markers)
sce.markers["pct.diff"]=sce.markers$pct.1-sce.markers$pct.2
sce.markers <- sce.markers[sce.markers$p_val_adj<0.05,]
table(sce.markers$cluster)
length(unique(sce.markers$gene))#762
head(sce.markers)
writeMatrix(sce.markers,'results/04.scRNA/scRNA_marker_gene.txt')
### 选择前5个marker基因
Top10 <- sce.markers %>% group_by(cluster) %>% slice_max(n =5, order_by = avg_logFC)  
# Top5 <- sce.markers %>% group_by(cluster) %>% slice_max(n =5, order_by = avg_logFC) 

sce1 <-ScaleData(sce, features =  rownames(sce))
table(sce1$cell_type)
sce1 <- subset(sce1, downsample = 200)

length(Top10$gene)
length(unique(Top10$gene))
###差异基因热图
marker.heatmap1=DotPlot(object = sce1, features = unique(Top10$gene),
                        cols=c("blue", "red"),scale = T)+ 
  RotatedAxis()+ ggtitle("Marker Genes")+ 
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab('')+ylab('')
marker.heatmap1
marker.heatmap2=DoHeatmap(sce1, features = Top10$gene,label=F,angle = 0) 
marker.heatmap2

#富集分析
library(clusterProfiler)
library(org.Hs.eg.db)
ids=bitr(sce.markers$gene,'SYMBOL','ENTREZID','org.Hs.eg.db') ## 将SYMBOL转成ENTREZID
sce.markers2=merge(sce.markers,ids,by.x='gene',by.y='SYMBOL')

## 函数split()可以按照分组因子，把向量，矩阵和数据框进行适当的分组。
## 它的返回值是一个列表，代表分组变量每个水平的观测。
gcSample=split(sce.markers2$ENTREZID, sce.markers2$cluster)
save(gcSample,file='results/04.scRNA/gcSample.RData')

load('results/04.scRNA/cluster_enrich.RData')
# sce.markers.enrich.res <- compareCluster(geneClusters = gcSample,
#                                          fun = "enrichKEGG",
#                                          organism = "hsa", pvalueCutoff = 0.05)
# 
# sce_cluster_enrich=sce.markers.enrich.res@compareClusterResult
cluster_enrichment_dotplot=enrichplot::dotplot(sce.markers.enrich.res)+
  ggtitle("KEGG Enrichment Analysis")+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
cluster_enrichment_dotplot
head(sce_cluster_enrich)
write.csv(sce_cluster_enrich,'results/04.scRNA/sce_cluster_enrich.csv',row.names = F)

save(sce,file='results/04.scRNA/sce.RData')



fig4=mg_merge_plot(mg_merge_plot(sc_marker_dotplot,fig4b,labels = c('A','B')),
                   marker.heatmap2,cluster_enrichment_dotplot,
                   nrow=3,heights = c(1,1.6,2),labels = c('','C','D'))
savePDF('results/04.scRNA/Fig4.pdf',fig4,height = 20,width = 15)