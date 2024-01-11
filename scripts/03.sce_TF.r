#########转录因子
##########05.基于scRNA-seq数据，计算ARS在细胞中的分布#############
load('results/04.scRNA/sce.RData')
colnames(sce@meta.data)
sce <- AddModuleScore(sce,
                      features = list(tcga.module.risk$module.gene),
                      ctrl = 10,
                      name = "AS")
colnames(sce@meta.data)
colnames(sce@meta.data)[13] <- 'AS'
range(sce@meta.data$AS)
sce$ARS=ifelse(sce$AS<0,'Low','High')
VlnPlot(sce,features = 'AS', pt.size = 0, adjust = 2,group.by = "cell_type2")

intersect(colnames(tcga.icg),rownames(sce@assays$RNA@data))
VlnPlot(sce,features = intersect(colnames(tcga.icg),rownames(sce@assays$RNA@data))[1:10],
        pt.size = 0, adjust = 2,group.by = "cell_type2")

##########提取肿瘤细胞
Idents(sce)<-'cell_type2'
tumor_sc=sce[, Idents(sce) %in% c('cancer cell')]
range(tumor_sc@meta.data$AS)



tumor_sc$AS_group=ifelse(tumor_sc$AS<0,'Low','High')
table(tumor_sc$AS_group)
dim(tumor_sc@meta.data)
colnames(tumor_sc@meta.data)


sc_df=data.frame(sce@meta.data, sce@reductions$tsne@cell.embeddings)
head(sc_df)
table(sc_df$AS_group)
dim(sc_df)
fig4d=ggplot(sc_df, aes(tSNE_1, tSNE_2, color=AS)) + 
  geom_point( size=1.5) + scale_color_viridis(option="A")  + 
  theme_light(base_size = 15)+#labs(title = "TNFA_SIGNALING_VIA_NFKB")+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+
  theme(plot.title = element_text(hjust = 0.5))
fig4d


##########06.转录因子在高和低ARS的差异激活#########
dir.create('results/05.sce_TF')
library(SCENIC)
library(RcisTarget)
cellInfo <- data.frame(tumor_sc@meta.data)
colnames(cellInfo)[which(colnames(cellInfo)=="orig.ident")] <- "sample"
colnames(cellInfo)[which(colnames(cellInfo)=="seurat_clusters")] <- "cluster"
colnames(cellInfo)[which(colnames(cellInfo)=="cell_type")] <- "cell_type"
colnames(cellInfo)[which(colnames(cellInfo)=="AS_group")] <- "AS_group"
cellInfo <- cellInfo[,c("sample","cluster","cell_type","AS_group")]
head(cellInfo)
table(cellInfo$AS_group)
saveRDS(cellInfo, file="int/cellInfo.Rds")



# #为了节省计算资源，随机抽取1000个细胞的数据子集
# set.seed(123)
# subcell <- sample(colnames(tumor_sc),1000)
# sce_tf <- tumor_sc[,subcell]
# exprMat <- as.matrix(sce_tf@assays$RNA@counts)
# dim(exprMat)
# exprMat[1:4,1:4] 
# colnames(sce_tf@meta.data)
# cellInfo <-  sce_tf@meta.data[,c(15,2,3)]
# colnames(cellInfo)=c('AS_group', 'nGene' ,'nUMI')
# head(cellInfo)
# table(cellInfo$AS_group)
# # High  Low 
# # 357 1713 
# # 保证cisTarget_databases 文件夹下面有下载好2个1G的文件
# dir.create("cisTarget_databases")
scenicOptions <- initializeScenic(org="hgnc",
                                  dbDir="cisTarget_databases", nCores=1)
# saveRDS(scenicOptions, file="int/scenicOptions.Rds") 
# 
# 
# ### Co-expression network
# #过滤标准是 基因表达量之和>细胞数*3%，且在1%的细胞中表达
# genesKept <- geneFiltering(exprMat, scenicOptions,
#                            minCountsPerGene = 10 * 0.01 * ncol(exprMat),
#                            minSamples = ncol(exprMat) * 0.15)
# exprMat_filtered <- exprMat[genesKept, ]
# exprMat_filtered[1:4,1:4]
# dim(exprMat_filtered)
# runCorrelation(exprMat_filtered, scenicOptions)
# range(exprMat_filtered)
# exprMat_filtered_log <- log2(exprMat_filtered+1)
# range(exprMat_filtered_log)
# dim(exprMat_filtered_log)
# ###这一步跑了10个小时
# runGenie3(exprMat_filtered_log, scenicOptions)


##推断共表达模块
runSCENIC_1_coexNetwork2modules(scenicOptions)
##推断转录调控网络（regulon）
runSCENIC_2_createRegulons(scenicOptions,coexMethod = c('top50','top50perTarget'))
#以上代码可增加参数coexMethod=c("w001", "w005", "top50", "top5perTarget", "top10perTarget", "top50perTarget"))
#默认6种方法的共表达网络都计算，可以少选几种方法以减少计算量

#==regulon活性评分与可视化==##
#regulons计算AUC值并进行下游分析
exprMat_all <- as.matrix(tumor_sc@assays$RNA@data)
range(exprMat_all)
dim(exprMat_all)
runSCENIC_3_scoreCells(scenicOptions, exprMat=exprMat_all)
library(rbokeh)
aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMat_all)
savedSelections <- shiny::runApp(aucellApp)


##导入原始regulonAUC矩阵
AUCmatrix <- readRDS("int/3.4_regulonAUC.Rds")
AUCmatrix <- AUCmatrix@assays@data@listData$AUC
AUCmatrix[1:5,1:5]
AUCmatrix <- data.frame(t(AUCmatrix), check.names=F)
AUCmatrix[1:5,1:5]
dim(AUCmatrix)

RegulonName_AUC <- colnames(AUCmatrix)
RegulonName_AUC <- gsub(' \\(','_',RegulonName_AUC)
RegulonName_AUC <- gsub('\\)','',RegulonName_AUC)
colnames(AUCmatrix) <- RegulonName_AUC
scRNAauc <- AddMetaData(tumor_sc, AUCmatrix)
scRNAauc@assays$integrated <- NULL
colnames(scRNAauc@meta.data)
saveRDS(scRNAauc,'scRNAauc.rds')


##利用Seurat可视化AUC
dir.create('scenic_seurat')
DefaultAssay(scRNAauc) <- "RNA"
Idents(scRNAauc)<-'AS_group'
table(Idents(scRNAauc))


library(pheatmap)
cellInfo <- readRDS("int/cellInfo.Rds")
celltype = subset(cellInfo,select = 'AS_group')
AUCmatrix <- t(AUCmatrix)


###########差异转录因子分布与富集分析#########

calcRSS <- function(AUC, cellAnnotation, cellTypes=NULL)
{
  if(any(is.na(cellAnnotation))) stop("NAs in annotation")
  if(any(class(AUC)=="aucellResults")) AUC <- getAUC(AUC)
  normAUC <- AUC/rowSums(AUC)
  if(is.null(cellTypes)) cellTypes <- unique(cellAnnotation)
  # 
  ctapply <- lapply
  if(require('BiocParallel')) ctapply <- bplapply
  
  rss <- ctapply(cellTypes, function(thisType)
    sapply(rownames(normAUC), function(thisRegulon)
    {
      pRegulon <- normAUC[thisRegulon,]
      pCellType <- as.numeric(cellAnnotation==thisType)
      pCellType <- pCellType/sum(pCellType)
      .calcRSS.oneRegulon(pRegulon, pCellType)
    })
  )
  rss <- do.call(cbind, rss)
  colnames(rss) <- cellTypes
  return(rss)
}


plotRSS <- function(rss, labelsToDiscard=NULL, zThreshold=1,
                    cluster_columns=FALSE, order_rows=TRUE, thr=0.01, varName="cellType",
                    col.low="grey90", col.mid="darkolivegreen3", col.high="darkgreen",
                    revCol=FALSE, verbose=TRUE)
{
  varSize="RSS"
  varCol="Z"
  if(revCol) {
    varSize="Z"
    varCol="RSS"
  }
  
  rssNorm <- scale(rss) # scale the full matrix...
  rssNorm <- rssNorm[,which(!colnames(rssNorm) %in% labelsToDiscard)] # remove after calculating...
  rssNorm[rssNorm < 0] <- 0
  
  ## to get row order (easier...)
  rssSubset <- rssNorm
  if(!is.null(zThreshold)) rssSubset[rssSubset < zThreshold] <- 0
  tmp <- .plotRSS_heatmap(rssSubset, thr=thr, cluster_columns=cluster_columns, order_rows=order_rows, verbose=verbose)
  rowOrder <- rev(tmp@row_names_param$labels)
  rm(tmp)
  
  
  ## Dotplot
  rss.df <- reshape2::melt(rss)
  head(rss.df)
  colnames(rss.df) <- c("Topic", varName, "RSS")
  rssNorm.df <- reshape2::melt(rssNorm)
  colnames(rssNorm.df) <- c("Topic", varName, "Z")
  rss.df <- base::merge(rss.df, rssNorm.df)
  
  rss.df <- rss.df[which(!rss.df[,varName] %in% labelsToDiscard),] # remove after calculating...
  if(nrow(rss.df)<2) stop("Insufficient rows left to plot RSS.")
  
  rss.df <- rss.df[which(rss.df$Topic %in% rowOrder),]
  rss.df[,"Topic"] <- factor(rss.df[,"Topic"], levels=rowOrder)
  p <- DoHeatmap (rss.df, 
                  var.x=varName, var.y="Topic", 
                  var.size=varSize, min.size=.5, max.size=5,
                  var.col=varCol, col.low=col.low, col.mid=col.mid, col.high=col.high)
  
  invisible(list(plot=p, df=rss.df, rowOrder=rowOrder))
}

#' @aliases plotRSS
#' @import ggplot2
#' @importFrom ggrepel geom_label_repel
#' @export 
plotRSS_oneSet <- function(rss, setName, n=5)
{
  library(ggplot2)
  library(ggrepel)
  
  rssThisType <- sort(rss[,setName], decreasing=TRUE)
  thisRss <- data.frame(regulon=names(rssThisType), rank=seq_along(rssThisType), rss=rssThisType)
  thisRss$regulon[(n+1):nrow(thisRss)] <- NA
  
  ggplot(thisRss, aes(x=rank, y=rss)) + 
    geom_point(color = "blue", size = 1) + 
    ggtitle(setName) + 
    geom_label_repel(aes(label = regulon),
                     box.padding   = 0.35, 
                     point.padding = 0.5,
                     segment.color = 'grey50',
                     na.rm=TRUE) +
    theme_classic()
}



## Internal functions:
.H <- function(pVect){
  pVect <- pVect[pVect>0] # /sum(pVect) ??
  - sum(pVect * log2(pVect))
}

# Jensen-Shannon Divergence (JSD)
calcJSD <- function(pRegulon, pCellType)
{
  (.H((pRegulon+pCellType)/2)) - ((.H(pRegulon)+.H(pCellType))/2)
}

# Regulon specificity score (RSS)
.calcRSS.oneRegulon <- function(pRegulon, pCellType)
{
  jsd <- calcJSD(pRegulon, pCellType)
  1 - sqrt(jsd)
}

.plotRSS_heatmap <- function(rss, thr=NULL, row_names_gp=gpar(fontsize=5), order_rows=TRUE, cluster_rows=FALSE, name="RSS", verbose=TRUE, ...)
{
  if(is.null(thr)) thr <- signif(quantile(rss, p=.97),2)
  
  library(ComplexHeatmap)
  rssSubset <- rss[rowSums(rss > thr)>0,]
  rssSubset <- rssSubset[,colSums(rssSubset > thr)>0]
  
  if(verbose) message("Showing regulons and cell types with any RSS > ", thr, " (dim: ", nrow(rssSubset), "x", ncol(rssSubset),")")
  
  if(order_rows)
  {
    maxVal <- apply(rssSubset, 1, which.max)
    rss_ordered <- rssSubset[0,]
    for(i in 1:ncol(rssSubset))
    {
      tmp <- rssSubset[which(maxVal==i),,drop=F]
      tmp <- tmp[order(tmp[,i], decreasing=FALSE),,drop=F]
      rss_ordered <- rbind(rss_ordered, tmp)
    }
    rssSubset <- rss_ordered
    cluster_rows=FALSE
  }
  
  Heatmap(rssSubset, name=name, row_names_gp=row_names_gp, cluster_rows=cluster_rows, ...)
} 


#####计算RSS与CSI########
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellInfo[colnames(regulonAUC), "AS_group"], )
dim(rss)
head(rss)
write.csv(rss,'results/05.sce_TF/riskgroup_RSS.CSV')
rss=crbind2DataFrame(rss)
rss=data.frame(TF=rownames(rss),rss)
rss=melt(rss)
head(rss)
table(rss$variable)
rss.high=rss[which(rss$variable=='High'),]
rss.high.order=rss.high[order(rss.high$value),]

fp_dat_high=data.frame(Samples=rev(1:nrow(rss.high)),rss.high.order)
tail(fp_dat_high)
fig5a=ggplot(fp_dat_high,aes(x=Samples,y=value))+
  geom_point(size=.5)+ylab('Regulon specific score(RSS)')+xlab('Regulons')+
  geom_point(aes(1, 0.1874943),color="red",size=.5) +
  geom_point(aes(2, 0.1872287),color="red",size=.5) +
  geom_point(aes(3, 0.1867311),color="red",size=.5) +
  geom_point(aes(4, 0.1862303),color="red",size=.5) +
  geom_point(aes(5, 0.1860812),color="red",size=.5) +
  ggtitle('High ARS')
fig5a

rss.low=rss[which(rss$variable=='Low'),]
rss.low.order=rss.low[order(rss.low$value),]

fp_dat_low=data.frame(Samples=rev(1:nrow(rss.low)),rss.low.order)
tail(fp_dat_low)
fig5b=ggplot(fp_dat_low,aes(x=Samples,y=value))+
  geom_point(size=.5)+ylab('Regulon specific score(RSS)')+xlab('Regulons')+
  geom_point(aes(1, 0.7066729),color="red",size=.5) +
  geom_point(aes(2, 0.7054738),color="red",size=.5) +
  geom_point(aes(3,  0.7029288),color="red",size=.5) +
  geom_point(aes(4,  0.7029288),color="red",size=.5) +
  geom_point(aes(5,  0.7029288),color="red",size=.5) +
  ggtitle('Low ARS')
fig5b


#####CSI计算
calculate_csi <- function(regulonAUC,
                          calc_extended = FALSE,
                          verbose = FALSE){
  
  compare_pcc <- function(vector_of_pcc,pcc){
    pcc_larger <- length(vector_of_pcc[vector_of_pcc > pcc])
    if(pcc_larger == length(vector_of_pcc)){
      return(0)
    }else{
      return(length(vector_of_pcc))
    }
  }
  
  calc_csi <- function(reg,reg2,pearson_cor){
    test_cor <- pearson_cor[reg,reg2]
    total_n <- ncol(pearson_cor)
    pearson_cor_sub <- subset(pearson_cor,rownames(pearson_cor) == reg | rownames(pearson_cor) == reg2)
    
    sums <- apply(pearson_cor_sub,MARGIN = 2, FUN = compare_pcc, pcc = test_cor)
    fraction_lower <- length(sums[sums == nrow(pearson_cor_sub)]) / total_n
    return(fraction_lower)
  }
  
  regulonAUC_sub <- regulonAUC@assays@data@listData$AUC
  
  if(calc_extended == TRUE){
    regulonAUC_sub <- subset(regulonAUC_sub,grepl("extended",rownames(regulonAUC_sub)))
  } else if (calc_extended == FALSE){
    regulonAUC_sub <- subset(regulonAUC_sub,!grepl("extended",rownames(regulonAUC_sub)))
  }
  
  regulonAUC_sub <- t(regulonAUC_sub)
  
  pearson_cor <- cor(regulonAUC_sub)
  pearson_cor_df <- as.data.frame(pearson_cor)
  pearson_cor_df$regulon_1 <- rownames(pearson_cor_df)
  pearson_cor_long <- pearson_cor_df %>%
    gather(regulon_2,pcc,-regulon_1) %>%
    mutate("regulon_pair" = paste(regulon_1,regulon_2,sep="_"))
  
  
  regulon_names <- unique(colnames(pearson_cor))
  num_of_calculations <- length(regulon_names)*length(regulon_names)
  
  csi_regulons <- data.frame(matrix(nrow=num_of_calculations,ncol = 3))
  
  colnames(csi_regulons) <- c("regulon_1",
                              "regulon_2",
                              "CSI")
  
  num_regulons <- length(regulon_names)
  
  f <- 0
  for(reg in regulon_names){
    ## Check if user wants to print info
    if(verbose == TRUE){
      print(reg)
    }
    for(reg2 in regulon_names){
      f <- f + 1
      
      fraction_lower <- calc_csi(reg,reg2,pearson_cor)
      
      csi_regulons[f,] <- c(reg,reg2,fraction_lower)
      
    }
  }
  csi_regulons$CSI <- as.numeric(csi_regulons$CSI)
  return(csi_regulons)
}


csi <- calculate_csi(regulonAUC, calc_extended = FALSE, verbose = FALSE)
head(csi)

table(csi$regulon_1)
csi_dat=dcast(csi,regulon_1~regulon_2)
csi_dat[1:5,1:5]
rownames(csi_dat)=csi_dat[,1]
csi_dat=csi_dat[,-1]
csi_dat[1:5,1:5]
dim(csi_dat)
write.csv(csi_dat,'results/05.sce_TF/Regulon_CSI.csv')

###特异性连接指数热图
fig5c=pheatmap(csi_dat,border_color =NA,name = 'CSI',
               #cutree_cols = 3,cutree_rows = 3,
               show_colnames = F,
               legend_labels = c("low","","high"),
               fontsize_row = 12
)
fig5c=as.ggplot(fig5c)
fig5c

fig5=mg_merge_plot(fig5c,mg_merge_plot(fig5a,fig5b,labels = c('B','C')),
                   nrow=2,heights = c(1.5,1),labels = c('A',''))
savePDF('results/05.sce_TF/Fig5.pdf',fig5,height = 20,width = 17)

############转录因子富集分析####
tf_regulon=readRDS("int/3.1_regulons_forAUCell.Rds")
names(tf_regulon)
class(tf_regulon)
write.csv(tf_regulon,'results/05.sce_TF/regulon_geneset.csv')

output_gmt<-function(geneset,file){
  sink(file) #将输出结果重定向到file
  lapply(names(geneset),function(i){
    #按照列名，将要连接的变量转化为向量型，用制表符连接成一个字符串
    cat(paste(c(i,'tmp',geneset[[i]]),collapse='\t'))
    cat('\n') #输出后新起一行
  })
  sink() #结束重定向
}
output_gmt(tf_regulon,'results/05.sce_TF/regulon_geneset.gmt')
library(topGO)
library(clusterProfiler)
library(Rgraphviz)
library(org.Hs.eg.db)


#high组
########"JUN_extended (47g)"
h1_enrich_res=mg_clusterProfiler(genes = tf_regulon[[tail(fp_dat_high)[5,'TF']]],keytype = 'SYMBOL')
head(h1_enrich_res$KEGG)
h1_enrich_res_GO=rbind(h1_enrich_res$GO_BP@result,h1_enrich_res$GO_CC@result,h1_enrich_res$GO_MF@result)
h1_enrich_res_GO$DB=c(rep('BP',nrow(h1_enrich_res$GO_BP@result)),
                      rep('CC',nrow(h1_enrich_res$GO_CC@result)),
                      rep('MF',nrow(h1_enrich_res$GO_MF@result)))

h1_enrich_res_GO.fit =h1_enrich_res_GO %>% group_by(DB) %>% slice_head(n =5)
table(h1_enrich_res_GO.fit$DB)


high.tf.p1=list()
high.tf.p1[[1]]=enrichplot::dotplot(h1_enrich_res$KEGG)+ggtitle('KEGG')+
  theme(legend.position = c(1, 0),legend.justification = c(1,0))
high.tf.p1[[2]]=ggplot(data = h1_enrich_res_GO.fit, aes(x = Count, y = Description))+
  facet_grid(DB ~ ., scales = "free", space = "free", margins = F)+
  geom_point(aes(size = Count,color = -log10(p.adjust))) +
  scale_color_gradient(low = 'blue', high = 'red')+theme_bw()+
  theme(axis.text.y = element_text(size = 10))


######"FOS_extended (44g)"
h2_enrich_res=mg_clusterProfiler(genes = tf_regulon[[tail(fp_dat_high)[4,'TF']]],keytype = 'SYMBOL')
head(h2_enrich_res$KEGG)
h2_enrich_res_GO=rbind(h2_enrich_res$GO_BP@result,h2_enrich_res$GO_CC@result,h2_enrich_res$GO_MF@result)
h2_enrich_res_GO$DB=c(rep('BP',nrow(h2_enrich_res$GO_BP@result)),
                      rep('CC',nrow(h2_enrich_res$GO_CC@result)),
                      rep('MF',nrow(h2_enrich_res$GO_MF@result)))

h2_enrich_res_GO.fit =h2_enrich_res_GO %>% group_by(DB) %>% slice_head(n =5)
table(h2_enrich_res_GO.fit$DB)


high.tf.p2=list()
high.tf.p2[[1]]=enrichplot::dotplot(h2_enrich_res$KEGG)+ggtitle('KEGG')+
  theme(legend.position = c(1, 0),legend.justification = c(1,0))
high.tf.p2[[2]]=ggplot(data = h2_enrich_res_GO.fit, aes(x = Count, y = Description))+
  facet_grid(DB ~ ., scales = "free", space = "free", margins = F)+
  geom_point(aes(size = Count,color = -log10(p.adjust))) +
  scale_color_gradient(low = 'blue', high = 'red')+theme_bw()+
  theme(axis.text.y = element_text(size = 10))






######"JUNB_extended (154g)"
h3_enrich_res=mg_clusterProfiler(genes = tf_regulon[[tail(fp_dat_high)[2,'TF']]],keytype = 'SYMBOL')
head(h3_enrich_res$KEGG)

h3_enrich_res_GO=rbind(h3_enrich_res$GO_BP@result,h3_enrich_res$GO_CC@result,h3_enrich_res$GO_MF@result)
h3_enrich_res_GO$DB=c(rep('BP',nrow(h3_enrich_res$GO_BP@result)),
                      rep('CC',nrow(h3_enrich_res$GO_CC@result)),
                      rep('MF',nrow(h3_enrich_res$GO_MF@result)))

h3_enrich_res_GO.fit =h3_enrich_res_GO %>% group_by(DB) %>% slice_head(n =5)
table(h3_enrich_res_GO.fit$DB)


high.tf.p3=list()
high.tf.p3[[1]]=enrichplot::dotplot(h3_enrich_res$KEGG)+ggtitle('KEGG')+
  theme(legend.position = c(1, 0),legend.justification = c(1,0))
high.tf.p3[[2]]=ggplot(data = h3_enrich_res_GO.fit, aes(x = Count, y = Description))+
  facet_grid(DB ~ ., scales = "free", space = "free", margins = F)+
  geom_point(aes(size = Count,color = -log10(p.adjust))) +
  scale_color_gradient(low = 'blue', high = 'red')+theme_bw()+
  theme(axis.text.y = element_text(size = 10))




#low组
### "BCLAF1_extended (2226g)"
l1_enrich_res=mg_clusterProfiler(genes = tf_regulon[[tail(fp_dat_low)[6,'TF']]],keytype = 'SYMBOL')
head(l1_enrich_res$KEGG)
l1_enrich_res_GO=rbind(l1_enrich_res$GO_BP@result,l1_enrich_res$GO_CC@result,l1_enrich_res$GO_MF@result)
l1_enrich_res_GO$DB=c(rep('BP',nrow(l1_enrich_res$GO_BP@result)),
                      rep('CC',nrow(l1_enrich_res$GO_CC@result)),
                      rep('MF',nrow(l1_enrich_res$GO_MF@result)))

l1_enrich_res_GO.fit =l1_enrich_res_GO %>% group_by(DB) %>% slice_head(n =5)
table(l1_enrich_res_GO.fit$DB)
l1_enrich_res_GO.fit=crbind2DataFrame(l1_enrich_res_GO.fit)

low.tf.p1=list()
low.tf.p1[[1]]=enrichplot::dotplot(l1_enrich_res$KEGG)+ggtitle('KEGG')+
  theme(legend.position = c(1, 0),legend.justification = c(1,0))
head(l1_enrich_res_GO.fit)
low.tf.p1[[2]]=ggplot(data = l1_enrich_res_GO.fit, aes(x = Count, y = Description))+
  facet_grid(DB ~ ., scales = "free", space = "free", margins = F)+
  geom_point(aes(size = Count,color = -log10(p.adjust))) +
  scale_color_gradient(low = 'blue', high = 'red')+theme_bw()+
  theme(axis.text.y = element_text(size = 10))



#####"BCLAF1 (956g)"
l2_enrich_res=mg_clusterProfiler(genes = tf_regulon[[tail(fp_dat_low)[5,'TF']]],keytype = 'SYMBOL')
head(l2_enrich_res$KEGG)

l2_enrich_res_GO=rbind(l2_enrich_res$GO_BP@result,l2_enrich_res$GO_CC@result,l2_enrich_res$GO_MF@result)
l2_enrich_res_GO$DB=c(rep('BP',nrow(l2_enrich_res$GO_BP@result)),
                      rep('CC',nrow(l2_enrich_res$GO_CC@result)),
                      rep('MF',nrow(l2_enrich_res$GO_MF@result)))

l2_enrich_res_GO.fit =l2_enrich_res_GO %>% group_by(DB) %>% slice_head(n =5)
table(l2_enrich_res_GO.fit$DB)


low.tf.p2=list()
low.tf.p2[[1]]=enrichplot::dotplot(l2_enrich_res$KEGG)+ggtitle('KEGG')+
  theme(legend.position = c(1, 0),legend.justification = c(1,0))
low.tf.p2[[2]]=ggplot(data = l2_enrich_res_GO.fit, aes(x = Count, y = Description))+
  facet_grid(DB ~ ., scales = "free", space = "free", margins = F)+
  geom_point(aes(size = Count,color = -log10(p.adjust))) +
  scale_color_gradient(low = 'blue', high = 'red')+theme_bw()+
  theme(axis.text.y = element_text(size = 10))




###"UQCRB (1160g)"
l3_enrich_res=mg_clusterProfiler(genes = tf_regulon[[tail(fp_dat_low)[4,'TF']]],keytype = 'SYMBOL')
head(l3_enrich_res$KEGG)
l3_enrich_res_GO=rbind(l3_enrich_res$GO_BP@result,l3_enrich_res$GO_CC@result,l3_enrich_res$GO_MF@result)
l3_enrich_res_GO$DB=c(rep('BP',nrow(l3_enrich_res$GO_BP@result)),
                      rep('CC',nrow(l3_enrich_res$GO_CC@result)),
                      rep('MF',nrow(l3_enrich_res$GO_MF@result)))

l3_enrich_res_GO.fit =l3_enrich_res_GO %>% group_by(DB) %>% slice_head(n =5)
table(l3_enrich_res_GO.fit$DB)

low.tf.p3=list()
low.tf.p3[[1]]=enrichplot::dotplot(l3_enrich_res$KEGG)+ggtitle('KEGG')+
  theme(legend.position = c(1, 0),legend.justification = c(1,0))
low.tf.p3[[2]]=ggplot(data = l3_enrich_res_GO.fit, aes(x = Count, y = Description))+
  facet_grid(DB ~ ., scales = "free", space = "free", margins = F)+
  geom_point(aes(size = Count,color = -log10(p.adjust))) +
  scale_color_gradient(low = 'blue', high = 'red')+theme_bw()+
  theme(axis.text.y = element_text(size = 10))


fig6_high_kegg=mg_merge_plot(high.tf.p1[[1]],high.tf.p2[[1]],high.tf.p3[[1]],ncol=1,nrow=3,labels = c(LETTERS[1:3]))
fig6_low_kegg=mg_merge_plot(low.tf.p1[[1]],low.tf.p2[[1]],low.tf.p3[[1]],ncol=1,nrow=3,labels = c(LETTERS[4:6]))

fig6_high_GO=mg_merge_plot(high.tf.p1[[2]],high.tf.p2[[2]],high.tf.p3[[2]],ncol=1,nrow=3)
fig6_low_GO=mg_merge_plot(low.tf.p1[[2]],low.tf.p2[[2]],low.tf.p3[[2]],ncol=1,nrow=3)

TF_high=mg_merge_plot(fig6_high_kegg,fig6_high_GO,ncol=2)
TF_low=mg_merge_plot(fig6_low_kegg,fig6_low_GO,ncol=2)

savePDF('results/05.sce_TF/TF_high.pdf',TF_high,height = 20,width = 22)
savePDF('results/05.sce_TF/TF_low.pdf',TF_low,height = 20,width = 22)

pdf('results/05.sce_TF/P714_FigS2.pdf',height = 40,width = 25,onefile = F)
mg_merge_plot(TF_high,TF_low,nrow=2)
dev.off()