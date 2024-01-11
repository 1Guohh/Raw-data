#####基于bulk RNA-seq数据，对肿瘤样本进行基于GSVA算法的GO/KEGG富集分析######
dir.create('results/04.ARS_TME')
h.gmt <- read.gmt("h.all.v7.5.1.symbols.gmt")
hallmark_geneset=read.xlsx('HALLMARK_geneset.xlsx',check.names = F)
table(hallmark_geneset$Process.category)
hall_immu_pathway=hallmark_geneset[which(hallmark_geneset$Process.category=='immune'),'Hallmark.name']
hall_immu_pathway

immu.hall.genesets.list=split(x= h.gmt, f=h.gmt$ont)
immu.hall.genesets.list=sapply(immu.hall.genesets.list, function(x){subset(x,select='gene',drop=TRUE)})
names(immu.hall.genesets.list)=gsub('HALLMARK_','',names(immu.hall.genesets.list))
names(immu.hall.genesets.list)
h.immu.geneset=immu.hall.genesets.list[hall_immu_pathway]

# tcga.immu.hall.ssGSEA <- gsva(as.matrix(tcga_tpm_log_T),
#                          h.immu.geneset,
#                          method='ssgsea',
#                          min.sz=10,
#                          max.sz=500,
#                          verbose=TRUE)
# save(tcga.immu.hall.ssGSEA,file = 'results/tcga.immu.hall.ssGSEA.RData')
load('results/tcga.immu.hall.ssGSEA.RData')
rownames(tcga.immu.hall.ssGSEA)=gsub('HALLMARK_','',rownames(tcga.immu.hall.ssGSEA))
tcga.immu.hall.ssGSEA[1:5,1:5]
dim(tcga.immu.hall.ssGSEA)
my_mutiboxplot(dat = t(tcga.immu.hall.ssGSEA[,tcga.risktype.cli$Samples]),
               group = tcga.risktype.cli$Risktype )


kegg_gmt=list.files('origin_datas/immune_gmt/',pattern = '.gmt')
kegg_gmt_list=read.gmt(paste0('origin_datas/immune_gmt/',kegg_gmt[1]))
dim(kegg_gmt_list)
#88 2
for (f in 1:length(kegg_gmt)) {
  gmt=read.gmt(paste0('origin_datas/immune_gmt/',kegg_gmt[f]))
  kegg_gmt_list=rbind(kegg_gmt_list,gmt)
}
kegg_gmt_list=kegg_gmt_list[-c(1:88),]
table(kegg_gmt_list$ont)
# 
kegg.genesets.list=split(x=kegg_gmt_list,f=kegg_gmt_list$ont)
kegg.genesets.list=sapply(kegg.genesets.list, function(x){subset(x,select='gene',drop=TRUE)})
# tcga.kegg.ssGSEA=ssGSEAScore_by_muti_group_genes(gene.exp = tcga_tpm_log_T
#                                                 ,genelist = kegg.genesets.list)
# save(tcga.kegg.ssGSEA,file='results/tcga.kegg.ssGSEA.RData')
load('results/tcga.kegg.ssGSEA.RData')
tcga.kegg.ssGSEA[1:5,1:5]
dim(tcga.kegg.ssGSEA)
rownames(tcga.kegg.ssGSEA)=gsub('KEGG_','',rownames(tcga.kegg.ssGSEA))
my_mutiboxplot(dat = t(tcga.kegg.ssGSEA[,tcga.risktype.cli$Samples]),
               group = tcga.risktype.cli$Risktype )


cli_anno=data.frame(Risktype=tcga.risktype.cli[order(tcga.risktype.cli$Riskscore),'Risktype'])
rownames(cli_anno)=tcga.risktype.cli[order(tcga.risktype.cli$Riskscore),'Samples']
head(cli_anno)


risktype.col.use=risktype.col
names(risktype.col.use)=c('High','Low')

pdf('results/07.ARS_immu_pathway/TCGA_Hallmark_gsva.pdf',height = 3,width = 8)
pheatmap(tcga.immu.hall.ssGSEA[,rownames(cli_anno)],
         name = 'ssGSEA score',scale = 'row',
         color =  colorRampPalette(c("navy", "white", "red"))(100),
         legend_breaks = c(-3,0,3),
         annotation_col = cli_anno,
         annotation_colors = list(Risktype=risktype.col.use),
         cluster_rows = F,cluster_cols = F,
         show_rownames = T,show_colnames = F)
dev.off()


#####ARS与免疫微环境########
load('results/tcga.est.RData')
fig6a=my_mutiboxplot(dat = tcga.est[tcga.risktype.cli$Samples,],
                     group = tcga.risktype.cli$Risktype,group_cols = risktype.col,
                     fill='Risktype',angle = 0,hjust = 0.5,size = 12)
fig6a


tcga.timer=immu_timer(tcga_tpm_log_T)
fig6b=my_mutiboxplot(dat = tcga.timer[tcga.risktype.cli$Samples,],
                     group = tcga.risktype.cli$Risktype,group_cols = risktype.col,
                     fill='Risktype',angle = 0,hjust = 0.5,size = 12)
fig6b
load('results/tcga.immu.ssgsea.RData')
fig6c=my_mutiboxplot(dat = tcga.immu.ssgsea[tcga.risktype.cli$Samples,],
                     group = tcga.risktype.cli$Risktype,fill='Risktype',
                     group_cols = risktype.col,
                     legend.position = 'none',size = 12)
fig6c
load('results/tcga.ciber.RData')
fig6d=my_mutiboxplot(dat = tcga.ciber[tcga.risktype.cli$Samples,1:22],
                     group = tcga.risktype.cli$Risktype,group_cols = risktype.col,
                     fill='Risktype',legend.position = 'none',size = 12)
fig6d


tcga.mcp=immu_MCPcounter(tcga_tpm_log_T)
my_mutiboxplot(dat = tcga.mcp[tcga.risktype.cli$Samples,],group = tcga.risktype.cli$Risktype )

fig6=mg_merge_plot(mg_merge_plot(fig6a,fig6b,widths = c(1,1.5),labels = c('A','B'),common.legend = T),
                   fig6c,fig6d,nrow=3,labels = c('','C','D'))

savePDF('results/07.ARS_immu_pathway/Fig6.pdf',fig6,height = 16,width = 20)


#####基于scRNA-seq数据，验证GO/KEGG富集分析########
Idents(tumor_sc)='AS_group'
table(Idents(tumor_sc))
expr <- AverageExpression(tumor_sc, assays = "RNA", slot = "data")[[1]]
expr <- expr[rowSums(expr)>0,]  #选取非零基因
expr <- as.matrix(expr)
dim(expr)
sc.gsva.res.1 <- gsva(expr, h.immu.geneset, method="ssgsea")
# gsva.df <- data.frame(Genesets=rownames(sc.gsva.res), sc.gsva.res, check.names = F)
pdf('results/07.ARS_immu_pathway/sc_hallmark.pdf',height = 4,width = 8)
pheatmap(sc.gsva.res.1, show_colnames = T, 
         scale = "row",name = 'score')
dev.off()

