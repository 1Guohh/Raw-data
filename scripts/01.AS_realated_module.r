
######01.建立AS风险模型#########
dir.create('results/01.AS_realated_module')
head(luad_as_cli)
write.csv(luad_as_cli,'results/01.AS_realated_module/TCGA_LUAD_cli.csv',quote = F)
#rownames(luad_as_cli)=luad_as_cli$Samples
range(luad_as_cli$AneuploidyScore)


#######肿瘤和正常组织差异基因分析
table(tcga_type$type)
tcga.limma=mg_limma_DEG(exp = tcga_exp[,tcga_type$Samples],group = tcga_type$type,ulab = 'Tumor',dlab = 'Normal')
tcga.limma$Summary
fig1a=my_volcano(dat = tcga.limma,p_cutoff = 0.05,fc_cutoff = 1)
fig1a
tcga.degs=tcga.limma$DEG[which(tcga.limma$DEG$adj.P.Val<0.05 & abs(tcga.limma$DEG$logFC)>1),]
write.csv(tcga.degs,'results/01.AS_realated_module/TCGA_DEGs.csv',quote = F)
dim(tcga.degs)
#3348

#计算差异基因与AS score相关性
AS_cor_gene=cbind.data.frame(AS=luad_as_cli$AneuploidyScore,
                             t(tcga_tpm_log_T[intersect(rownames(tcga.degs),rownames(tcga_tpm_log_T)),luad_as_cli$Samples]))
dim(AS_cor_gene)
AS_gene_cor_res=Hmisc::rcorr(as.matrix(AS_cor_gene),type = 'spearman')
AS_gene_cor_res$P[is.na(AS_gene_cor_res$P)] <- 0
AS_gene_cor_res$r[1:5,1:5]
AS_gene_cor_res_fit=data.frame(gene=rownames(AS_gene_cor_res$n[-1,]),
                               cor=AS_gene_cor_res$r[-1,1],
                               p=AS_gene_cor_res$P[-1,1])
head(AS_gene_cor_res_fit)
table(AS_gene_cor_res_fit$cor>0.3,AS_gene_cor_res_fit$p<0.05)
AS_gene_cor_res_fit=AS_gene_cor_res_fit[which(AS_gene_cor_res_fit$cor>0.3),]
dim(AS_gene_cor_res_fit)
#182
write.csv(AS_gene_cor_res_fit,'results/01.AS_realated_module/AS_cor_DEGs.csv',quote = F,row.names = F)

samples_anno=data.frame(Categore=tcga_type[order(tcga_type$type),'type'])
rownames(samples_anno)=tcga_type[order(tcga_type$type),'Samples']
head(samples_anno)
fig1b=pheatmap(tcga_exp[AS_gene_cor_res_fit$gene,rownames(samples_anno)],
               name = 'Gene Expression',
               scale = 'row',breaks = c(-3,0,3),
               annotation_col = samples_anno,
               cluster_cols = F,
               show_colnames = F,show_rownames = F)
library(ggplotify)
fig1b = as.ggplot(fig1b)
fig1b
#######AS相关基因功能富集分析
as_cor.enrich=mg_clusterProfiler(genes = rownames(AS_gene_cor_res_fit),keytype = 'SYMBOL')
fig1c=enrichplot::dotplot(as_cor.enrich$KEGG)+ggtitle('KEGG')
fig1d=enrichplot::dotplot(as_cor.enrich$GO_BP)+ggtitle('Biological Process')
fig1e=enrichplot::dotplot(as_cor.enrich$GO_CC)+ggtitle('Cellular Component')
fig1f=enrichplot::dotplot(as_cor.enrich$GO_MF)+ggtitle('Molecular Function')

fig1=mg_merge_plot(mg_merge_plot(fig1a,fig1b,labels = c('A','B')),
                   mg_merge_plot(fig1c,fig1d,fig1e,fig1f,ncol=2,nrow=2,labels = LETTERS[3:6]),
                   nrow=2,heights = c(1,2))
savePDF('results/01.AS_realated_module/Fig1.pdf',fig1,height = 18,width = 16)

AS_cor_DEGs_enrich=list()
AS_cor_DEGs_enrich$KEGG=as_cor.enrich$KEGG@result
AS_cor_DEGs_enrich$GO_BP=as_cor.enrich$GO_BP@result
AS_cor_DEGs_enrich$GO_CC=as_cor.enrich$GO_CC@result
AS_cor_DEGs_enrich$GO_MF=as_cor.enrich$GO_MF@result

write.xlsx(AS_cor_DEGs_enrich,'results/01.AS_realated_module/AS_cor_DEGs_enrich.xlsx',overwrite = T)

############TCGA###########
pre.genes=AS_gene_cor_res_fit$gene
length(pre.genes)
tcga_model_data=t(tcga_tpm_log_T[pre.genes,luad_as_cli$Samples])
colnames(tcga_model_data)=gsub('-','__',colnames(tcga_model_data))
tcga_model_data=merge(data.frame(Samples=luad_as_cli$Samples,OS=luad_as_cli$OS,OS.time=luad_as_cli$OS.time),
                      data.frame(Samples=rownames(tcga_model_data),tcga_model_data),by='Samples')
rownames(tcga_model_data)=tcga_model_data$Samples
tcga_model_data=tcga_model_data[,-1]
tcga_model_data=crbind2DataFrame(tcga_model_data)
dim(tcga_model_data)

tcga.cox=cox_batch(dat = tcga_tpm_log_T[pre.genes,luad_as_cli$Samples],time = luad_as_cli$OS.time,event = luad_as_cli$OS)
tcga.cox=na.omit(tcga.cox)
tcga.cox$coef=log(tcga.cox$HR)
head(tcga.cox)
rownames(tcga.cox)=gsub('-','__',rownames(tcga.cox))
p_cutoff=0.01
table(tcga.cox$p.value<p_cutoff)
tcga.cox.fit=tcga.cox[tcga.cox$p.value<p_cutoff,]
#125
write.csv(tcga.cox.fit,'results/01.AS_realated_module/tcga.cox.fit.csv')

tcga.lasso=get_riskscore.lasso(dat = tcga_model_data[luad_as_cli$Samples,rownames(tcga.cox.fit)],
                               os = tcga_model_data[luad_as_cli$Samples,]$OS,
                               os.time = tcga_model_data[luad_as_cli$Samples,]$OS.time)
length(tcga.lasso$lasso.gene)#5
fig2a=tcga.lasso$plot


tcga.module.risk=get_riskscore(dat = tcga_model_data[luad_as_cli$Samples,tcga.lasso$lasso.gene],
                               os = tcga_model_data[luad_as_cli$Samples,]$OS,
                               os.time = tcga_model_data[luad_as_cli$Samples,]$OS.time,
                               step = F,direction = c("both", "backward", "forward")[1])
length(tcga.module.risk$module.gene)
#5
tcga.module.risk$model


###fig2画coefficient bar图
cox <- coxph(as.formula(paste0("Surv(OS.time, OS) ~"
                               ,paste0(tcga.module.risk$module.gene,collapse = '+'))),
             data =as.data.frame(tcga_model_data))
coef(cox)
tcga.model.cox=data.frame(gene=names(coef(cox)),coef=as.numeric(coef(cox)))

fig2b=ggplot(tcga.model.cox,aes(x=reorder(gene, coef),y=coef,fill=coef)) + 
  geom_bar(stat="identity",width = 0.8) +
  theme(legend.position="none")+
  coord_flip()+xlab('Gene')+ylab('coefficient')
fig2b

pdf('results/01.AS_realated_module/gene_forest_muti.pdf',height =4,width = 8,onefile = F)
ggforest(cox,data = tcga_model_data)
dev.off()

tcga.cox.fit[names(coef(cox)),]

bioForest=function(rt=null,col){
  #读取输入文件
  gene <- rownames(rt)
  hr <- sprintf("%.3f",rt[,2])
  hrLow  <- sprintf("%.3f",rt[,3])
  hrHigh <- sprintf("%.3f",rt[,4])
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt[,1]<0.001, "<0.001", sprintf("%.3f", rt[,1]))
  
  #输出图形
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2.5))
  
  #绘制森林图左边的临床信息
  xlim = c(0,3)
  par(mar=c(4,2.5,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
  text(3.1,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3.1,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1)
  
  #绘制右边森林图
  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
  abline(v=1,col="black",lty=2,lwd=2)
  boxcolor = ifelse(as.numeric(hr) > 1, col[2], col[1])
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.5)
  axis(1)
}
pdf('results/01.AS_realated_module/gene_forest_sig.pdf',height =4,width = 8,onefile = F)
bioForest(tcga.cox.fit[names(coef(cox)),-5],col=c('red,blue'))
dev.off()

fig2c=ggplotTimeROC(time = tcga.module.risk$result$time/365
                    ,status = tcga.module.risk$result$status
                    ,score = tcga.module.risk$result$riskscore
                    ,mks = c(1,3,5))
fig2c
tcga.risktype.cli=data.frame(luad_as_cli,
                             Riskscore=tcga.module.risk$result[luad_as_cli$Samples,]$riskscorez)
tcga.risktype.cli$Risktype=ifelse(tcga.risktype.cli$Riskscore>0,'High','Low')
risktype.col=ggsci::pal_nejm()(10)[2:3]
tcga.km=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ Risktype,
                                data = tcga.risktype.cli),
                   data=tcga.risktype.cli,
                   conf.int = F,pval = T,fun = "pct",risk.table = T, size = 0.7,
                   title='TCGA-LUAD',ggtheme=custom_theme(),
                   linetype = c("solid", "dashed","strata")[1],
                   palette = risktype.col,
                   legend = c('top', 'bottom', 'left', 'right', 'none')[1],
                   #legend = c(0.8,0.75), # 指定图例位置
                   legend.title = "",
                   legend.labs = c("High","Low"))
fig2d=mg_merge_plot(tcga.km$plot,tcga.km$table,ncol = 1,nrow = 2,heights = c(3.5,1),align = 'v')
fig2d

fig2e=plotMutiBar(table(tcga.risktype.cli$Status,tcga.risktype.cli$Risktype))
fig2e



#####GSE31210 #######
GSE31210_model_data=t(GSE31210_exp[intersect(pre.genes,rownames(GSE31210_exp)),GSE31210_cli$Samples])
colnames(GSE31210_model_data)=gsub('-','__',colnames(GSE31210_model_data))
GSE31210_model_data=merge(data.frame(Samples=GSE31210_cli$Samples,OS=GSE31210_cli$OS,OS.time=GSE31210_cli$OS.time),
                          data.frame(Samples=rownames(GSE31210_model_data),GSE31210_model_data),by='Samples')
rownames(GSE31210_model_data)=GSE31210_model_data$Samples
GSE31210_model_data=GSE31210_model_data[,-1]
GSE31210_model_data=crbind2DataFrame(GSE31210_model_data)
dim(GSE31210_model_data)

GSE31210.module.risk=get_riskscore(dat = GSE31210_model_data[GSE31210_cli$Samples,intersect(colnames(GSE31210_model_data),tcga.module.risk$module.gene)],
                                   #dat=GSE31210_model_data[GSE31210_cli$Samples,rownames(GSE31210.cox.fit)],
                                   os = GSE31210_model_data[GSE31210_cli$Samples,]$OS,
                                   os.time = GSE31210_model_data[GSE31210_cli$Samples,]$OS.time,
                                   step = F,direction = c("both", "backward", "forward")[1])
length(GSE31210.module.risk$module.gene)

fig2f=ggplotTimeROC(time = GSE31210.module.risk$result$time/365
                    ,status = GSE31210.module.risk$result$status
                    ,score = GSE31210.module.risk$result$riskscore
                    ,mks = c(1,3,5))
fig2f
GSE31210.risktype.cli=data.frame(GSE31210_cli,
                                 Riskscore=GSE31210.module.risk$result[GSE31210_cli$Samples,]$riskscorez)
GSE31210.risktype.cli$Risktype=ifelse(GSE31210.risktype.cli$Riskscore>0,'High','Low')

GSE31210.km=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ Risktype,
                                    data = GSE31210.risktype.cli),
                       data=GSE31210.risktype.cli,
                       conf.int = F,pval = T,fun = "pct",risk.table = T, size = 0.7,
                       title='GSE31210',ggtheme=custom_theme(),
                       linetype = c("solid", "dashed","strata")[1],
                       palette = risktype.col,
                       legend = c('top', 'bottom', 'left', 'right', 'none')[1],
                       #legend = c(0.8,0.75), # 指定图例位置
                       legend.title = "",
                       legend.labs = c("High","Low"))
fig2g=mg_merge_plot(GSE31210.km$plot,GSE31210.km$table,ncol = 1,nrow = 2,heights = c(3.5,1),align = 'v')
fig2g

fig2h=plotMutiBar(table(GSE31210.risktype.cli$Status,GSE31210.risktype.cli$Risktype))
fig2h


####GSE30219 #######
GSE30219_model_data=t(GSE30219_exp[intersect(pre.genes,rownames(GSE30219_exp)),GSE30219_cli$Samples])
colnames(GSE30219_model_data)=gsub('-','__',colnames(GSE30219_model_data))
GSE30219_model_data=merge(data.frame(Samples=GSE30219_cli$Samples,OS=GSE30219_cli$OS,OS.time=GSE30219_cli$OS.time),
                          data.frame(Samples=rownames(GSE30219_model_data),GSE30219_model_data),by='Samples')
rownames(GSE30219_model_data)=GSE30219_model_data$Samples
GSE30219_model_data=GSE30219_model_data[,-1]
GSE30219_model_data=crbind2DataFrame(GSE30219_model_data)
dim(GSE30219_model_data)

GSE30219.module.risk=get_riskscore(dat = GSE30219_model_data[GSE30219_cli$Samples,intersect(colnames(GSE30219_model_data),tcga.module.risk$module.gene)],
                                   #dat=GSE30219_model_data[GSE30219_cli$Samples,rownames(GSE30219.cox.fit)],
                                   os = GSE30219_model_data[GSE30219_cli$Samples,]$OS,
                                   os.time = GSE30219_model_data[GSE30219_cli$Samples,]$OS.time,
                                   step = F,direction = c("both", "backward", "forward")[1])
length(GSE30219.module.risk$module.gene)

fig2i=ggplotTimeROC(time = GSE30219.module.risk$result$time/365
              ,status = GSE30219.module.risk$result$status
              ,score = GSE30219.module.risk$result$riskscore
              ,mks = c(1,3,5))

GSE30219.risktype.cli=data.frame(GSE30219_cli,
                                 Riskscore=GSE30219.module.risk$result[GSE30219_cli$Samples,]$riskscorez)
GSE30219.risktype.cli$Risktype=ifelse(GSE30219.risktype.cli$Riskscore>0,'High','Low')

GSE30219.km=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ Risktype,
                                    data = GSE30219.risktype.cli),
                       data=GSE30219.risktype.cli,
                       conf.int = F,pval = T,fun = "pct",risk.table = T, size = 0.7,
                       title='GSE30219',ggtheme=custom_theme(),
                       linetype = c("solid", "dashed","strata")[1],
                       palette =risktype.col,
                       legend = c('top', 'bottom', 'left', 'right', 'none')[1],
                       #legend = c(0.8,0.75), # 指定图例位置
                       legend.title = "",
                       legend.labs = c("High","Low"))
fig2j=mg_merge_plot(GSE30219.km$plot,GSE30219.km$table,ncol = 1,nrow = 2,heights = c(3.5,1),align = 'v')
fig2k=plotMutiBar(table(GSE30219.risktype.cli$Status,GSE30219.risktype.cli$Risktype))

fig2=mg_merge_plot(mg_merge_plot(fig2a,fig2b,widths = c(2,1),labels = c('','C')),
                   mg_merge_plot(fig2c,fig2d,fig2e,fig2f,fig2g,fig2h,fig2i,fig2j,fig2k,
                                 ncol=3,nrow=3,labels = LETTERS[4:12]),
                   nrow=2,heights = c(1,3))
savePDF('results/01.AS_realated_module/Fig2.pdf',fig2,height = 20,width = 17)

library(tinyarray)
pdf('results/01.AS_realated_module/heatmap1.pdf',height = 4,width = 5.5)
draw_heatmap(t(tcga_model_data[order(tcga.risktype.cli$Riskscore),tcga.module.risk$module.gene]),
             group_list = factor(tcga.risktype.cli$Risktype[order(tcga.risktype.cli$Riskscore)]),
             cluster_cols = F,
             show_rownames = T,
             legend = T,
             annotation_legend = T,
             n_cutoff = 2,
             color_an =risktype.col)
dev.off()
pdf('results/01.AS_realated_module/heatmap2.pdf',height = 4,width = 5.5)
draw_heatmap(GSE31210_exp[GSE31210.module.risk$module.gene,order(GSE31210.risktype.cli$Riskscore)],
             group_list = factor(GSE31210.risktype.cli$Risktype[order(GSE31210.risktype.cli$Riskscore)]),
             cluster_cols = F,
             show_rownames = T,
             legend = T,
             annotation_legend = T,
             n_cutoff = 2,
             color_an = risktype.col)
dev.off()
pdf('results/01.AS_realated_module/heatmap3.pdf',height = 4,width = 5.5)
draw_heatmap(GSE30219_exp[GSE30219.module.risk$module.gene,order(GSE30219.risktype.cli$Riskscore)],
             group_list = factor(GSE30219.risktype.cli$Risktype[order(GSE30219.risktype.cli$Riskscore)]),
             cluster_cols = F,
             show_rownames = T,
             legend = T,
             annotation_legend = T,
             n_cutoff = 2,
             color_an = risktype.col)
dev.off()

######临床特征列线图#######
tcga_cox_datas=tcga.risktype.cli
colnames(tcga_cox_datas)
table(tcga_cox_datas$Age1)
table(tcga_cox_datas$T.stage)
tcga_cox_datas$T.stage[tcga_cox_datas$T.stage=='T1'|tcga_cox_datas$T.stage=='T2']<-'T1+T2'
tcga_cox_datas$T.stage[tcga_cox_datas$T.stage=='T3'|tcga_cox_datas$T.stage=='T4']<-'T3+T4'

table(tcga_cox_datas$N.stage)
tcga_cox_datas$N.stage[tcga_cox_datas$N.stage=='N1'|tcga_cox_datas$N.stage=='N2'|tcga_cox_datas$N.stage=='N3']<-'N1+N2+N3'

table(tcga_cox_datas$M.stage)

table(tcga_cox_datas$Stage)
tcga_cox_datas$Stage[tcga_cox_datas$Stage=='I'|tcga_cox_datas$Stage=='II']<-'I+II'
tcga_cox_datas$Stage[tcga_cox_datas$Stage=='III'|tcga_cox_datas$Stage=='IV']<-'III+IV'

#Age
tcga_cox_datas=crbind2DataFrame(tcga_cox_datas)
Age_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Age,
                             data=tcga_cox_datas))
Age_sig_cox_dat <- data.frame(Names=rownames(Age_sig_cox[[8]]),
                              HR = round(Age_sig_cox[[7]][,2],3),
                              lower.95 = round(Age_sig_cox[[8]][,3],3),
                              upper.95 = round(Age_sig_cox[[8]][,4],3),
                              p.value=round(Age_sig_cox[[7]][,5],3))
Age_sig_cox_dat

#Gender
Gender_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Gender,
                                data=tcga_cox_datas))
Gender_sig_cox_dat <- data.frame(Names=rownames(Gender_sig_cox[[8]]),
                                 HR = round(Gender_sig_cox[[7]][,2],3),
                                 lower.95 = round(Gender_sig_cox[[8]][,3],3),
                                 upper.95 = round(Gender_sig_cox[[8]][,4],3),
                                 p.value=round(Gender_sig_cox[[7]][,5],3))
Gender_sig_cox_dat

#T.stage
T.stage_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~T.stage,
                                 data=tcga_cox_datas))
T.stage_sig_cox_dat <- data.frame(Names=rownames(T.stage_sig_cox[[8]]),
                                  HR = round(T.stage_sig_cox[[7]][,2],3),
                                  lower.95 = round(T.stage_sig_cox[[8]][,3],3),
                                  upper.95 = round(T.stage_sig_cox[[8]][,4],3),
                                  p.value=round(T.stage_sig_cox[[7]][,5],3))
T.stage_sig_cox_dat

#N.stage
N.stage_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~N.stage,
                                 data=tcga_cox_datas))
N.stage_sig_cox_dat <- data.frame(Names=rownames(N.stage_sig_cox[[8]]),
                                  HR = round(N.stage_sig_cox[[7]][,2],3),
                                  lower.95 = round(N.stage_sig_cox[[8]][,3],3),
                                  upper.95 = round(N.stage_sig_cox[[8]][,4],3),
                                  p.value=round(N.stage_sig_cox[[7]][,5],3))
N.stage_sig_cox_dat

#M.stage
M.stage_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~M.stage,
                                 data=tcga_cox_datas))
M.stage_sig_cox_dat <- data.frame(Names=rownames(M.stage_sig_cox[[8]]),
                                  HR = round(M.stage_sig_cox[[7]][,2],3),
                                  lower.95 = round(M.stage_sig_cox[[8]][,3],3),
                                  upper.95 = round(M.stage_sig_cox[[8]][,4],3),
                                  p.value=round(M.stage_sig_cox[[7]][,5],3))
M.stage_sig_cox_dat

#Stage
Stage_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Stage,
                               data=tcga_cox_datas))
Stage_sig_cox_dat <- data.frame(Names=rownames(Stage_sig_cox[[8]]),
                                HR = round(Stage_sig_cox[[7]][,2],3),
                                lower.95 = round(Stage_sig_cox[[8]][,3],3),
                                upper.95 = round(Stage_sig_cox[[8]][,4],3),
                                p.value=round(Stage_sig_cox[[7]][,5],3))
Stage_sig_cox_dat


#riskscore
riskscore_sig_cox<-summary(coxph(formula=Surv(OS.time, OS)~Riskscore,
                                 data=tcga_cox_datas))
riskscore_sig_cox_dat <- data.frame(Names=rownames(riskscore_sig_cox[[8]]),
                                    HR = round(riskscore_sig_cox[[7]][,2],3),
                                    lower.95 = round(riskscore_sig_cox[[8]][,3],3),
                                    upper.95 = round(riskscore_sig_cox[[8]][,4],3),
                                    p.value=round(riskscore_sig_cox[[7]][,5],3))
riskscore_sig_cox_dat

sig_cox_dat <- rbind(Age_sig_cox_dat,
                     Gender_sig_cox_dat,
                     # T.stage_sig_cox_dat,
                     # N.stage_sig_cox_dat,
                     # M.stage_sig_cox_dat,
                     Stage_sig_cox_dat,
                     #Grade_sig_cox_dat,
                     riskscore_sig_cox_dat)
data.sig <- data.frame(Names=sig_cox_dat$Names,
                       p.value=sig_cox_dat$p.value,
                       sig_cox_dat$HR,
                       sig_cox_dat$lower.95,
                       sig_cox_dat$upper.95)
data.sig <- crbind2DataFrame(data.sig)
data.sig$Names
rownames(data.sig) <- c("Age",
                        "Gender",
                        # "T.stage",
                        # "N.stage",
                        # "M.stage",
                        "Stage",
                        #"Grade",
                        "RiskScore")
data.sig$Names <- rownames(data.sig)
data.sig
data.sig$p.value=ifelse(data.sig$p.value<0.001,'<0.001',data.sig$p.value)
pdf('results/01.AS_realated_module/Fig3a.pdf',height = 7,width = 7,onefile = F)
mg_forestplot_v2(data.sig,xlog = T,colgap =20,lineheight = 10,zero = 1,
                 boxsize = 0.2,lwd.zero=1,lwd.ci=1,lwd.xaxis=1
                 ,box_col='green',summary_col="black",lines_col='black',zero_col='grey'
                 ,xlab='Hazard Ratio',lty.ci = 5,graph.pos = 4)
dev.off()

#######多因素
#T.stage +  N.stage +M.stage
muti_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Age+Gender+Stage+Riskscore, 
                              data=tcga_cox_datas))
muti_cox_dat <- data.frame(Names=rownames(muti_sig_cox[[8]]),
                           HR = round(muti_sig_cox[[7]][,2],3),
                           lower.95 = round(muti_sig_cox[[8]][,3],3),
                           upper.95 = round(muti_sig_cox[[8]][,4],3),
                           p.value=round(muti_sig_cox[[7]][,5],3))
data.muti <- data.frame(Names=muti_cox_dat$Names,
                        p.value=muti_cox_dat$p.value,
                        muti_cox_dat$HR,
                        muti_cox_dat$lower.95,
                        muti_cox_dat$upper.95)
data.muti <- crbind2DataFrame(data.muti)
data.muti
rownames(data.muti) <- c("Age",
                         "Gender",
                         #"T.stage",
                         "Stage",
                         #"Grade",
                         "RiskScore")
data.muti$Names <- rownames(data.muti)
data.muti
data.muti$p.value=ifelse(data.muti$p.value<0.001,'<0.001',data.muti$p.value)
pdf('results/01.AS_realated_module/Fig3b.pdf',height =7,width = 7,onefile = F)
mg_forestplot_v2(data.muti,xlog = T,colgap =20,lineheight = 10,zero = 1,
                 boxsize = 0.2,lwd.zero=1,lwd.ci=1,lwd.xaxis=1
                 ,box_col='green',summary_col="black",lines_col='black',zero_col='grey'
                 ,xlab='Hazard Ratio',lty.ci = 5,graph.pos = 4)
dev.off()


######GSE30219临床特征
GSE30219_cox_datas=GSE30219.risktype.cli
colnames(GSE30219_cox_datas)

table(GSE30219_cox_datas$Stage)
GSE30219_cox_datas$Stage[GSE30219_cox_datas$Stage=='I'|GSE30219_cox_datas$Stage=='II']<-'I+II'
GSE30219_cox_datas$Stage[GSE30219_cox_datas$Stage=='III'|GSE30219_cox_datas$Stage=='IV']<-'III+IV'

#Age
GSE30219_cox_datas=crbind2DataFrame(GSE30219_cox_datas)
GSE30219_Age_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Age,
                                      data=GSE30219_cox_datas))
GSE30219_Age_sig_cox_dat <- data.frame(Names=rownames(GSE30219_Age_sig_cox[[8]]),
                                       HR = round(GSE30219_Age_sig_cox[[7]][,2],3),
                                       lower.95 = round(GSE30219_Age_sig_cox[[8]][,3],3),
                                       upper.95 = round(GSE30219_Age_sig_cox[[8]][,4],3),
                                       p.value=round(GSE30219_Age_sig_cox[[7]][,5],3))
GSE30219_Age_sig_cox_dat

#Gender
GSE30219_Gender_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Gender,
                                         data=GSE30219_cox_datas))
GSE30219_Gender_sig_cox_dat <- data.frame(Names=rownames(GSE30219_Gender_sig_cox[[8]]),
                                          HR = round(GSE30219_Gender_sig_cox[[7]][,2],3),
                                          lower.95 = round(GSE30219_Gender_sig_cox[[8]][,3],3),
                                          upper.95 = round(GSE30219_Gender_sig_cox[[8]][,4],3),
                                          p.value=round(GSE30219_Gender_sig_cox[[7]][,5],3))
GSE30219_Gender_sig_cox_dat



#Stage
GSE30219_Stage_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Stage,
                                        data=GSE30219_cox_datas))
GSE30219_Stage_sig_cox_dat <- data.frame(Names=rownames(GSE30219_Stage_sig_cox[[8]]),
                                         HR = round(GSE30219_Stage_sig_cox[[7]][,2],3),
                                         lower.95 = round(GSE30219_Stage_sig_cox[[8]][,3],3),
                                         upper.95 = round(GSE30219_Stage_sig_cox[[8]][,4],3),
                                         p.value=round(GSE30219_Stage_sig_cox[[7]][,5],3))
GSE30219_Stage_sig_cox_dat


#riskscore
GSE30219_riskscore_sig_cox<-summary(coxph(formula=Surv(OS.time, OS)~Riskscore,
                                          data=GSE30219_cox_datas))
GSE30219_riskscore_sig_cox_dat <- data.frame(Names=rownames(GSE30219_riskscore_sig_cox[[8]]),
                                             HR = round(GSE30219_riskscore_sig_cox[[7]][,2],3),
                                             lower.95 = round(GSE30219_riskscore_sig_cox[[8]][,3],3),
                                             upper.95 = round(GSE30219_riskscore_sig_cox[[8]][,4],3),
                                             p.value=round(GSE30219_riskscore_sig_cox[[7]][,5],3))
GSE30219_riskscore_sig_cox_dat

GSE30219_sig_cox_dat <- rbind(GSE30219_Age_sig_cox_dat,
                              GSE30219_Gender_sig_cox_dat,
                              GSE30219_Stage_sig_cox_dat,
                              GSE30219_riskscore_sig_cox_dat)
GSE30219_data.sig <- data.frame(Names=GSE30219_sig_cox_dat$Names,
                                p.value=GSE30219_sig_cox_dat$p.value,
                                GSE30219_sig_cox_dat$HR,
                                GSE30219_sig_cox_dat$lower.95,
                                GSE30219_sig_cox_dat$upper.95)
GSE30219_data.sig <- crbind2DataFrame(GSE30219_data.sig)
GSE30219_data.sig$Names
rownames(GSE30219_data.sig) <- c("Age",
                                 "Gender",
                                 "Stage",
                                 "RiskScore")
GSE30219_data.sig$Names <- rownames(GSE30219_data.sig)
GSE30219_data.sig
GSE30219_data.sig$p.value=ifelse(GSE30219_data.sig$p.value<0.001,'<0.001',GSE30219_data.sig$p.value)
pdf('results/01.AS_realated_module/Fig3c.pdf',height = 7,width = 7,onefile = F)
mg_forestplot_v2(GSE30219_data.sig,xlog = T,colgap =20,lineheight = 10,zero = 1,
                 boxsize = 0.2,lwd.zero=1,lwd.ci=1,lwd.xaxis=1
                 ,box_col='green',summary_col="black",lines_col='black',zero_col='grey'
                 ,xlab='Hazard Ratio',lty.ci = 5,graph.pos =4)
dev.off()




##多因素
GSE30219_muti_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Age+Gender+Stage+Riskscore, 
                                       data=GSE30219_cox_datas))
GSE30219_muti_cox_dat <- data.frame(Names=rownames(GSE30219_muti_sig_cox[[8]]),
                                    HR = round(GSE30219_muti_sig_cox[[7]][,2],3),
                                    lower.95 = round(GSE30219_muti_sig_cox[[8]][,3],3),
                                    upper.95 = round(GSE30219_muti_sig_cox[[8]][,4],3),
                                    p.value=round(GSE30219_muti_sig_cox[[7]][,5],3))
GSE30219.data.muti <- data.frame(Names=GSE30219_muti_cox_dat$Names,
                                 p.value=GSE30219_muti_cox_dat$p.value,
                                 GSE30219_muti_cox_dat$HR,
                                 GSE30219_muti_cox_dat$lower.95,
                                 GSE30219_muti_cox_dat$upper.95)
GSE30219.data.muti <- crbind2DataFrame(GSE30219.data.muti)
GSE30219.data.muti
rownames(GSE30219.data.muti) <- c("Age",
                                  "Gender",
                                  "Stage",
                                  "RiskScore")
GSE30219.data.muti$Names <- rownames(GSE30219.data.muti)
GSE30219.data.muti
GSE30219.data.muti$p.value=ifelse(GSE30219.data.muti$p.value<0.001,'<0.001',GSE30219.data.muti$p.value)
pdf('results/01.AS_realated_module/Fig3d.pdf',height =7,width = 7,onefile = F)
mg_forestplot_v2(GSE30219.data.muti,xlog = T,colgap =20,lineheight = 10,zero = 1,
                 boxsize = 0.2,lwd.zero=1,lwd.ci=1,lwd.xaxis=1
                 ,box_col='green',summary_col="black",lines_col='black',zero_col='grey'
                 ,xlab='Hazard Ratio',lty.ci = 5,graph.pos =4)
dev.off()



########列线图
pdf('results/01.AS_realated_module/Fig3e.pdf', width = 12, height = 10)
nom.plot=mg_nomogram(data.frame(RiskScore=tcga_cox_datas$Riskscore,
                                #Age=tcga_cox_datas$Age,
                                Stage=tcga_cox_datas$Stage),
                     os = tcga_cox_datas$OS.time,
                     status = tcga_cox_datas$OS,
                     mks = c(1,3,5)
)
dev.off()
mg_nomogram_buti(nom.plot$Mod,cut.time = c(1*365,3*365,5*365))

pdf('results/01.AS_realated_module/Fig3f.pdf', width = 12, height = 10)
nom.plot=mg_nomogram(data.frame(RiskScore=GSE30219_cox_datas$Riskscore,
                                #Age=tcga_cox_datas$Age,
                                Stage=GSE30219_cox_datas$Stage),
                     os = GSE30219_cox_datas$OS.time,
                     status = GSE30219_cox_datas$OS,
                     mks = c(1,3,5)
)
dev.off()
mg_nomogram_buti(nom.plot$Mod,cut.time = c(1*365,3*365,5*365))

