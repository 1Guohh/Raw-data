####04.潜在治疗策略【免疫治疗+耐药】#######
dir.create('results/06.immu_treatment')
##########IMV210队列#############
library("IMvigor210CoreBiologies")
data(cds)
pheno<-pData(cds)
head(pheno)
exper_tpm=mg_get_immu_pd1_treament_exp()
exper_id=exper_tpm$tpm
exper_id$symbol=rownames(exper_id)
rownames(exper_id)<-exper_id$symbol
exper_id$symbol<-NULL
range(exper_id)
exper_id_use<-log2(exper_id+1)
dim(exper_id_use)
# 31085   348
range(exper_id_use)
#rownames(exper_id_use)=gsub('-','__',rownames(exper_id_use))
exper_id_use[1:5,1:5]


IMvigor210_model_data=data.frame(OS=pheno$censOS,OS.time = pheno$os,
                                 t(exper_id_use[intersect(pre.genes,rownames(exper_id_use)),rownames(pheno)]))
head(IMvigor210_model_data)


imv210.module.risk=get_riskscore(dat = IMvigor210_model_data[,intersect(tcga.module.risk$module.gene,rownames(exper_id_use))],
                                 os = IMvigor210_model_data$OS,
                                 os.time = IMvigor210_model_data$OS.time,
                                 step = F,direction = c("both", "backward", "forward")[1])
length(imv210.module.risk$module.gene)
imv210.module.risk$model

imv210.module.risk$result$Risktype=ifelse(imv210.module.risk$result$riskscore>median(imv210.module.risk$result$riskscore),'High','Low')

ggplotTimeROC(time = imv210.module.risk$result$time/12
              ,status = imv210.module.risk$result$status
              ,score = imv210.module.risk$result$riskscore
              ,mks = c(1.5,1,2))

imv210_km=ggsurvplot(fit = survfit(Surv(time,status)~Risktype,
                                   data =imv210.module.risk$result),
                     data=imv210.module.risk$result,fun = "pct", 
                     palette = ggsci::pal_jama()(9)[1:2],
                     risk.table = F, title='IMvigor210',pval = T,
                     legend = c(0.8,0.75),legend.labs = c("High","Low"))
imv210_km
imv.risk<-cbind.data.frame(IMvigor210_model_data[rownames(pheno),c('OS','OS.time')],
                           Riskscore=imv210.module.risk$result[rownames(pheno),]$riskscore,
                           Risktype=imv210.module.risk$result[rownames(pheno),]$Risktype,
                           data.frame(binaryResponse=pheno$binaryResponse,
                                      Response=pheno$`Best Confirmed Overall Response`,
                                      IC=pheno$`IC Level`,
                                      TC=pheno$`TC Level`,
                                      IP=pheno$`Immune phenotype`,
                                      Stage=pheno$`TCGA Subtype`))
head(imv.risk)


imv.risk1=imv.risk
imv.risk1=crbind2DataFrame(imv.risk1)
table(imv.risk1$Stage)
imv.risk1$Stage[imv.risk1$Stage=='I'|imv.risk1$Stage=='II']='I+II'
imv.risk1$Stage[imv.risk1$Stage=='III'|imv.risk1$Stage=='IV']='III+IV'
#stage早期和晚期KM
imv210_km1=ggsurvplot(fit=survfit( Surv(OS.time/12, OS) ~ Risktype,
                                   data = imv.risk1[which(imv.risk1$Stage=='I+II'),]),
                      data=imv.risk1[which(imv.risk1$Stage=='I+II'),],
                      #surv.median.line = "hv",
                      conf.int = F,pval = T,fun = "pct",risk.table = F, 
                      title='IMvigor210 Stage I+II',
                      linetype = c("solid", "dashed","strata")[1],
                      #legend = c('top', 'bottom', 'left', 'right', 'none')[5],
                      legend = c(0.8,0.75), # 指定图例位置
                      legend.title = "",
                      legend.labs = c("High","Low"))
imv210_km1

imv210_km2=ggsurvplot(fit=survfit( Surv(OS.time/12, OS) ~ Risktype,
                                   data = imv.risk1[which(imv.risk1$Stage=='III+IV'),]),
                      data=imv.risk1[which(imv.risk1$Stage=='III+IV'),],
                      conf.int = F,pval = T,fun = "pct",risk.table = F,
                      title='IMvigor210 Stage III+IV',
                      linetype = c("solid", "dashed","strata")[1],
                      #legend = c('top', 'bottom', 'left', 'right', 'none')[5],
                      legend = c(0.8,0.75), # 指定图例位置
                      legend.title = "",
                      legend.labs = c("High","Low"))
imv210_km2

table(imv.risk$binaryResponse)
imv210_boxplot=imv.risk[which(imv.risk$binaryResponse!='NA'),] %>%
  ggplot(aes(x=binaryResponse, y=Riskscore,fill = binaryResponse)) +
  #geom_violin()+  
  scale_fill_manual(values = pal_nejm()(10)[3:4])+
  geom_boxplot()+
  theme_classic(base_size = 20)+
  ggpubr::stat_compare_means(aes(group=binaryResponse), label = "p.format", method = 'wilcox.test')+
  theme_classic()+
  theme(legend.position = 'none',axis.text = element_text(color = 'black'),
        axis.text.x = element_text(size = 12),axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 12),axis.title.y = element_text(size = 15))
imv210_boxplot


#柱状图
paste0('-log10(p.value)=',round(-log10(chisq.test(table(imv.risk$binaryResponse,imv.risk$Risktype))$p.value),2))
compaired <- list(c("High", "Low"))
results=prop.table(table(imv.risk$binaryResponse,imv.risk$Risktype),margin=2)
results1=reshape2::melt(results)
colnames(results1)<-c("binaryResponse","Senescore","Percentage")
results1$Percentage<-round(results1$Percentage,digits=2)

imv210_bar=ggplot(results1,aes(x=Senescore,y=Percentage,fill=binaryResponse))+
  geom_bar(position = "fill",stat="identity")+
  scale_fill_manual(values = pal_nejm()(10)[3:4])+
  #geom_signif(comparisons = compaired, map_signif_level = T,test = 'chisq.test')+
  theme_bw()+labs(x="Risktype", y = "Percentage",
                  title = paste0('P value=',round(chisq.test(table(imv.risk$binaryResponse,imv.risk$Risktype))$p.value,5)))+
  geom_text(aes(label = Percentage),position=position_stack(vjust =0.5),size = 5)+
  theme(legend.position = 'top')
imv210_bar

imv210_fig=mg_merge_plot(imv210_km$plot,imv210_boxplot,imv210_bar,ncol=3,labels = LETTERS[1:3])
imv210_fig


#####TIDE+IC50########
# tcga_tide_dat <- t(scale(t(tcga_tpm_log_T[,luad_as_cli$Samples]),scale = F))
# dim(tcga_tide_dat)
# write.table(tcga_tide_dat,file = 'results/tcga_tide_dat.txt',quote = F, sep = '\t')

tcga_tide_res<-read.csv('results/tcga.tide.res.csv',row.names = 1,stringsAsFactors = F)
head(tcga_tide_res)
tcga_tide_cli=cbind(tcga.risktype.cli,tcga_tide_res[tcga.risktype.cli$Samples,])
tide_sel=c('TIDE','IFNG','Exclusion','Dysfunction','MDSC')
#pdf('results/04.subtype.immmu/Fig4d.pdf',height = 6,width = 6,onefile = F)
tide.p1=my_mutiboxplot(dat = tcga_tide_res[tcga.risktype.cli$Samples,tide_sel[1:4]],
                       group = tcga.risktype.cli$Risktype,fill = 'Risktype',
                       group_cols = risktype.col,
                       ylab = 'Score',angle = 0,hjust = 0.5)
tide.p1

tide.results=prop.table(table(tcga_tide_cli$Responder,tcga_tide_cli$Risktype),margin=2)
tide.results=reshape2::melt(tide.results)
colnames(tide.results)<-c("Responder","Risktype","Percentage")
tide.results$Percentage<-round(tide.results$Percentage,digits=2)

chisq.test(table(tcga_tide_cli$Responder,tcga_tide_cli$Risktype))
tide.p2=ggplot(tide.results,aes(x=Risktype,y=Percentage,fill=Responder))+
  geom_bar(position = "fill",stat="identity")+
  scale_fill_manual(values = pal_nejm()(10)[6:7])+
  #geom_signif(comparisons = compaired, map_signif_level = T,test = 'chisq.test')+
  theme_bw()+labs(x="Risktype", y = "Percentage",
                  title = paste0('P value=1.267e-12'))+
  geom_text(aes(label = Percentage),position=position_stack(vjust =0.5),size = 5)+
  theme(legend.position = 'top')
tide.p2

tide_fig=mg_merge_plot(tide.p1,tide.p2,widths = c(2,1),labels = LETTERS[4:5])


fig8=mg_merge_plot(imv210_fig,tide_fig,nrow=2)
savePDF('results/08.immu_treatment/Fig8.pdf',fig8,height = 12,width = 16)

########IC50 药物敏感性
library(pRRophetic)
#library(ggplot2)
############### Cisplatin,顺铂
# set.seed(12345)
# predictedPtype_Cisplatin <- pRRopheticPredict(as.matrix(tcga_tpm_log_T)
#                                               , "Cisplatin"
#                                               , selection=1
#                                               ,dataset = "cgp2016")
# predictedPtype_Cisplatin <- data.frame(predictedPtype_Cisplatin)
# 
# tcga_durg_ic50_res <- predictedPtype_Cisplatin

drugs <- c("Cisplatin","Erlotinib","Rapamycin","Sunitinib","PHA-665752","MG-132","Paclitaxel","Cyclopamine","AZ628","Sorafenib","VX-680","Imatinib","TAE684","Crizotinib","Saracatinib","S-Trityl-L-cysteine","Z-LLNle-CHO","Dasatinib","GNF-2","CGP-60474","CGP-082996","A-770041","WH-4-023","WZ-1-84","BI-2536","BMS-509744","CMK","Pyrimethamine","JW-7-52-1","A-443654","GW843682X","MS-275","Parthenolide","KIN001-135","TGX221","Bortezomib","XMD8-85","Roscovitine","Salubrinal","Lapatinib","Vinorelbine","NSC-87877","QS11","CP466722","Midostaurin","Shikonin","AKT inhibitor VIII","Embelin","Bexarotene","Bleomycin","Phenformin")
length(drugs)
# for (drug in drugs) {
#   print(drug)
#   set.seed(12345)
#   tmpic50 <- pRRopheticPredict(as.matrix(tcga_tpm_log_T)
#                                , drug
#                                , selection=1
#                                , dataset = "cgp2016")
#   tmpic50 <- data.frame(tmpic50)
#   colnames(tmpic50) <- drug
#   tcga_durg_ic50_res <- cbind(tcga_durg_ic50_res, tmpic50)
# }
#save(tcga_durg_ic50_res,file='results/tcga_durg_ic50.RData')
load('results/tcga_durg_ic50.RData')
head(tcga_durg_ic50_res)
dim(tcga_durg_ic50_res)
# tcga_durg_ic50_res=tcga_durg_ic50_res[,-1]


tcga.ic50.dat=cbind.data.frame(Riskscore=tcga.risktype.cli$Riskscore,
                               tcga_durg_ic50_res[tcga.risktype.cli$Samples,])

ic50.cor.RS=Hmisc::rcorr(as.matrix(tcga.ic50.dat),type = 'spearman')
ic50.cor.RS.res=data.frame(Names=names(ic50.cor.RS$r['Riskscore',]),
                           cor=as.numeric(ic50.cor.RS$r['Riskscore',]),
                           p.val=as.numeric(ic50.cor.RS$P['Riskscore',]))
ic50.cor.RS.res=ic50.cor.RS.res[-1,]
head(ic50.cor.RS.res)
colnames(ic50.cor.RS.res)=c('Drugs','cor','pvalue')
rownames(ic50.cor.RS.res)=ic50.cor.RS.res$Drugs
ic50.cor.RS.res$Drugs=factor(ic50.cor.RS.res$Drugs,
                             levels = ic50.cor.RS.res$Drugs[order(ic50.cor.RS.res$cor,decreasing = T)], ordered=TRUE)
ic50.cor.RS.res$pvalue=ifelse(ic50.cor.RS.res$pvalue==0,1e-16,ic50.cor.RS.res$pvalue)
head(ic50.cor.RS.res)

high.suit.drug=as.character(ic50.cor.RS.res$Drugs[which(ic50.cor.RS.res$pvalue<0.05 & ic50.cor.RS.res$cor < -0.5)])
length(high.suit.drug)
low.suit.drug=as.character(ic50.cor.RS.res$Drugs[which(ic50.cor.RS.res$pvalue<0.05 & ic50.cor.RS.res$cor  > 0.3)])
length(low.suit.drug)

#pdf('results/10.risktype.immu/Fig10d.pdf',height = 6,width = 14,onefile = F)

fig9a=ic50.cor.RS.res[high.suit.drug,] %>%
  ggplot(aes(x=Drugs, y=cor)) +
  geom_point(aes(fill=-log10(pvalue)),size = 3, pch = 21) +
  geom_segment(aes(x=Drugs, xend=Drugs,y=0 ,yend=cor),
               size=0.5,linetype=2)+ylim(0,-0.75)+coord_flip()

fig9b=my_mutiboxplot(dat = tcga_durg_ic50_res[tcga.risktype.cli$Samples,high.suit.drug],
                     group = tcga.risktype.cli$Risktype,#legend.position = 'top',
                     group_cols =risktype.col,ylab = 'IC50')+
  theme_classic2()+theme(legend.position = 'top')


fig9c=ic50.cor.RS.res[low.suit.drug,] %>%
  ggplot(aes(x=Drugs, y=cor)) +
  geom_point(aes(fill=-log10(pvalue)),size = 3, pch = 21) +
  geom_segment(aes(x=Drugs, xend=Drugs,y=0 ,yend=cor),
               size=0.5,linetype=2)+coord_flip()
fig9d=my_mutiboxplot(dat = tcga_durg_ic50_res[tcga.risktype.cli$Samples,low.suit.drug],
                     group = tcga.risktype.cli$Risktype,#legend.position = 'top',
                     group_cols = risktype.col,ylab = 'IC50')+
  theme_classic2()+theme(legend.position = 'top')




