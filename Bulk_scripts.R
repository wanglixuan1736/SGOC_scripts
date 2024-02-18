##两个数据集中表达相反的114代谢通路
##热图
#加载数据
#加载这些通路
library(dplyr)
library(tidyverse)
library(pheatmap)
path<-read.csv("4/only_CA/sc_TCGA_DIFF_PATHWAY.csv")
deg_path<-read.csv("4/only_CA/DEA.ssgsea_NC2018_HNSC_T_N_ssgsea.csv")

##加载sc的通路对比数据

pathways<-path$Pathway
deg_path1<-deg_path %>%
  filter(Pathways %in% pathways)
deg_path1<-deg_path1[order(deg_path1$t,decreasing = T),]
pa<-deg_path1$Pathways
#加载tcga的代谢通路表达矩阵
load("4/only_CA/gsva_NC2018_HNSC_N_T_ssgsea.Rdata")
#加载sc的代谢通路表达矩阵
load("4/noLP新分析/epma_gsvameta_114_cbind_metadata.Rdata")
##构建热图需要的数据
epma_gsvameta_114<-epma_gsvameta_114[order(epma_gsvameta_114$cell_type_me),]
annotation_col <- as.data.frame(epma_gsvameta_114[,1])
rownames(annotation_col)<-rownames(epma_gsvameta_114)
colnames(annotation_col)<-"group"
epma_gsvameta_114_1<-epma_gsvameta_114[,-c(1:3)]
epma_gsvameta_114_1<-as.data.frame(t(epma_gsvameta_114_1))
epma_gsvameta_114_1<-epma_gsvameta_114_1 %>%
  rownames_to_column("Pathways")
pa_1<-as.data.frame(pa)
pa_1$order<-c(1:40)
colnames(pa_1)[1]<-"Pathways"
test<-epma_gsvameta_114_1[1:10,1:10]
epma_gsvameta_114_1<-merge(pa_1,epma_gsvameta_114_1,by="Pathways")
epma_gsvameta_114_1<-epma_gsvameta_114_1[order(epma_gsvameta_114_1$order),]
test<-epma_gsvameta_114_1[1:10,1:10]
rownames(epma_gsvameta_114_1)<-epma_gsvameta_114_1[,1]
epma_gsvameta_114_1<-epma_gsvameta_114_1[,-c(1,2)]
hetadata<-epma_gsvameta_114_1
##加载分组信息
load("clinical/HNSC_exp_clin.Rdata")
clin<-HNSC_exp_clin[,c(1,15)]

### 画热图
gsva_NC2018<-as.data.frame(gsva_NC2018)
gsva_NC2018<-subset(gsva_NC2018,rownames(gsva_NC2018) %in% pa)
gsva_NC2018<-gsva_NC2018 %>%
  rownames_to_column("Pathways")
gsva_NC2018<-merge(deg_path1,gsva_NC2018,by="Pathways")
gsva_NC2018<-gsva_NC2018[order(gsva_NC2018$t,decreasing = T),]
rownames(gsva_NC2018)<-gsva_NC2018[,1]
gsva_NC2018<-gsva_NC2018[,-c(1:8)]
gsva_NC2018<-as.data.frame(t(gsva_NC2018))
gsva_NC2018<-gsva_NC2018 %>%
  rownames_to_column("TCGA_id") 
gsva_NC2018<-merge(clin,gsva_NC2018,by="TCGA_id")
gsva_NC2018<-gsva_NC2018[order(gsva_NC2018$sample),]

rownames(gsva_NC2018)<-gsva_NC2018[,1]
gsva_NC20181<-gsva_NC2018[,-c(1,2)]
hetadata<-as.data.frame(t(gsva_NC20181))
### 创建分组信息
annotation_col <- gsva_NC2018[,c(1,2)]

annotation_col2<-as.data.frame(annotation_col[,-1])
rownames(annotation_col2)<-annotation_col[,1]
colnames(annotation_col2)[1]<-"group"
hetadata<-as.data.frame(hetadata)

###
pheatmap(hetadata)
###
library(pheatmap)
ann_colors = list( group = c(Epithelial_cells="#24708B", Malignant_cells="#CDC673"))
pheatmap(hetadata, #热图的数据
         cluster_rows = F,#行聚类
         cluster_cols = F,#列聚类，可以看出样本之间的区分度
         annotation_col =annotation_col, #标注样本分类
         annotation_colors = ann_colors,# 列表指定注释行和注释列的颜色
         annotation_legend=TRUE, # 显示注释
         show_rownames = T,# 显示行名
         show_colnames = F,# 显示行名
         scale = "row", #以行来标准化
         clustering_method = "ward.D2",
         color =colorRampPalette(colors = c("#191970","white","#EEB422"))(100),#调色
         #filename = "heatmap_F.pdf",#是否保存
         cellwidth = 0.05, cellheight = 8,# 格子比例
         fontsize =6)


load(file = "clinical/HNSC_exp_clin.Rdata")
exprSet1<-HNSC_exp_clin
test1<-exprSet1[1:20,1:20]
exprSet1<-exprSet1[-153,]
rownames(exprSet1)<-exprSet1[,1]
exprSet1<-exprSet1[,-1]

library(tibble)
library(dplyr)
library(tidyr)
gene_p53<-exprSet1[,c( "sample","TP53")]
load(file = "output/metadata_tumor_riskgroup_0.4_bestcutpoint.Rdata")
meta<-metadata
meta$ID<-substr(metadata$TCGA_id,1,15)
gene_p53<-rownames_to_column(gene_p53)
colnames(gene_p53)[1]<-"ID"
gene_p53<-gene_p53[,-2]
data<-merge(gene_p53,meta,by="ID")
rownames(data)<-data[,1]
data<-data[,-c(1,3)]
data<-rownames_to_column(data)
####所以在此处添加上了sample
genelist<-exprSet1[,c( "sample","SERINC1","PLOD2","HPRT1","TBPL1","SLC44A4",
                       "GIMAP7","MTHFS","PKLR","MTHFR")]

genelist<-rownames_to_column(genelist)
gene_TRIM17<-rownames_to_column(gene_TRIM17)
colnames(gene_TRIM17)[1]<-"ID"
gene_TRIM17$TCGA_id<-substr(gene_TRIM17$ID,1,15)
gene_TRIM17$sample <- ifelse(substring(gene_TRIM17$TCGA_id,14,15)=="01","Tumor","Normal")
data<-rownames_to_column(gene_TRIM17)
#################################################################
####################################################################
###韦恩图画图
rm(list=ls())
#install.packages("VennDiagram") ##下载并加载包
library(VennDiagram)
load(file = "R_data/0.1/HNSC_TvsN_resLFC_protein_gene.Rdata")
genelist<-read.csv(file = "resource/oneC_genelist.csv")
res<-res_1 %>% 
  filter(adj.P.Val<0.05)
A<-genelist$oneC
B<-res$gene
venn.plot <- venn.diagram(
  list(NAM=A,differgene=B),
  filename = "NAM_DEGsvenn.tiff",##韦恩图的名字
  lty = 1,
  lwd = 1,
  col = "black",  ##圈的颜色
  fill = c("red",  "blue"),##对应每个圈的颜色，有几个数据集，就需要有相应数量的颜色
  alpha = 0.60,
  cat.col = "black",##此处设置每个数据集的名称颜色，也可以使用c（）函数输入三种颜色
  cat.cex = 0.8,
  cat.fontface = "bold",
  margin = 0.07,
  cex = 0.8
)



#以2个分组为例
#指定统计的分组列，并设置作图颜色、字体样式等
venn_list <- list(NAM=A,differgene=B)

venn.diagram(venn_list, filename = 'venn2.png', imagetype = 'png', 
             fill = c('red', 'blue'), alpha = 0.50, cat.col = rep('black', 2), 
             col = 'black', cex = 1.5, fontfamily = 'serif', 
             cat.cex = 1.5, cat.fontfamily = 'serif')


###单因素COX
rm(list = ls())
top9<-as.data.frame(data.table::fread("output/0.1/0.01.txt"))
rownames(top9)<-top9[,1]
top9<-top9[,-1]
HR=gsub("[\\(\\)]","",top9$`HR (95% CI for HR)`)
HR=gsub("-"," ",HR)
HR=as.data.frame(do.call(cbind,strsplit(HR," ")),stringsAsFactors=F)
names(HR)=rownames(top9)
#################################
#开始绘图，直接保存到pdf文件中,第一种绘图方式
#################################
pdf(file="output/5/univariate_forest_0.01.pdf",width=7)
#左边和右边边距稍微留多一点来写变量名称，pvalue和HR
par(mar=c(5,6,4,13))
#先用小方块画出HR
plot(as.numeric(HR[1,]),1:dim(HR)[2],
     pch=15,cex=2,col="blue",bty='n',yaxt='n',ylab=NA,xlab="Hazard Ratio",
     xlim=range(as.numeric(unlist(HR)))
)
#添加中线
abline(v=1,col="grey",lwd=2,lty=2)
for(i in 1:ncol(HR)){
  x=as.numeric(HR[2:3,i])
  #循环画出CI
  lines(x,c(i,i),col="blue")
  #添加变量名
  text(0.2,i,rownames(top9)[i],xpd=T,adj = c(0,0))
  #添加p值
  text(2.1,i,as.numeric(top9[i,1]),xpd=T,adj = c(0,0))
  #添加HR和CI
  text(2.7,i,as.character(top9[i,2]),xpd=T,adj = c(0,0))
}
#添加标题
text(2.1,ncol(HR)+0.5,"pvalue",xpd=T,adj = c(0,0))
text(2.7,ncol(HR)+0.5,"HR(CI)",xpd=T,adj = c(0,0))
library(export)
dev.off()

#################################################################
###cofficients画图
coff<-coff[,1:2]
colnames(coff)<-c("gene","cofficient")
library(ggplot2)
library(hrbrthemes)
library(data.table)
coff$col<-ifelse(coff$cofficient>0,1,2)
coff$col<-as.factor(coff$col)
ggplot(coff, aes(x = cofficient, y = gene,fill=col)) + 
  geom_bar(stat = "identity",   
           show.legend = F,   
           width = .5) + 
  xlab("cofficient") + 
  ylab("gene") +  
  theme(panel.background = element_blank(),
        axis.line.x = element_line(colour = "black")) # 去掉背景格子


###train生存分组
rm(list = ls())
#加载这两个R包
library("survival")
library("survminer")
library(dplyr)
library(tidyverse)
library(tidyr)
load(file = "clinical/HNSC_exp_clin_train_set.Rdata")
test<-train_set[1:10,1:20]

exprSet1<-as.data.frame(t(train_set))
test2<-exprSet1[1:10,1:20]
data.clin<-train_set[,1:14]
data.clin$times<-data.clin$OS.time/365
colnames(train_set)[6]<-"status"

aaa2<-as.data.frame(data.table::fread("4/6/lasso_cox_folum_4GENE.txt"))
rownames(aaa2)<-aaa2[,1]
aaa2<-aaa2[,-1]

tcgaexp <- exprSet1[c(which(rownames(exprSet1)%in%rownames(aaa2))),]
tcgaexp <- tcgaexp[order(rownames(tcgaexp)),]
aaa2<- aaa2[order(rownames(aaa2)),]
###########################################3
#将tcgaexp中数字由“chr”变成“num”
tcgaexp<-as.data.frame(lapply(tcgaexp, as.numeric))
rownames(tcgaexp)<-c(rownames(aaa2))
colnames(tcgaexp)<-c(colnames(exprSet1))
tcgascore <- aaa2$coef*tcgaexp
tcgascore <- t(colSums(tcgascore))
tcgascore <- t(tcgascore)
colnames(tcgascore) <- c('score')

data.clin2 <- cbind(data.clin,tcgascore)
data.clin2<-rownames_to_column(data.clin2)
data.sur <- data.clin2[,c(1,17,7,16)]
names(data.sur) <- c('ID','score','status','time')
b <- surv_cutpoint(data.sur,time = 'time',event = 'status',variables = 'score',minprop = 0.5)
cutpoint <- summary(b)$cutpoint
high <- length(which(data.sur$score > cutpoint))
data.sur2 <- data.sur[order(data.sur$score,decreasing = T),]
rownames(data.sur2) <- data.sur2$ID
data.sur2 <- cbind(data.sur2,group = rep(c('High','Low'),c(high,250-high)))

#######################################################################
####test
rm(list = ls())
#加载这两个R包
library("survival")
library("survminer")
library(dplyr)
library(tidyverse)
library(tidyr)
#加载TCGA临床数据
load("clinical/HNSC_exp_clin_test_set.Rdata")
exprSet<-test_set
exprSet$time<-exprSet$OS.time/365
colnames(exprSet)[6]<-"status"


aaa2<-as.data.frame(data.table::fread("4/5/lasso_cox_folum_top15_4GENE.txt"))
rownames(aaa2)<-aaa2[,1]
aaa2<-aaa2[,-1]
exprSet1<-as.data.frame(t(exprSet))

tcgaexp <- exprSet1[c(which(rownames(exprSet1)%in%rownames(aaa2))),]
tcgaexp <- tcgaexp[order(rownames(tcgaexp)),]
aaa2<- aaa2[order(rownames(aaa2)),]
###########################################3
#将tcgaexp中数字由“chr”变成“num”
tcgaexp<-as.data.frame(lapply(tcgaexp, as.numeric))
rownames(tcgaexp)<-c(rownames(aaa2))
colnames(tcgaexp)<-c(colnames(exprSet1))
tcgascore <- aaa2$coef*tcgaexp
tcgascore <- t(colSums(tcgascore))
tcgascore <- t(tcgascore)
colnames(tcgascore) <- c('score')

data.clin<-exprSet[,1:14]
data.clin$times<-data.clin$OS.time/365

data.clin2 <- cbind(data.clin,tcgascore)
data.clin2<-rownames_to_column(data.clin2)
data.sur <- data.clin2[,c(1,17,7,16)]
names(data.sur) <- c('ID','score','status','time')
data.sur$status<-as.logical(data.sur$status)
#b <- surv_cutpoint(data.sur,time = 'time',event = 'status',variables = 'score',minprop = 0.4)
#cutpoint <- summary(b)$cutpoint
cutpoint<-8.969836
high <- length(which(data.sur$score > cutpoint))
data.sur2 <- data.sur[order(data.sur$score,decreasing = T),]
rownames(data.sur2) <- data.sur2$ID
data.sur2 <- cbind(data.sur2,group = rep(c('High','Low'),c(high,249-high)))

###################################################################
#total_TCGA数据验证
rm(list = ls())
setwd("E://0/1/2/3/TCGA/")
#加载这两个R包
library("survival")
library("survminer")
library(dplyr)
library(tidyverse)
library(tidyr)
#加载TCGA临床数据
load("clinical/HNSC_exp_clin_tumor_500.Rdata")
test1<-exprSet[1:20,1:20]
rownames(exprSet)<-exprSet[,1]
exprSet<-exprSet[,-1]
exprSet$time<-exprSet$OS.time/365
exprSet1<-as.data.frame(t(exprSet))
test2<-exprSet1[1:10,1:20]
data.clin<-exprSet[,1:14]
data.clin$times<-data.clin$OS.time/365
colnames(exprSet)[6]<-"status"

aaa2<-as.data.frame(data.table::fread("4/5/lasso_cox_folum_top15_4GENE.txt"))
rownames(aaa2)<-aaa2[,1]
aaa2_1<-as.data.frame(aaa2[,-1])
rownames(aaa2_1)<-rownames(aaa2)
tcgaexp <- exprSet1[c(which(rownames(exprSet1)%in%rownames(aaa2))),]
tcgaexp <- tcgaexp[order(rownames(tcgaexp)),]
aaa2<- aaa2[order(rownames(aaa2)),]
###########################################3
#将tcgaexp中数字由“chr”变成“num”
tcgaexp<-as.data.frame(lapply(tcgaexp, as.numeric))
rownames(tcgaexp)<-c(rownames(aaa2))
colnames(tcgaexp)<-c(colnames(exprSet1))
tcgascore <- aaa2$coef*tcgaexp
tcgascore <- t(colSums(tcgascore))
tcgascore <- t(tcgascore)
colnames(tcgascore) <- c('score')



data.clin2 <- cbind(data.clin,tcgascore)
data.clin2<-rownames_to_column(data.clin2)
data.sur <- data.clin2[,c(1,17,7,16)]
names(data.sur) <- c('ID','score','status','time')
b <- surv_cutpoint(data.sur,time = 'time',event = 'status',variables = 'score',minprop = 0.3)
cutpoint <- summary(b)$cutpoint
cutpoint<-8.96836
high <- length(which(data.sur$score > cutpoint))
data.sur2 <- data.sur[order(data.sur$score,decreasing = T),]
rownames(data.sur2) <- data.sur2$ID
data.sur2 <- cbind(data.sur2,group = rep(c('High','Low'),c(high,499-high)))

ibrary(dplyr)
library(survival) #核心分析

### 准备数据：去掉小于30天的,将天数变为年;去掉正常组织
data.sur2$OS.time<-data.sur2$time*365
coxdata <- data.sur2 %>% 
  filter(OS.time >= 30) 

### Surv()：创建一个生存对象，时间/事件
rt <- data.frame(coxdata[,c(4,3)],riskScore=coxdata[,5])
mySurv <- Surv(rt$time, rt$status) 
group<-rt$riskScore
### 使用survfit函数拟合一条生存曲线
### 使用summary函数查看模型汇总结果
sfit <- survfit(formula = mySurv~group, data = rt)
summary<-summary(sfit)

### 画Kaplan-Meier曲线
#   ggsurvplot_list() 绘制多个对象
#   ggsurvplot_facet() 分面到多个panels
#   ggsurvplot_group_by() 一幅图中多个分组
#   ggsurvplot_add_all() 总合所有的情况
#   ggsurvplot_combine() 一个图中结合多个survfit对象

library(survminer) 
ggsurvplot(sfit, 
           conf.int=TRUE, pval=TRUE, risk.table=TRUE, 
           #legend.labs=c("H1", "H2","L1","L2","M"), 
           #legend.title="1C",  
           #palette=c("dodgerblue2", "orchid2","green", "red","black"), 
           title="Kaplan-Meier Curve", 
           risk.table.height=.3)
library(export)
ggsurvplot(sfit, 
           size=1.25,#曲线粗细
           legend = c(0.8,0.8), # 指定图例strata位置
           #surv.plot.height=0.25,# 生存图的高度，默认为0.75；# 当risk.table = FALSE时忽略
           ggtheme = theme_classic(base_size = 20), #想要网格就运行这行
           conf.int = F, #不画置信区间，想画置信区间就把F改成T
           #conf.int.style = "step",#置信区间的类型，还可改为ribbon
           censor = F, #不显示观察值所在的位置
           #surv.median.line = "hv", #添加中位生存曲线
           linetype = "strata",# 改变曲线的类型
           censor.shape="+", censor.size = 4,#更改删失点形状及大小，默认为"+", 可选"|"
           #palette = c("#E7B800", "#2E9FDF"), #线的颜色对应高、低
           ylab="Overall survival",xlab = " Time (years)", # xlab, ylab 分别指x轴和y轴标签更改横纵坐标
           xlim = c(0,10),# xlim, ylim # 指定x轴和y轴的范围，如xlim = c(0,30), ylim = c(0,1) 
           break.x.by = 5,#横坐标间隔
           #risk.table = TRUE,tables.height = 0.32,#加risk table
           font.legend = c(14, "bold","black"),#图例的字体大小
           font.title = c(14,"bold"),font.x = c(14,"bold"),font.y = c(14,"bold"),#设置其他字体大小
           font.xtickslab = c(14, "bold","black"),font.ytickslab = c(14, "bold","black"),#坐标轴刻度线字体
           #在左下角标出pvalue、HR、95% CI
           #太小的p value标为p < 0.001
           pval = TRUE)

### 整体比较
survdiff(mySurv~group, data = rt)

rm(list=ls())
#BiocManager::install("GSVA")
library("GSVA")
library("tidyverse")
gene_sets <- as.matrix(t(data.table::fread("E:/0/1/7.1/7/xinfenxi/malignant_markergene_logfc2.5.csv")))
#  Generate a list that contains genes in genesets
#gene_sets<-as.data.frame(data.table::fread("resource/NUM_genelist.csv"))
#gene_sets<-as.character(gene_sets)
#rownames(gene_sets)<-"fe"
gs <- list()
for (i in 1:nrow(gene_sets)){
  a <- as.vector(gene_sets[i,1:ncol(gene_sets)])
  a <- na.omit(a)
  a <- a[a != ""]
  a <- matrix(a, ncol = 1)
  gs[[length(gs)+1]] <- a
  rm(a,i)
}
names(gs) <- rownames(gene_sets)
load(file = "clinical/HNSC_exp_clin.Rdata")
exprSet1<-HNSC_exp_clin
test1<-exprSet1[1:20,1:20]
exprSet1<-exprSet1[-153,]
rownames(exprSet1)<-exprSet1[,1]
exprSet1<-exprSet1[,-1]
expr_T<-exprSet1[exprSet1$sample=="Tumor",]
exprSet2<-expr_T[,-c(1:14)]
class(exprSet2)

test<-exprSet2[1:20,1:20]
exprSet2<-as.data.frame(t(exprSet2))

HNSC_ssgsva <- gsva(as.matrix(exprSet2),gs)
HNSC_ssgsva <- as.data.frame(t(HNSC_ssgsva))
HNSC_ssgsva <- rownames_to_column(HNSC_ssgsva,var = 'Id') 

library("survival")
library("survminer")
library(dplyr)
library(tidyverse)
library(tidyr)
clin<-exprSet1[,c(1:13)]
clin<-rownames_to_column(clin)
colnames(clin)[1]<-"Id"
ma_clin<-merge(clin,HNSC_ssgsva,by="Id")
ma_survival<-ma_clin[,c(1,7,8,15:21)]
rownames(ma_survival)<-ma_survival[,1]
ma_survival<-ma_survival[,-1]
pFilter=1 #设一个p值标准，后面用
outResult=data.frame() #建一个空白数据框，后面for循环输出用
sigGenes=c("OS","OS.time") #建一个向量，后面for循环输出用，因为后面还要用到surstat及surtime，所以先放在向量里
for(i in colnames(ma_survival[,3:ncol(ma_survival)])){ #从第3列开始循环，因为1列2列不是gene，是surstat和surtime
  tdcox <- coxph(Surv(OS.time, OS) ~ ma_survival[,i], data = ma_survival)#开始逐一循环cox分析
  tdcoxSummary = summary(tdcox) #summary命令对tdcox总结，方面后面提取数据
  pvalue=tdcoxSummary$coefficients[,"Pr(>|z|)"] #提取p值，这个主要是后面提取有意义的gene用
  if(pvalue<pFilter){ # 这里我们不需要提取所有基因的数据，只需要有意义的gene的结果，所以设置如果pvalue<0.05才提取数据
    sigGenes=c(sigGenes,i)
    outResult=rbind(outResult,#合并行，实际上是对循环结果的合并，前面设置的空白数据框outResult这里用，循环必须有个开始
                    cbind(id=i,#合并列，是每个基因的统计数据
                          HR=tdcoxSummary$conf.int[,"exp(coef)"],#提取单个基因的HR
                          L95CI=tdcoxSummary$conf.int[,"lower .95"],#提取单个基因的HR的95%CI低值
                          H95CI=tdcoxSummary$conf.int[,"upper .95"],#提取单个基因的HR的95%CI高值
                          pvalue=tdcoxSummary$coefficients[,"Pr(>|z|)"])#提取单个基因的p值
    )
  }
}


##反卷积
#安装MuSiC和Biobase
#BiocManager::install("TOAST")
#devtools::install_github('xuranw/MuSiC')
#BiocManager::install("Biobase",force = TRUE)
rm(list = ls())
library(Biobase)
library(MuSiC)
library(tidyverse)
##加载count数据和临床信息
load(file = "clinical/HNSC_exp_clin.Rdata")
test<-HNSC_exp_clin[1:153,1:20]
surv<-HNSC_exp_clin[,c(1:14)]
surv<-surv[-153,]

load(file = "R_data/0.1/HNSC_RNASEQ_exprdf_counts.Rdata")

#### class 查看属性
class(expr_df)
## 需要用as.data.frame来转换
exprSet <- as.data.frame(expr_df)
class(exprSet)
### 查看TCGA_id分组意义
### https://mp.weixin.qq.com/s/Ph1O6V5RkxkyrKpVmB5ODA
### 样本名称
TCGA_id <- colnames(exprSet)[-1]
table(substring(TCGA_id,14,15))
exprSet<-exprSet[,-c(269,324)]
### 创建分组信息
sample <- ifelse(substring(TCGA_id,14,15)=="01","Tumor","Normal")
sample <- factor(sample,levels = c("Tumor","Normal"),ordered = F)
### 获取配对信息，如果不是配对样本，就不需要这个信息
#paire_info <- as.factor(as.numeric(as.factor(substring(TCGA_id,1,12))))
### 创建metadata
data1 <- data.frame(TCGA_id,sample) 
test<-exprSet1[1:10,1:10]
exprSet1<-exprSet
rownames(exprSet1)<-exprSet1[,1]
exprSet1<-exprSet1[,-1]
exprSet1<-as.data.frame(t(exprSet1))
exprSet1<-rownames_to_column(exprSet1)
colnames(exprSet1)[1]<-"TCGA_id"
exprSet1<-merge(data1,exprSet1,by="TCGA_id")
group_list<-exprSet1[,"sample"]
exprSet2 =exprSet1[group_list=="Tumor",]
dim(exprSet2)
rownames(exprSet2)<-exprSet2[,1]
exprSet2<-exprSet2[,-c(1,2)]
test<-exprSet2[1:10,1:10]

exprSet3<-as.data.frame(t(exprSet2))
exprSet3<-rownames_to_column(exprSet3)
colnames(exprSet3)[1]<-"gene_id"
test<-exprSet3[1:10,1:10]
#探针转换
load(file = "resource/gtf_df.Rdata")

expr_df<-exprSet3
test<-expr_df[1:10,1:10]
### 提取编码基因(当然也可以提取非编码RNA)
library(dplyr)
res_2<- gtf_df %>% 
  ## 筛选gene,和编码指标
  dplyr::filter(type=="gene") %>%
  #dplyr::filter(type=="gene",gene_type=="protein_coding") %>%
  ## 选出基因名称，和ensemble id这两列
  dplyr::select(c(gene_name,gene_id)) %>% 
  ## 和表达量的数据交叉合并，等同于merge
  dplyr::inner_join(expr_df,by ="gene_id") %>% 
  ## 去掉多余列
  dplyr::select(-"gene_id") %>% 
  ## 以下是为了删除多于的行
  ## 增加一列
  mutate(rowMean = rowMeans(.[,-1])) %>% 
  ## 排序
  arrange(desc(rowMean)) %>% 
  ## 去重
  distinct(gene_name,.keep_all = T) %>% 
  ## 删除多余列
  dplyr::select(-rowMean)
expr_count_T<-res_2
text<-expr_count_T[1:10,1:10]
rownames(expr_count_T)<-expr_count_T[,1]
expr_count_T<-expr_count_T[,-1]
expr_count_T<-as.data.frame(t(expr_count_T))
expr_count_T<-expr_count_T%>%
  rownames_to_column("ID")
expr_count_T$TCGA_id<-substr(expr_count_T$ID,1,15)
data<-merge(surv,expr_count_T,by="TCGA_id")
test<-data[1:10,1:20]
surv<-data[,c(1:15)]
bulkexp<-data[,-c(2:15)]
rownames(bulkexp)<-bulkexp[,1]
bulkexp<-bulkexp[,-1]
bulkexp<-as.data.frame(t(bulkexp))
test<-bulkexp[1:2,1:2]
identical(surv$TCGA_id,colnames(bulkexp))
rownames(surv)=surv[,1]
surv<-surv[,-1]
surv1<-surv[,c(1:7,14)]
surv1<-surv1[,-6]
surv1$age<-as.numeric(surv1$age)
surv1$OS.time<-as.numeric(surv1$OS.time)
metadata <- data.frame(
  labelDescription =
    colnames(surv1),
  row.names = colnames(surv1))
surv1<-as.data.frame(surv1)
Anno <- new("AnnotatedDataFrame", data = surv1,varMetadata = metadata)
all(rownames(surv1) == colnames(bulkexp))
bulkexp<-as.matrix(bulkexp)
bulk <- new("ExpressionSet", exprs = bulkexp, phenoData = Anno)
bulk
bulk.mtx <- exprs(bulk)
test<-bulk.mtx[1:2,1:2]
bulk.mtx <- as.data.frame(bulk.mtx)



#开始跑反卷积
load(file = "4/noLP新分析/反卷积/MuSic_tcga_data.Rdata")
sce<-readRDS("4/noLP新分析/反卷积/sce_malig0.1_formusic.rds")

c.genes <- sort(intersect(rownames(bulk.mtx), rownames(sce_malig)))

sce_malig <- sce_malig[c.genes,]

bulk.mtx <- as.matrix(bulk.mtx[c.genes,])
# case.mtx <- as.matrix(case.mtx[c.genes,])

dim(sce_malig)
dim(bulk.mtx)
# dim(case.mtx)

Est.prop.TCGA <- music_prop(bulk.mtx = bulk.mtx, 
                            sc.sce = sce, 
                            clusters = "cell_type_me", 
                            #samples = "sample", 
                            verbose = F)###进行解卷积

Est.prop.TCGA <- music_prop(bulk.mtx = bulk.mtx,sc.sce = sce_malig,clusters = "SCT_snn_res.0.02",samples="tissue.type",verbose = F)
NNLS <- Est.prop.TCGA$Est.prop.allgene###MuSiC的结果有两种，源自于两种方法，NNLS和MuSiC
MuSiC <- Est.prop.TCGA$Est.prop.weighted###这两种结果都可以inner_join到临床数据中进行预后分析。

rm(list = ls())
library(pheatmap)
load(file = "4/only_CA/MuSiC/MuSiC_ONLY_CA_SCT_snn_res.0.02.Rdata")
load(file = "4/only_CA/MuSiC/NNLS_ONLY_CA_SCT_snn_res.0.02.Rdata")
load(file = "4/only_CA/MuSiC/MuSiC_SCT_snn_res.0.02_results.Rdata")

##一审补充
load(file = "E:/0/投稿/2_IMPDH/投稿IJOS/一审修回/Bulk部分补充/MuSiC/MuSiC_ONLY_CA_SCT_snn_res.0.02.Rdata")
load(file = "E:/0/投稿/2_IMPDH/投稿IJOS/一审修回/Bulk部分补充/MuSiC/MuSiC_SCT_snn_res.0.02_results.Rdata")
load(file = "E:/0/投稿/2_IMPDH/投稿IJOS/一审修回/Bulk部分补充/MuSiC/NNLS_ONLY_CA_SCT_snn_res.0.02.Rdata")
library(pheatmap)
p1 <- pheatmap(MuSiC, show_rownames = F,
               scale = "row", angle_col = "45",
               cluster_rows = T,cluster_cols = T,
               color = colorRampPalette(c("navy","white","firebrick"))(50))

p2 <- pheatmap(NNLS, show_rownames = F,
               scale = "row", angle_col = "45",
               cluster_rows = T,cluster_cols = T,
               color = colorRampPalette(c("navy","white","firebrick"))(50))


jitter.fig = Jitter_Est(list(data.matrix(Est.prop.TCGA$Est.prop.weighted),
                             data.matrix(Est.prop.TCGA$Est.prop.allgene)),
                        method.name = c('MuSiC', 'NNLS'), title = 'Jitter plot of Est Proportions')
jitter.fig
##丰度图
# 处理数据，排序
MuSiC1<-as.data.frame(NNLS)
colnames(MuSiC1)<-c("SG2","SG0","SG1")
MuSiC1<- MuSiC1 %>%  
  rownames_to_column("ID")
MuSiC1$SG3<-rep(0,499)
MuSiC1$SG4<-rep(0,499)
df<-MuSiC1[,c(1,2)]
df = df %>%
  arrange(df$SG2) %>%
  mutate("x" = row_number())

# 绘图
ggplot(df,aes(x = x,y = SG2))+
  geom_point(color="#FFA6AA")+ #"#BF3EEE"，,#7FFF00，#aee7f7，#FF8C00
  theme_bw() + theme(panel.grid=element_blank())+
  theme(axis.ticks = element_blank(),
        axis.text.y = element_blank())+
  #scale_y_continuous(expand = c(0,0))+
  coord_cartesian(ylim = c(0,1))

MuSiC2<-MuSiC1
rownames(MuSiC2)<-MuSiC2[,1]
MuSiC2<-MuSiC2[,-1]
colSums(MuSiC2)
colMeans(MuSiC2)
##按照malignant的cluster分组
rm(list = ls())
MuSiC<-as.data.frame(MuSiC)
load(file = "4/noLP新分析/反卷积/MuSiC_malignant_snnres_0.02.Rdata")
load(file = "4/noLP新分析/反卷积/NNLS_malignant_snnres_0.02.Rdata")
#load(file = "4/noLP新分析/反卷积/Est.prop.TCGA_malignant_cluster23.10.2.Rdata")
p1 <- pheatmap(MuSiC, show_rownames = F,
               scale = "row", angle_col = "45",
               cluster_rows = T,cluster_cols = T,
               color = colorRampPalette(c("navy","white","firebrick"))(50))

p2 <- pheatmap(NNLS, show_rownames = F,
               scale = "row", angle_col = "45",
               cluster_rows = T,cluster_cols = T,
               color = colorRampPalette(c("navy","white","firebrick"))(50))


jitter.fig = Jitter_Est(list(data.matrix(Est.prop.TCGA$Est.prop.weighted),
                             data.matrix(Est.prop.TCGA$Est.prop.allgene)),
                        method.name = c('MuSiC', 'NNLS'), title = 'Jitter plot of Est Proportions')
jitter.fig


##添加预后信息，计算单因素cox
load(file = "clinical/HNSC_exp_clin.Rdata")
test<-HNSC_exp_clin[1:153,1:20]
surv<-HNSC_exp_clin[,c(1:15)]
surv<-surv[-153,]
surv<-subset(surv,sample=="Tumor")
identical(surv$TCGA_id,rownames(MuSiC))
su_malig<-cbind(surv,MuSiC)
su_malig1<-su_malig[,-c(2:6,15)]

library("survival")
library("survminer")
library(dplyr)
library(tidyverse)
library(tidyr)
colnames(su_malig1)[10:12]<-c("C2","C0","C1")
covariates <- c("C2","C0","C1")
rownames(su_malig1)<-su_malig1[,1]

#分别对每一个变量，构建生存分析的公式

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(OS.time, OS)~', x)))
#对每一个特征做cox回归分析
univ_models <- lapply( univ_formulas, function(x){coxph(x, data = su_malig1)})
#提取HR，95%置信区间和p值
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         #获取p值
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         #获取HR
                         HR <-signif(x$coef[2], digits=2);
                         #获取95%置信区间
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(p.value,HR)
                         names(res)<-c("p.value","HR (95% CI for HR)")
                         return(res)
                       })
#转换成数据框，并转置

res_0.1 <- t(as.data.frame(univ_results, check.names = FALSE))
res_0.1 <-as.data.frame(res_0.1,stringsAsFactors=F)
write.table(file="4/5/univariate_cox_result.txt",as.data.frame(res),quote=F,sep="\t")
write.table(file="E:/0/投稿/2_IMPDH/投稿IJOS/一审修回/Bulk部分补充/MuSiC/univariate_cox_result_MuSiC.txt",as.data.frame(res_0.1),quote=F,sep="\t")

#单因素cox画图

HR=gsub("[\\(\\)]","",res_0.1$`HR (95% CI for HR)`)
HR=gsub("-"," ",HR)
HR=as.data.frame(do.call(cbind,strsplit(HR," ")),stringsAsFactors=F)
names(HR)=rownames(res_0.1)

#左边和右边边距稍微留多一点来写变量名称，pvalue和HR
par(mar=c(5,6,4,13))
#先用小方块画出HR
plot(as.numeric(HR[1,]),1:dim(HR)[2],
     pch=15,cex=2,col="#000080",bty='n',yaxt='n',ylab=NA,xlab="Hazard Ratio",
     xlim=range(as.numeric(unlist(HR)))
)
#添加中线
abline(v=1,col="grey",lwd=2,lty=2)
for(i in 1:ncol(HR)){
  x=as.numeric(HR[2:3,i])
  #循环画出CI
  lines(x,c(i,i),col="#000080")
  #添加变量名
  text(0,i,rownames(res_0.1)[i],xpd=T,adj = c(0,0))
  #添加p值
  text(5,i,as.numeric(res_0.1[i,1]),xpd=T,adj = c(0,0))
  #添加HR和CI
  text(9,i,as.character(res_0.1[i,2]),xpd=T,adj = c(0,0))
}
#添加标题
text(5,ncol(HR)+0.2,"pvalue",xpd=T,adj = c(0,0))
text(9,ncol(HR)+0.2,"HR(CI)",xpd=T,adj = c(0,0))

##多因素COX
res.cox <- coxph(Surv(OS.time, OS) ~ C2 + C0, data =su_malig1)
x <- summary(res.cox)
pvalue=signif(as.matrix(x$coefficients)[,5],2)
HR=signif(as.matrix(x$coefficients)[,2],2)
HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
HR <- paste0(HR, " (", 
             HR.confint.lower, "-", HR.confint.upper, ")")

res_muti=data.frame(p.value=pvalue,
                     HR=HR,
                     stringsAsFactors = F)
names(res_muti)<-c("p.value","HR (95% CI for HR)")
##以上多因素COX没有结果,不做了
su_malig1$time<-su_malig1$OS.time/365
data.sur<-su_malig1[,c(1,11,2,13)]
names(data.sur) <- c('ID','score','status','time')
b <- surv_cutpoint(data.sur,time = 'time',event = 'status',variables = 'score',minprop = 0.1)
cutpoint <- summary(b)$cutpoint
high <- length(which(data.sur$score > cutpoint))
data.sur2 <- data.sur[order(data.sur$score,decreasing = T),]
rownames(data.sur2) <- data.sur2$ID
data.sur2 <- cbind(data.sur2,group = rep(c('High','Low'),c(high,499-high)))

data.sur2$OS.time<-data.sur2$time*365
coxdata <- data.sur2 %>% 
  filter(OS.time >= 30) 

### Surv()：创建一个生存对象，时间/事件
rt <- data.frame(coxdata[,c(4,3)],riskScore=coxdata[,5])
mySurv <- Surv(rt$time, rt$status) 
group<-rt$riskScore
### 使用survfit函数拟合一条生存曲线
### 使用summary函数查看模型汇总结果
sfit <- survfit(formula = mySurv~group, data = rt)
summary<-summary(sfit)

### 画Kaplan-Meier曲线
#   ggsurvplot_list() 绘制多个对象
#   ggsurvplot_facet() 分面到多个panels
#   ggsurvplot_group_by() 一幅图中多个分组
#   ggsurvplot_add_all() 总合所有的情况
#   ggsurvplot_combine() 一个图中结合多个survfit对象

library(survminer) 
ggsurvplot(sfit, 
           conf.int=TRUE, pval=TRUE, risk.table=TRUE, 
           #legend.labs=c("H1", "H2","L1","L2","M"), 
           #legend.title="1C",  
           #palette=c("dodgerblue2", "orchid2","green", "red","black"), 
           title="Kaplan-Meier Curve", 
           risk.table.height=.15)
library(export)
ggsurvplot(sfit, 
           size=1.25,#曲线粗细
           legend = c(0.8,0.8), # 指定图例strata位置
           #surv.plot.height=0.25,# 生存图的高度，默认为0.75；# 当risk.table = FALSE时忽略
           ggtheme = theme_classic(base_size = 20), #想要网格就运行这行
           conf.int = F, #不画置信区间，想画置信区间就把F改成T
           #conf.int.style = "step",#置信区间的类型，还可改为ribbon
           censor = F, #不显示观察值所在的位置
           #surv.median.line = "hv", #添加中位生存曲线
           linetype = "strata",# 改变曲线的类型
           censor.shape="+", censor.size = 4,#更改删失点形状及大小，默认为"+", 可选"|"
           #palette = c("#E7B800", "#2E9FDF"), #线的颜色对应高、低
           ylab="Overall survival",xlab = " Time (years)", # xlab, ylab 分别指x轴和y轴标签更改横纵坐标
           xlim = c(0,5),# xlim, ylim # 指定x轴和y轴的范围，如xlim = c(0,30), ylim = c(0,1) 
           break.x.by = 1,#横坐标间隔
           #risk.table = TRUE,tables.height = 0.32,#加risk table
           font.legend = c(14, "bold","black"),#图例的字体大小
           font.title = c(14,"bold"),font.x = c(14,"bold"),font.y = c(14,"bold"),#设置其他字体大小
           font.xtickslab = c(14, "bold","black"),font.ytickslab = c(14, "bold","black"),#坐标轴刻度线字体
           #在左下角标出pvalue、HR、95% CI
           #太小的p value标为p < 0.001
           pval = TRUE)


###limma分析
#install.packages("E:/0/edgeR_3.42.4.zip",
                 #repos = NULL,type = "win.binary")
library(dplyr)
library(tidyverse)
library(tidyr)
load(file = "R_data/0.1/HNSC_RNASEQ_exprdf_counts.Rdata")
## 需要用as.data.frame来转换
exprSet <- as.data.frame(expr_df)
class(exprSet)
### 查看TCGA_id分组意义
### https://mp.weixin.qq.com/s/Ph1O6V5RkxkyrKpVmB5ODA
### 样本名称
TCGA_id <- colnames(exprSet)[-1]
table(substring(TCGA_id,14,15))
exprSet<-exprSet[,-c(269,324)]

load(file = "R_data/HNSC_RNASEQ_exprSet_counts_499T.Rdata")
test<-exprSet[1:10,1:10]
### 创建分组信息
sample <- ifelse(substring(TCGA_id,14,15)=="01","Tumor","Normal")
sample <- factor(sample,levels = c("Tumor","Normal"),ordered = F)
####筛选肿瘤样本
### 获取配对信息，如果不是配对样本，就不需要这个信息
#paire_info <- as.factor(as.numeric(as.factor(substring(TCGA_id,1,12))))
### 创建metadata
data1 <- data.frame(TCGA_id,sample) 
test<-exprSet1[1:10,1:10]
exprSet1<-exprSet
rownames(exprSet1)<-exprSet1[,1]
exprSet1<-exprSet1[,-1]
exprSet1<-as.data.frame(t(exprSet1))
exprSet1<-rownames_to_column(exprSet1)
colnames(exprSet1)[1]<-"TCGA_id"
exprSet1<-merge(data1,exprSet1,by="TCGA_id")
group_list<-exprSet1[,"sample"]
exprSet2 =exprSet1[group_list=="Tumor",]
dim(exprSet2)
rownames(exprSet2)<-exprSet2[,1]
exprSet2<-exprSet2[,-c(1,2)]
test<-exprSet2[1:10,1:10]
exprSet3<-as.data.frame(t(exprSet2))
exprSet3<-rownames_to_column(exprSet3)
colnames(exprSet3)[1]<-"gene_id"
#exprSet3<-exprSet3[,-154]
test<-exprSet3[1:10,1:10]
load(file = "IJOS投稿/HNSC_cox_riskscore_total_TCGA_SGOC_5gene.Rdata")
load(file="4/only_CA/重新5/HNSC_cox_riskscore_total_TCGA_SGOC.Rdata")
exprSet2<-rownames_to_column(exprSet2)
exprSet2$ID<-substr(exprSet2$rowname,1,15)
exprSet4<-merge(data.sur2,exprSet2,by="ID")
test<-exprSet4[1:10,1:10]
metadata<-exprSet4[,c(6,5)]
colnames(metadata)<-c("TCGA_id","sample")
metadata$sample<-as.factor(metadata$sample)

load("4/metadata_tumor_3Gene_9.15.Rdata")
library(DESeq2)
exprSet3<-exprSet4[,-c(1:5)]
test<-exprSet3[1:10,1:10]
rownames(exprSet3)<-exprSet3[,1]
exprSet3<-exprSet3[,-1]
exprSet3<-as.data.frame(t(exprSet3))
exprSet3<-exprSet3%>%
  rownames_to_column("gene")
dds <-DESeqDataSetFromMatrix(countData=exprSet3, 
                             colData=metadata, 
                             design=~sample,
                             tidy=TRUE)
### 获取行数
nrow(dds)
### 过滤
### 如果一个基因在所有样本中的counts数小于等于1，我们就把他删掉
dds <- dds[rowSums(counts(dds))>1,]
### 获取行数
nrow(dds)
### 最困难的一步来了，DESeq2主程序
dds <- DESeq(dds)

normalized_counts <- as.data.frame(counts(dds, normalized=TRUE))


load(file = "resource/gtf_df.Rdata")
test <- expr_df[1:10,1:10]
expr_df<-normalized_counts
expr_df<-rownames_to_column(expr_df)
colnames(expr_df)[1]<-"gene_id"

### 提取编码基因(当然也可以提取非编码RNA)
library(dplyr)
res_2<- gtf_df %>% 
  ## 筛选gene,和编码指标
  #dplyr::filter(type=="gene") %>%
  dplyr::filter(type=="gene",gene_type=="protein_coding") %>%
  ## 选出基因名称，和ensemble id这两列
  dplyr::select(c(gene_name,gene_id)) %>% 
  ## 和表达量的数据交叉合并，等同于merge
  dplyr::inner_join(expr_df,by ="gene_id") %>% 
  ## 去掉多余列
  dplyr::select(-"gene_id") %>% 
  ## 以下是为了删除多于的行
  ## 增加一列
  mutate(rowMean = rowMeans(.[,-1])) %>% 
  ## 排序
  arrange(order(rowMean,decreasing = T)) %>% 
  ## 去重
  distinct(gene_name,.keep_all = T) %>% 
  ## 删除多余列
  dplyr::select(-rowMean) 

expr_count<-res_2
rownames(expr_count)<-expr_count[,1]
expr_count<-expr_count[,-1]
test<-expr_count[1:10,1:10]







library(limma)
library(edgeR)
### 加载数据，注意解决报错
#load(file = "output/exprSet_rmdup.Rdata")

### 1.创建分组
### 这一步根据样本来就行，原则就是: 跟样本匹配，取决于样本的排序
metadata$sample<-as.character(metadata$sample)
group <- metadata[,2]
##group <- c("con","con","treat","con","treat","treat") 
### 分组变成向量，并且限定leves的顺序
### levels里面，把对照组放在前面
group <- factor(group,levels = c("Low","High"),ordered = F)
### 构建比较矩阵
design <- model.matrix(~group)
### 比较矩阵命名
colnames(design) <- levels(group)
rownames(design)=colnames(expr_count)
deg<-DGEList(counts = expr_count)
deg<-calcNormFactors(deg)
#design
de<-voom(deg,design,plot = TRUE,normalize="quantile")
### 2.线性模型拟合
fit <- lmFit(de,design)
### 3.贝叶斯检验
fit2 <- eBayes(fit)
### 4.输出差异分析结果,其中coef的数目不能操过design的列数
### 此处的2代表的是design中第二列和第一列的比较
DEGene=topTable(fit2,adjust='fdr',coef=2,number=Inf) 
DEGene_1<-rownames_to_column(DEGene)
write.csv(DEGene,file = "IJOS投稿/HNSC_High_Low_protein_5G_allDiff.csv")

write.csv(DEGene,file = "4/5/HNSC_High_Low_protein_TOP15_0.45_4gene_allDiff.csv")
load(file = "4/6/HNSC_High_Low_protein_4G_allDiff_9.16.Rdata")


################################################

### 本节任务：实操GO分析，KEGG分析
rm(list = ls())
### 加载差异基因列表
load(file = "4/6/HNSC_High_Low_protein_4G_allDiff_9.16.Rdata")
load(file = "4/6/tougao_lasso_marker_2_logfc0.3.Rdata")
load(file = "4/only_CA/H_L_markergene_30%.Rdata")
load(file = "Fig4/H_L_markergene_30%.Rdata")

### 筛选差异基因
library(dplyr)
library(tibble)
diffgene_H <- H_L_markergene1 %>% 
  filter(p_val_adj < 0.05) %>% 
  filter(avg_log2FC> 0)
write.csv(diffgene,"4/6/1327ge-diffgene_logFC_0.59_HvsL.csv")

library(clusterProfiler)
library(org.Hs.eg.db)
#load(file = "output/选取表达有意义基因5/P0.01分析/diffgene_logFC0.59_bestcutpoint_0.4.Rdata")
### 这个分析需要什么数据？
### 获得基因列表
diffgene<-rownames_to_column(H_L_markergene1)
colnames(H_L_markergene1)[1]<-"gene"
gene <- diffgene_ALL$gene
#基因名称转换，返回的是数据框
gene = bitr(gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
head(gene)

#################################################
#GO分析三个大类，细胞组分，cellular compartment

ego_CC <- enrichGO(gene = gene$ENTREZID,
                   OrgDb= org.Hs.eg.db,
                   ont = "CC",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.01,
                   readable = TRUE)


cc<-as.data.frame(ego_CC)
write.table(file="output/GO_CC_differgene_logFC_9.527.txt",as.data.frame(cc),quote=F,sep="\t")
#
#################################################
#GO分析三个大类，生物过程BP,biological process
ego_BP_L <- enrichGO(gene = gene$ENTREZID,
                   OrgDb= org.Hs.eg.db,
                   ont = "BP",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.01,
                   readable = TRUE)
bp_high<-as.data.frame(ego_BP_H)
bp_L<-as.data.frame(ego_BP_L)
write.table(file="4/only_CA/GO_BP_L_scRNA_differgene_logFC_0.3.txt",as.data.frame(bp_L),quote=F,sep="\t")
#

##########################################################################
###选取免疫画图
ego_result_BP <- as.data.frame(ego_BP)[c(8,9,17,18,20,25,27,30:32,36,37,41:43,53,56,59:61), ]
ego_result_BP <- as.data.frame(ego_BP)[c(1,3,5,6,10,12,13,15,17,18), ]
BP<-as.data.frame(data.table::fread("4/6/GO_BP_immune_10.txt",header = T,sep="\t"))
rownames(BP)<-BP[,1]
BP<-BP[,-1]
go_enrich_df <- data.frame(
  ID=BP$ID, 
  Description=BP$Description,
  GeneNumber=BP$Count,
  type=factor(rep("biological process", 10)), 
  levels="biological process")

###横着的柱状图
BP<-BP[order(BP$Count),]
rownames(BP) <- 1:nrow(BP)

ego_result_BP1<-BP[1:10,]
ego_result_BP1<-ego_result_BP1[order(ego_result_BP1$Count),]
ego_result_BP1$order=factor(rev(as.integer(rownames(ego_result_BP1))),labels = rev(ego_result_BP1$Description))
library(ggplot2)
ggplot(ego_result_BP1,aes(y=order,x=Count,fill=p.adjust))+
  geom_bar(stat = "identity",width=0.7)+####柱子宽度
  #coord_flip()+##颠倒横纵轴
  scale_fill_gradient(low = "red",high ="blue" )+#颜色自己可以换
  #labs(title = "biological process",
       #x = "Gene numbers", 
       #y = "Pathways")+
  #theme(axis.title.x = element_text(face = "bold",size = 16),
        #axis.title.y = element_text(face = "bold",size = 16),
        #legend.title = element_text(face = "bold",size = 16))+
  theme_bw()

#################################################
#GO分析三个大类，分子功能MF, Molecular function
ego_MF <- enrichGO(gene = gene$ENTREZID,
                   OrgDb= org.Hs.eg.db,
                   ont = "MF",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.01,
                   readable = TRUE)

mf<-as.data.frame(ego_MF)

GOplotIn_BP<-ego_BP_L[c(15,23,33,68,90,120,130,170,194,198),c(1,2,6,8)]
library(stringr)
GOplotIn_BP$geneID <-str_replace_all(GOplotIn_BP$geneID,'/',',') #把GeneID列中的’/’替换成‘,’
names(GOplotIn_BP)<-c('ID','term','adj_pval','genes')#修改列名,后面弦图绘制的时候需要这样的格式
GOplotIn_BP$category = "BP"#分类信息
GOplotIn_BP<-rownames_to_column(GOplotIn_BP)
GOplotIn_BP<-GOplotIn_BP[,-1]
genedata<-data.frame(ID=diffgene_L$gene,logFC=-diffgene_L$avg_log2FC)
#colnames(genedata)<-c("ID","logFC")
circ_BP<-GOplot::circle_dat(GOplotIn_BP,genedata) 
#circ_BP <- circle_dat(GOplotIn_BP,genedata)#GOplot导入数据格式整理
GOCircle(circ_BP,lfc.col = "red") #弦表图

#high-risk group
GOplotIn_BP<-ego_BP_H[c(9,54,66,69,115,125,163,224,250,303),c(2,3,7,9)]
library(stringr)
GOplotIn_BP$geneID <-str_replace_all(GOplotIn_BP$geneID,'/',',') #把GeneID列中的’/’替换成‘,’
names(GOplotIn_BP)<-c('ID','term','adj_pval','genes')#修改列名,后面弦图绘制的时候需要这样的格式
GOplotIn_BP$category = "BP"#分类信息
GOplotIn_BP<-rownames_to_column(GOplotIn_BP)
GOplotIn_BP<-GOplotIn_BP[,-1]
genedata<-data.frame(ID=diffgene_H$rowname,logFC=diffgene_H$avg_log2FC)
#colnames(genedata)<-c("ID","logFC")
circ_BP<-GOplot::circle_dat(GOplotIn_BP,genedata) 
#circ_BP <- circle_dat(GOplotIn_BP,genedata)#GOplot导入数据格式整理
GOCircle(circ_BP,lfc.col = "red") #弦表图
##############################################################
#**KEGG分析**
##############################################################
library(KEGG.db)
EGG <- enrichKEGG(gene = gene$ENTREZID,
                  organism = 'hsa',
                  keyType = 'kegg',
                  pvalueCutoff = 0.05,
                  use_internal_data =T)
#head(EGG,2)
kk=setReadable(EGG,OrgDb=org.Hs.eg.db,keyType  ="ENTREZID")
egg<-as.data.frame(EGG)
write.table(file="4/6/KEGG_scRNA_differgene_logFC0.3.txt",as.data.frame(egg),quote=F,sep="\t")
#


#KEGG弦图
kegg=data.frame(kk)[,c("ID","Description","geneID","p.adjust" )]
#kegg<-kk[c(1:3,9,15,20,26,29,30),c(1,2,6,8)]
#KEGG富集分析结果中没有category这里一列，这里我们自己造一个
kegg$category=rep("kegg",nrow(kegg))

#重新给列命名，满足circle_dat函数的要求
names(kegg)=c("ID" ,"Term", "Genes","adj_pval","Category")
#将/隔开的基因名字转换成,隔开
kegg$Genes=gsub("/",",",kegg$Genes)
#将差异表达分析结果的第一列symbol改成ID，满足circle_dat函数的要求
genedata<-data.frame(ID=diffgene_ALL$gene,logFC=diffgene_ALL$avg_log2FC)



#创建circ对象，第一个参数为KEGG富集分析结果，第二个参数为差异表达分析结果
circ <- circle_dat(kegg, genedata)

#根据差异表达基因的logFC排序
genes=genedata$logFC
names(genes)=genedata$ID
genes=sort(genes)
#挑选logFC最大和最小的各100各基因，100这个数字可以自己修改
#也可以自己挑选基因，需要带上logFC
chord_genes=c(head(genes,200),tail(genes,200))
#构建数据框，第一列为基因名字，第二列为对应的logFC
chord_genes_df=data.frame(ID=names(chord_genes),logFC=chord_genes)
#构建chord对象，第一个参数为circ对象，第二个参数为带有logFC的基因名字，
#第三个是想展示的KEGG通路的名字
egg<-rownames_to_column(egg)
process<-kegg[c(1:3,6,13,21,23,20,30,19),]
process<-process$Term
chord <- chord_dat(circ, chord_genes_df,process)
library(DOSE)
#创建一个pdf文件来保存图片
pdf("4/6/KEGG_chord.pdf",height = 22,width = 20)
source("Helper.R")
KEGGChord(chord,   #chord对象
          limit=c(0,0), #第一个数每个基因至少需要根几个term相连，第二个数每个term至少需要根几个基因相连
          space = 0.01,  #右侧色块之间的间距
          gene.order = 'logFC',   #基因展示顺序根据logFC来
          gene.space = 0.25,  #基因名字和色块之间的距离
          gene.size = 4 ) #基因名字大小
          


library(stringr)
eggplotIn$geneID <-str_replace_all(eggplotIn$geneID,'/',',') #把GeneID列中的’/’替换成‘,’
names(eggplotIn)<-c('ID','term','adj_pval','genes')#修改列名,后面弦图绘制的时候需要这样的格式
eggplotIn$category = "BP"#分类信息
eggplotIn<-rownames_to_column(eggplotIn)
eggplotIn<-eggplotIn[,-1]
#colnames(genedata)<-c("ID","logFC")
circ_kegg<-GOplot::circle_dat(eggplotIn,genedata) 
chord_kegg<-chord_dat(data = circ_kegg,genes = genedata) #生成含有选定基因的数据框
#circ_BP <- circle_dat(GOplotIn_BP,genedata)#GOplot导入数据格式整理
GOChord(data = circ_kegg,#弦图
        title = 'GO-Biological Process',space = 0.01,#GO Term间距
        limit = c(1,1),gene.order = 'logFC',gene.space = 0.25,gene.size = 5,
        lfc.col = c('red','white','blue'), #上下调基因颜色
        process.label = 10) #GO Term字体大小

#提取数据；
library(stringr)
library(dplyr)
library(tidyr)
library(reshape2)
dt <- as_tibble(eggplotIn[c(2,4)])
#转成因子；
dt$Description <- as.factor(dt$Description)
#根据逗号进行分隔；
df <- str_split(dt$geneID,",",n=7,simplify=TRUE)
#生成数据框；
dt2 <- data.frame(dt[,1],df)
dt2
#将空值转成NA;
for (i in 2:7) {
  dt2[,i] <- str_replace_all(dt2[,i], pattern="^$", NA_character_)
}
dt2<-dt2[,c(1:7)]
#转成关系对数据；
dt3 <- pivot_longer(dt2,!Description, names_to = "index",
                    values_to = "gene",
                    values_drop_na = TRUE)
dt3
merge<-merge(dt3,gene,by.x="gene",by.y="ENTREZID")
#############################################################################
## 这个图片可以订制
## rich factor
## 定制一个图片
## 内置的函数可以转换为数据框
df = as.data.frame(EGG)
dd = EGG@result
dd<-as.data.frame(data.table::fread("4/6/KEGG_immune_10.txt",header = T,sep="\t"))

## 计算富集分数,并增加为列，两个数相除
head(dd$BgRatio)
## 提取/前面的字符串
as.numeric(sub("/\\d+", "", dd$BgRatio))
dd$richFactor =dd$Count / as.numeric(sub("/\\d+", "", dd$BgRatio))
## 提取p值小于0.05 的
dd <- dd[dd$p.adjust < 0.05,]
dd1<-dd#[c(2,8,10,12,15,25,27,29,34,37),]
library(ggplot2)
## 正式画图
ggplot(dd1,aes(richFactor,forcats::fct_reorder(Description, richFactor))) + 
  ## 画横线
  geom_segment(aes(xend=0, yend = Description)) +
  ## 画点
  geom_point(aes(color=p.adjust, size = Count)) +
  ## 调整颜色的区间,begin越大，整体颜色越明艳
  scale_color_viridis_c(begin = 0.3, end = 1) +
  ## 调整泡泡的大小
  scale_size_continuous(range=c(2, 10)) +
  theme_bw() + 
  xlab("Rich factor") +
  ylab(NULL) + 
  ggtitle("")



####################################################################
###################################################################
###################制作所有免疫评分
###immune pathway score
#制作gmt文件
rm(list = ls())
library(tidyverse)
library(tidyr)
library(dplyr)
rt <- read.table(file = "resource/iES.txt",header = T,sep = "\t")
name <- unique(rt$type)
description <- rep(NA,length(name))
names(description) <- name
genes <- lapply(name, function(name){
  as.vector(rt[rt$type == name,"name"])
})
names(genes) <- name

gmtinput <- list(name=name,description=description,genes=genes)
get_gmt <- function(gmtinput,filename){
  output <- file(filename, open="wt")
  lapply(gmtinput[["name"]],function(name){
    outlines = paste0(c(name, gmtinput[["description"]][[name]],
                        gmtinput[["genes"]][[name]]),collapse='\t')
    writeLines(outlines, con=output)
  })
  close(output)
}
get_gmt(gmtinput=gmtinput,filename="R_data/immune_subtype.gmt")

### 本节任务: GSVA
rm(list = ls())
library(dplyr)
library(tibble)
library(RColorBrewer)
library(GSVA)
library(GSEABase)
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(pheatmap)
library(TCGAbiolinks)
library(limma)
library(edgeR)
library(readr)

### 加载表达数据，这是vst标准化后的数据
load(file = "E:/0/1/2/3/TCGA/clinical/HNSC_exp_clin_tumor_500.Rdata")
rownames(exprSet)<-exprSet[,1]
exprSet<-exprSet[,-c(1:15)]
test<-exprSet[1:10,1:10]
exprSet1<-as.data.frame(t(exprSet))

### 加载基因list,整理成标准的gmt格式文件
#################################################################################################################
### 自制基因集GMT文件制作，IMMUNE
pathway <- read_delim("R_data/immune_subtype.gmt", "\t", escape_double = FALSE, 
                      trim_ws = TRUE,col_names = F,show_col_types = FALSE)
pathway <- as.data.frame(pathway)

rownames(pathway)<-pathway$X1
pathway <- data.frame(t(pathway))
pathway <- pathway[3:202,]
pathway_list <- vector("list",length(pathway))
pathway_list <- lapply(pathway, function(x) {
  unique(na.omit(x)) 
})

gsva_pathway <- gsva(as.matrix(exprSet1), pathway_list,method='ssgsea',
                     kcdf='Gaussian',abs.ranking=TRUE)

gsva_pathway<-as.data.frame(t(gsva_pathway))
save(gsva_pathway,file="1/data/gsva_immune_subtype_immunity2018_HNSC_high_low_ssgsea_9.527.Rdata")
####################################################################
###ESTIMATE score
rm(list = ls())
library(dplyr)
library(tibble)
library(tidyverse)
library(tidyr)
###制作TXT文件
##加载tpm数据
load("R_data/exprSet_protein_coding_gene_exprssion_tpm_499.Rdata")
sample_input_estimate<-as.data.frame(t(sample))
test<-sample_input_estimate[1:10,1:10]
sample_input_estimate<-rownames_to_column(sample_input_estimate)
colnames(sample_input_estimate)[1]<-"GeneSymbol"
write.table(sample_input_estimate, 
            file = '1/data/sample_input_estimate.txt', sep="\t",row.names =F, quote = F)
test<-sample_input_estimate[1:10,1:10]
#安装
#library(utils)
#rforge <- "http://r-forge.r-project.org"
#install.packages("estimate", repos=rforge, dependencies=TRUE)
library(estimate)
help(package="estimate")

rm(list = ls())
library(estimate)
in.file <- '1/data/sample_input_estimate.txt' #输入文件
outfile2E <- 'ESTIMATE_input.gct' #生成ESTIMATE 的输入文件
outputGCT(in.file, outfile2E) #该函数以GCT格式写入输入文件

filterCommonGenes(input.f= in.file, output.f= outfile2E, id="GeneSymbol")

# 该功能将每个平台的不同数量的基因与10412个普通基因相结合。
### code chunk number 2: estimate
#这个功能计算基质，免疫，并估计得分每个样本使用基因表达数据。

estimateScore("ESTIMATE_input.gct", "ESTIMATE_score_4g_9.16.gct")
plotPurity(scores="ESTIMATE_score_4g_9.16.gct", platform="affymetrix")#,samples="s516")
#根据ESTIMATE score绘制肿瘤纯度。

#将评分保存为txt格式 
ESTIMATE_score <- read.table("ESTIMATE_score_4g_9.16.gct", skip = 2,#前两行跳过 
                             header = TRUE,row.names = 1) 
ESTIMATE_score <- ESTIMATE_score[,2:ncol(ESTIMATE_score)] 
ESTIMATE_score 
write.table(data.frame(scores=rownames(ESTIMATE_score),ESTIMATE_score,
                       check.names = F),file = "1/data/ESTIMATE_score_4g.txt",
            row.names =F,quote = F,sep = "\t")

####################################################################
###X-cell score
rm(list = ls())
library(xCell)
library(dplyr)
library(tibble)
library(tidyverse)
library(tidyr)
###############制作xcell需要的数据
##加载tpm数据(9.527与之前一样)
load(file = "R_data/exprSet_protein_coding_gene_exprssion_tpm_499.Rdata")
sample_input_xcell<-as.data.frame(t(sample))
test<-sample_input_xcell[1:10,1:10]
save(sample_input_xcell,file = "R_data/sample_input_xcell.Rdata")

#####
load(file = "R_data/sample_input_xcell.Rdata")
exprSet<-sample_input_xcell
xCell_RNAseq_score<- xCellAnalysis(exprSet,rnaseq = T)
xCell_RNAseq_score<-as.data.frame(xCell_RNAseq_score)
save(xCell_RNAseq_score,file = "4/data/xCell_RNAseq_score.Rdata")
immune_xcell_score<-xCell_RNAseq_score[c(65,66,67),]
immune_xcell_score<-as.data.frame(t(immune_xcell_score))
save(immune_xcell_score,file = "1/data/immune_xcell_score.Rdata")
xCell_RNAseq_score1<-rownames_to_column(xCell_RNAseq_score)
xCELL_immunecell<-xCell_RNAseq_score[c(4,6:14,17,19,21,28,32:35,38,46,49,50,57,61:64),]
save(xCELL_immunecell,file = "4/data/xCELL_immunecell_heatmap.Rdata")
####################################################################
###MHC CYT score
rm(list = ls())
library(dplyr)
library(tibble)
library(tidyverse)
library(tidyr)
##加载vst数据
load(file = "E:/0/1/2/3/TCGA/clinical/HNSC_exp_clin_tumor_500.Rdata")
rownames(exprSet)<-exprSet[,1]
exprSet<-exprSet[,-c(1:15)]
exprSet1<-as.data.frame(t(exprSet))
test<-exprSet[1:10,1:10]
MHC_gene<-exprSet1[c("HLA.A","HLA.B","HLA.C","TAP1","TAP2",
                    "NLRC5","PSMB9","PSMB8","B2M"),]
cyt_gene<-exprSet1[c("GZMA","PRF1"),]
MNC_meaan<-t(colMeans(MHC_gene))
MNC_meaan<-t(MNC_meaan)
colnames(MNC_meaan)<-"MHC"
MHC_meaan<-as.data.frame(MNC_meaan)
save(MHC_meaan,file = "1/data/MHC_mean_clindata_4g.Rdata")
cyt_meaan<-t(colMeans(cyt_gene))
cyt_meaan<-t(cyt_meaan)
colnames(cyt_meaan)<-"CYT"
cyt_meaan<-as.data.frame(cyt_meaan)
save(cyt_meaan,file = "1/data/cyt_mean_clindata_4g.Rdata")
###############################################################
###data combind
rm(list = ls())
library(dplyr)
library(tibble)
library(tidyverse)
library(tidyr)
#######################################
#加载分组数据
load(file = "1/metadata_tumor_4Gene_9.16.Rdata")
colnames(metadata)[1]<-"ID"
metadata$TCGA_id<-substr(metadata$ID,1,15)
#加载ESTIMATE数据
estimate<-read.table(file = "1/data/ESTIMATE_score_4g.txt")
#加载xcell——score数据
load(file = "1/data/immune_xcell_score.Rdata")
#加载ssGSEA-iES数据
load(file = "1/data/gsva_immune_subtype_immunity2018_HNSC_high_low_ssgsea_9.527.Rdata")
#加载MHC数据
load(file = "1/data/MHC_mean_clindata_4g.Rdata")
#加载CYT数据
load(file = "1/data/cyt_mean_clindata_4g.Rdata")
##############
#将这些数据合并
identical(metadata[,1], rownames(cyt_meaan))

colnames(estimate)<-estimate[1,]
estimate<-estimate[-1,]
rownames(estimate)<-estimate[,1]
estimate<-estimate[,-1]
estimate<-as.data.frame(t(estimate))
save(estimate,file = "1/data/estimate_HNSC_immune_4g.Rdata")
data1<-cbind(metadata,estimate)
data2<-cbind(data1,immune_xcell_score)
gsva_pathway<-gsva_pathway %>% 
  rownames_to_column("TCGA_id")
data3<-merge(data2,gsva_pathway,by="TCGA_id")
colnames(data3)[8]<-"ImmuneScore_xcell"
identical(rownames(cyt_meaan),data3[,1])
data4<-cbind(data3,cyt_meaan)
data5<-cbind(data4,MHC_meaan)
boxpolt_data<-data5
save(boxpolt_data,file = "1/data/boxpolt_data.Rdata")




#############immunedeconv
#install.packages("remotes")
#remotes::install_github("icbi-lab/immunedeconv")

# 导入所需R包
library('dplyr')
library('ggplot2')
library('tidyr')
library('immunedeconv')
library('tibble')
load(file = "R_data/exprSet_protein_coding_gene_exprssion_tpm_499.Rdata")
tpm_tumor<-as.data.frame(t(sample))
test<-tpm_tumor[1:10,1:10]
logtpm<-log(tpm_tumor)
deconvolution_methods
res_quantiseq = deconvolute(logtpm, "quantiseq", tumor = TRUE)
save(res_quantiseq,file = "R_data/res_quantiseq_heatmap.Rdata",tumor = TRUE)
res_mcp_counter = deconvolute(logtpm, "mcp_counter",tumor = TRUE)
save(res_mcp_counter,file = "output/0.7/res_mcp_counter_heatmap.Rdata")
res_epic = deconvolute(logtpm, "epic",tumor = TRUE)
save(res_epic,file = "output/0.7/res_epic_heatmap.Rdata")

############################################################################
###TMB CNV
#下载突变数据并整理
rm(list = ls())
library(tidyr)
library(dplyr)
library(tidyverse)
#BiocManager::install("maftools")
library(maftools)
library(TCGAbiolinks)
query_SNV <- GDCquery(project = "TCGA-HNSC",
                      data.category = "Simple Nucleotide Variation",
                      data.type = "Masked Somatic Mutation",
                      workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking")
GDCdownload(query_SNV)


library(maftools)
library(tidyverse)
#下载的数据有一个是不可以读进去的，所以删除|].
mafFilePath2 = dir(path = "GDCdata/TCGA-HNSC",pattern = "masked.maf.gz$",full.names = T,recursive=T)
mafdata2 <- lapply(mafFilePath2, function(x){read.maf(x,isTCGA=TRUE)})
snv_data = merge_mafs(mafdata2)
save(snv_data,file = "output/0.7/snv_data.Rdata")
#不可读da<-read.maf("GDCdata/TCGA-HNSC/harmonized/Simple_Nucleotide_Variation/Masked_Somatic_Mutation/62ac65bf-fa95-4c46-93e4-16bdd8e0543d/ae7541bc-7566-4be9-921a-b8bdb7049513.wxs.aliquot_ensemble_masked.maf.gz",isTCGA = TRUE)
load(file = "output/0.7/snv_data.Rdata")
HNSC_snv_data<-snv_data@data
#看整体分布
plotmafSummary(maf = snv_data, rmOutlier = TRUE, addStat = 'median')
library(export)

#看瀑布图
oncoplot(maf = snv_data, top = 10) # 高频突变的前10个基因



#临床信息
load(file = "clinical/clin_data_HNSC_deal_WLX.Rdata")
load(file = "1/metadata_tumor_4Gene_9.16.Rdata")
metadata$ID<-substr(metadata$TCGA_id,1,16)
colnames(clin_data1)[1]<-"ID"
merge_data<-merge(metadata,clin_data1,by="ID")
clindata<-merge_data[,c(2,3,5,34:40)]
colnames(clindata)[2]<-"group"
colnames(clindata)[1]<-"Tumor_Sample_Barcode"
clindata$bcr_patient_barcode<-substr(clindata$Tumor_Sample_Barcode,1,12)
clindata$type<-"HNSC"

save(clindata,file = "1/data/clindata.Rdata")

maf_df = snv_data@data #提取data数据框
#save(laml,maf_df,file = "maf.Rdata")
phe<-clindata[,c(11,12,2:10)]

phe$new_stage = ifelse(phe$clin_T_stage %in% c("T1"),'s1',
                       ifelse(phe$clin_T_stage %in% c("T2"),'s2',
                              ifelse(phe$clin_T_stage %in% c("T3"),'s3',
                                     ifelse(phe$clin_T_stage %in% c("T4"),'s4','other'
                                     ) ) ) )

table(phe$new_stage)
save(phe,file = "1/data/clin_maf.Rdata")
table(clindata$type)
cp_list=split(phe$bcr_patient_barcode,phe$type)
cg_tmp=snv_data@data
lapply(1:length(cp_list), function(i){
  cp=cp_list[[i]]
  pro=names(cp_list)[i]
  cg_cp_tmp=cg_tmp[substring(cg_tmp$Tumor_Sample_Barcode,1,12) %in% cp,]
  laml = read.maf(maf = cg_cp_tmp)
  save(laml,file= paste0('maftools-',pro,'.Rdata') ) 
})
load(file = 'maftools-HNSC.Rdata') 
oncoplot(maf = laml, top = 10)
#加上临床信息
laml@clinical.data
pos=match(substring(as.character(laml@clinical.data$Tumor_Sample_Barcode),1,12),phe$bcr_patient_barcode)
laml@clinical.data$stage=phe$new_stage[pos]
table(laml@clinical.data$stage)
#other s1 s2 s3 s4 
#24 173 590 220 19 
laml@clinical.data$group=phe$group[pos]
table(laml@clinical.data$group)
laml@clinical.data$gender=phe$gender[pos]
table(laml@clinical.data$gender)
laml@clinical.data$age=phe$age[pos]
table(laml@clinical.data$age)
save(laml,file = "1/data/laml_clin.Rdata")
maf_TMB = tmb(maf = laml,
          captureSize = 50,
          logScale = TRUE)
save(maf_TMB,file = "1/data/maf_TMB.Rdata")

meta<-metadata
meta$sid<-substr(meta$TCGA_id,1,12)
tmb_meta<-merge(meta,maf_TMB,by.x="sid",by.y="Tumor_Sample_Barcode")
TMB_data<-tmb_meta[,c(2,3,7)]
data2 <- TMB_data %>% 
  pivot_longer(cols=-c(1,2),
               names_to= "immunity",
               values_to = "score")

### 高低对比作图
library(ggplot2)
library(ggpubr)
ggboxplot(data2,x="immunity",y="score",fill = "sample",
          ylab = "immunity score",
          xlab = "",
          palette = c("#d7191c","#2b83ba"),
          label.select = list(top.up = 5000, top.down = -5000))+
  rotate_x_text(45)+
  stat_compare_means(mapping = aes(group=sample),
                     label = "p.signif",
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")),
                     method = "wilcox.test")
#以上为ns
#####################################################################
#######Timer数据库免疫细胞浸润
rm(list = ls())
ic<-read.csv("R_data/infiltration_estimation_for_tcga.csv",stringsAsFactors = F)
load(file = "R_data/metadata_tumor_riskgroup_0.4_bestcutpoint.Rdata")
metadata$cell_type<-substr(metadata$TCGA_id,1,15)
ic<-ic[,c(1,grep("TIMER",colnames(ic)))]
ic<-merge(metadata,ic,by="cell_type")
save(ic,file = "ic_TIMER_heatmap.Rdata")


###########################################################
##pheatmap画图
rm(list = ls())
#加载各部分免疫细胞表达
load(file = "output/0.7/cibersort_immunecell_heatmap.Rdata")
cibersort<-cell2
load(file = "output/0.7/ic_TIMER_heatmap.Rdata")
TIMER<-ic
load(file = "output/0.7/res_epic_heatmap.Rdata")
EPIC<-res_epic
load(file = "output/0.7/res_mcp_counter_heatmap.Rdata")
mcp_counter<-res_mcp_counter
load(file = "output/0.7/res_quantiseq_heatmap.Rdata")
QUANTISEQ<-res_quantiseq
load(file = "output/0.7/xCELL_immunecell_heatmap.Rdata")
xCELL<-xCELL_immunecell


cibersort<-cibersort[,-c(17,18,21,22)]#1
colnames(cibersort)<-paste(colnames(cibersort),"CIBERSORT",sep="_")
rownames(TIMER)<-TIMER[,2]
TIMER<-TIMER[,-c(1:3)]#2
xCELL<-as.data.frame(t(xCELL))
colnames(xCELL)<-paste(colnames(xCELL),"xCELL",sep="_")#3
QUANTISEQ<-as.data.frame(QUANTISEQ)
rownames(QUANTISEQ)<-QUANTISEQ[,1]
QUANTISEQ<-QUANTISEQ[,-1]
QUANTISEQ<-as.data.frame(t(QUANTISEQ))
colnames(QUANTISEQ)<-paste(colnames(QUANTISEQ),"QUANTISEQ",sep="_")#4
EPIC<-as.data.frame(EPIC)
rownames(EPIC)<-EPIC[,1]
EPIC<-EPIC[,-1]
EPIC<-as.data.frame(t(EPIC))
colnames(EPIC)<-paste(colnames(EPIC),"EPIC",sep="_")#5
mcp_counter<-as.data.frame(mcp_counter)
rownames(mcp_counter)<-mcp_counter[,1]
mcp_counter<-mcp_counter[,-1]
mcp_counter<-as.data.frame(t(mcp_counter))
colnames(mcp_counter)<-paste(colnames(mcp_counter),"MCP_COUNTER",sep="_")
mcp_counter[mcp_counter=="-Inf"] = 0#6

cibersort<-rownames_to_column(cibersort)
EPIC<-rownames_to_column(EPIC)
QUANTISEQ<-rownames_to_column(QUANTISEQ)
TIMER<-rownames_to_column(TIMER)
xCELL<-rownames_to_column(xCELL)
mcp_counter<-rownames_to_column(mcp_counter)
dat1<-merge(cibersort,EPIC,by="rowname")
dat2<-merge(dat1,QUANTISEQ,by="rowname")
dat3<-merge(dat2,TIMER,by="rowname")
dat4<-merge(dat3,xCELL,by="rowname")
dat5<-merge(dat4,mcp_counter,by="rowname")
data6<-dat5
data6<-data6[,-c(11,14,24,27,31,34,38,44,57,67,68,74,77,80,81)]
data6<-data6[,c(1:18,20:23,64,67,65,66,61:63,24:49,54:60)]
rownames(data6)<-data6[,1]
data6<-data6[,-1]
###进行标准化处理
#install.packages("vegan")
library(vegan)
data7<-as.data.frame(t(data6))
data7_nom<- decostand(data7,"standardize",MARGIN = 1)
data7_nom1<-as.data.frame(scale(data6))
heatdata_1<-as.data.frame(t(data7_nom1))
save(heatdata,file = "output/0.7/pheatdata_immunecells.Rdata")
###制作注释信息
load(file = "R_data/HNSC_cox_riskscore_total_TCGA_0.4.Rdata")
col<-data.sur2[,c(2,5)]
col<-rownames_to_column(col)
colnames(col)<-c("ID","riskscore","group")
row<-as.data.frame(rownames(heatdata_1))
meta<-as.data.frame(colnames(heatdata))
meta$ID<-substr(meta$`colnames(heatdata)`,1,15)
colmeta<-merge(col,meta,by="ID")
colmeta<-colmeta[,c(4,2,3)]
annotation_col<-colmeta[order(colmeta$riskscore),]
rownames(annotation_col)<-annotation_col[,1]
annotation_col<-annotation_col[,-1]
save(annotation_col,file = "output/0.7/annotation_col_data.Rdata")
row$methods<-c(rep("CIBERSORT",16),rep("EPIC",5),rep("QUANTISEQ",8),rep("TIMER",5),rep("xCELL",20))
rownames(row)<-row[,1]
row<-row[-1]
annotation_row<-row
save(annotation_row,file = "output/0.7/annotation_row_data.Rdata")

library(pheatmap)
annotation_col1<-rownames_to_column(annotation_col)
colnames(annotation_col1)[1]<-"TCGA_id"
heatdata1<-as.data.frame(t(heatdata_1))
heatdata1<-rownames_to_column(heatdata1)
colnames(heatdata1)[1]<-"TCGA_id"
heatdata2<-merge(annotation_col1,heatdata1,by="TCGA_id")
heatdata2<-heatdata2[order(heatdata2$riskscore),]
rownames(heatdata2)<-heatdata2[,1]
heatdata3<-heatdata2[,-c(1:3)]
heatdata4<-as.data.frame(t(heatdata3))

pheatmap(heatdata4, annotation_col = annotation_col,
         annotation_row = annotation_row, 
         cluster_rows = FALSE, cluster_cols = FALSE, 
         gaps_row = c(16,21,28,36,41), gaps_col = c(246),
         cellwidth = 1.5, cellheight = 5,
         color = colorRampPalette(colors = c("blue","white","red"))(100),
         show_colnames = F
         )

####################################################################
###immunity scores
rm(list = ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(reshape2)
load(file = "1/data/boxpolt_data.Rdata")
boxpolt_data<-boxpolt_data[,-2]
data_1<-boxpolt_data[,c(1,2,6)]#[,c(1,2,7:9)]
class(data_1$immune.xcell.score)
class(data_1$TumorPurity)
data_1$ESTIMATEScore<-as.numeric(data_1$ESTIMATEScore)
data_1$TumorPurity<-as.numeric(data_1$TumorPurity)
data_1$ImmuneScore_xcell<-as.numeric(data_1$ImmuneScore_xcell)
data_1$StromaScore<-as.numeric(data_1$StromaScore)
data_1$StromalScore<-as.numeric(data_1$StromalScore)
data_1$ImmuneScore<-as.numeric(data_1$ImmuneScore)
data2 <- data_1 %>% 
  pivot_longer(cols=-c(1,2),
               names_to= "immunity",
               values_to = "score")

### 高低对比作图

ggboxplot(data2,x="immunity",y="score",fill = "sample",
          ylab = "immunity score",
          xlab = "",
          palette = c("#d7191c","#2b83ba"),
          label.select = list(top.up = 5000, top.down = -5000))+
  rotate_x_text(45)+
  stat_compare_means(mapping = aes(group=sample),
                     label = "p.signif",
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")),
                     method = "t.test")

##############
###checkpoints
rm(list = ls())
load(file = "E:/0/1/2/3/TCGA/clinical/HNSC_exp_clin_tumor_500.Rdata")

test1<-exprSet[1:153,1:20]
rownames(exprSet)<-exprSet[,1]
exprSet<-exprSet[,-1]
test<-exprSet[1:10,1:10]
library(tibble)
library(dplyr)
library(tidyr)

####所以在此处添加上了sample
checkpoints<-exprSet[,c( "TCGA_id","TIGIT","CD40LG","HHLA2","TNFRSF13B","TNFRSF13C",
                          "CD8A","CD8B","PDCD1","CD276","TGFB1","TGFBR1","CD274","PDCD1LG2")]
load(file = "1/metadata_tumor_4Gene_9.16.Rdata")
me<-metadata
me$ID<-substr(me$TCGA_id,1,15)
checkpoints<-exprSet[,c( "sample","CCL5","CXCL9","CXCL10","GZMA","GZMB","PRF1",
                         "CD8A","CD8B","PDCD1","CD274","PDCD1LG2","CTLA4","HAVCR2","IDO1","LAG3","TIGIT")]

chec<-rownames_to_column(checkpoints)
chec<-checkpoints
colnames(chec)[1]<-"ID"
immun_check<-merge(me,chec,by="ID")
#rownames(immun_check)<-immun_check[,2]
immun_check<-immun_check[,-1]
#colnames(immun_check)[1]<-"sample"
save(data,file = "1/data/HNSC_immune_checkpoints_expdata.Rdata")
###
data2 <- immun_check %>% 
  pivot_longer(cols=-c(1,2),
               names_to= "immune_checkpoints",
               values_to = "expression")
### 高低对比作图

ggboxplot(data2,x="immune_checkpoints",y="expression",fill = "sample",
          ylab = "expression",
          xlab = "",
          palette = c("#d7191c","#2b83ba"),
          label.select = list(top.up = 5000, top.down = -5000))+
  rotate_x_text(45)+
  stat_compare_means(mapping = aes(group=sample),
                     label = "p.signif",
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")),
                     method = "wilcox.test")
#######################################################################
###immune cell expression boxplot-cibersort
rm(list = ls())
load(file = "R_data/exprSet_protein_coding_gene_exprssion_tpm_499.Rdata")
load(file = "1/data/boxpolt_data.Rdata")
tpm_tumor<-as.data.frame(t(sample))
test<-tpm_tumor[1:10,1:10]

sample2<-boxpolt_data[,c(1:3)]


source("CIBERSORT.R")
logtpm<-log(tpm_tumor)
test<-tpm_tumor[1:10,1:10]
test1<-logtpm[1:10,1:10]
exp2 = as.data.frame(logtpm)
exp2 = rownames_to_column(exp2)
test<-exp2[1:10,1:10]
write.table(exp2,file = "exp_4G_9.16.txt",row.names = F,quote = F,sep = "\t")
if(F){TME.results = CIBERSORT("LM22.txt", "exp_4G_9.16.txt" , 
                              perm = 1000, 
                              QN = T)
save(TME.results,file = "ciber_HNSC_4g_9.16.Rdata")}
TME.results<-CIBERSORT("LM22.txt","exp_4G_9.16.txt",perm = 1000, QN = T)
save(TME.results,file = "1/data/TME_results_cibersort.Rdata")
#load("R_data/ciber_HNSC.Rdata")
# load(file = "output/0.7/TME_results_cibersort.Rdata")
TME.results[1:4,1:4]
re <- as.matrix(t(TME.results))
cell = re[,rownames(sample)]#免疫
cell <- cell[-c(23:25),]
cell2<-as.data.frame(t(cell))
save(cell2,file = "1/data/cibersort_immunecell_heatmap.Rdata")
#箱式图#
cell2<-rownames_to_column(cell2)
cell3<-merge(cell2,sample2,by.x="rowname",by.y="ID")
cell3<-cell3[,-24]
data2 <- cell3 %>% 
  pivot_longer(cols=-c(1,24),
               names_to= "IMMU_cell",
               values_to = "score")
cell4<-cell3[,c(1,2,4,5,9,24)]
data2 <- cell4 %>% 
  pivot_longer(cols=-c(1,6),
               names_to= "IMMU_cell",
               values_to = "score")
ggboxplot(data2,x="IMMU_cell",y="score",fill = "sample",
          ylab = "immunity score",
          xlab = "",
          palette = c("#d7191c", "#2b83ba"),
          label.select = list(top.up = 5000, top.down = -5000))+
  rotate_x_text(45)+
  stat_compare_means(mapping = aes(group=sample),
                     label = "p.signif",
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")),
                     method = "wilcox.test")


####
library(tidyverse)
colour <- c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00","#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")

my_theme <- function(){
  theme(panel.grid = element_blank(),       # 网格线
        panel.border = element_blank(),     # 面板边框
        legend.position="right",            # legend位置
        legend.text = element_text(size=8), # legend内容大小
        legend.title = element_text(size=8),# legend标题大小
        axis.line = element_line(size=1),   # 坐标轴线
        text = element_text(family="Times"),# 文本字体
        axis.text.y = element_text(size = 8,face='bold',color='black'),# y轴标签样式
        axis.text.x = element_text(size = 8,face='bold',color='black',angle=90,hjust=1),        # x轴标签样式，angle=45 倾斜 45 度
        axis.title = element_text(size=10,face="bold"),  # 轴标题
        plot.title = element_text(hjust=0.5,size=10))    # 距，设置绘图区域距离边的据类，上、右、下、左
}  

p1 <- TME.results[,1:22] %>% reshape2::melt() %>%
  ggplot(aes(x=Var1,y=value,fill=Var2)) +
  geom_bar(stat='identity') +
  coord_flip()  +
  scale_fill_manual(values =colour ) +
  theme_bw()+ theme(panel.border = element_blank(),panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+my_theme()

pdf("cibersort.pdf")
p1
dev.off()

###################################################################
###相关性作图
rm(list = ls())
#load(file = "R_data/boxpolt_data.Rdata")
load(file = "1/HNSC_cox_riskscore_total_TCGA_4g.Rdata")
#boxpolt_data$ID<-substr(boxpolt_data$TCGA_id,1,15)
corrdata<-merge(boxpolt_data,data.sur2,by.x = "TCGA_id",by.y ="ID" )
corrdata2<-corrdata[,c(1,7:9,12,17)]
rownames(corrdata2)<-corrdata2[,1]
corrdata2<-corrdata2[,-1]
#corrdata2<-boxpolt_data[,c(3,4,6,7,10:12,15,16,2)]
corrdata2$score<-as.numeric(corrdata2$score)
corrdata2$ImmuneScore_xcell<-as.numeric(corrdata2$ImmuneScore_xcell)
corrdata2$StromaScore<-as.numeric(corrdata2$StromaScore)
corrdata2$MicroenvironmentScore<-as.numeric(corrdata2$MicroenvironmentScore)
corrdata2$Lymphocyoe_infiltration<-as.numeric(corrdata2$Lymphocyoe_infiltration)
corrdata2$ImmuneScore<-as.numeric(corrdata2$ImmuneScore)
corrdata2$ESTIMATEScore<-as.numeric(corrdata2$ESTIMATEScore)
corrdata2$TumorPurity<-as.numeric(corrdata2$TumorPurity)
corrdata2<-corrdata2[,c(5,1:4)]
scores_cor<-cor(corrdata2,method = 'spearman')
library(corrplot)
col<- colorRampPalette(c("#1C86EE", "white", "#FF6347"))(20)
corrplot(scores_cor, type="upper", order="hclust", col=col)
corrplot(scores_cor, type="upper", order="hclust", tl.col="black", tl.srt=45)
#如果相关性是非显著的不想显示或用不同的符号显示要怎么办呢？
#首先我们要得到一张相应的pvalue的表
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
# matrix of the p-value of the correlation
p.mat <- cor.mtest(corrdata2)
#不显著的显示叉叉,筛选的标准为0.01
corrplot(scores_cor, type="lower", tl.srt=45,
         p.mat = p.mat, sig.level = 0.05,tl.col="black",
         #insig = "blank"
         col=col)
#不显著的为空白,筛选的标准为0.01
corrplot(scores_cor, type="upper", #order="hclust",
         tl.srt=45,#修改字体 
         p.mat = p.mat, sig.level = 0.05,
         tl.col="black",
         #insig = "blank"
         col=col)


#####################################################################
###
phe<-clindata[,c(11,12,2:10)]

phe$new_stage = ifelse(phe$clin_T_stage %in% c("T1"),'s1',
                       ifelse(phe$clin_T_stage %in% c("T2"),'s2',
                              ifelse(phe$clin_T_stage %in% c("T3"),'s3',
                                     ifelse(phe$clin_T_stage %in% c("T4"),'s4','other'
                                     ) ) ) )

table(phe$new_stage)
save(phe,file = "clinical/clin_maf.Rdata")
#load maf data
load(file = "R_data/HNSC_snv_data.Rdata")
snv_data<-HNSC_snv_data
table(clindata$type)
cp_list=split(phe$bcr_patient_barcode,phe$type)
cg_tmp=HNSC_snv_data
lapply(1:length(cp_list), function(i){
  cp=cp_list[[i]]
  pro=names(cp_list)[i]
  cg_cp_tmp=cg_tmp[substring(cg_tmp$Tumor_Sample_Barcode,1,12) %in% cp,]
  laml = read.maf(maf = cg_cp_tmp)
  save(laml,file= paste0('maftools-',pro,'.Rdata') ) 
})

load(file = 'maftools-HNSC.Rdata') 
oncoplot(maf = laml, top = 10)
#加上临床信息
laml@clinical.data
pos=match(substring(as.character(laml@clinical.data$Tumor_Sample_Barcode),1,12),phe$bcr_patient_barcode)
laml@clinical.data$stage=phe$new_stage[pos]
table(laml@clinical.data$stage)
#other s1 s2 s3 s4 
#24 173 590 220 19 
laml@clinical.data$group=phe$group[pos]
table(laml@clinical.data$group)
laml@clinical.data$gender=phe$gender[pos]
table(laml@clinical.data$gender)
laml@clinical.data$age=phe$age[pos]
table(laml@clinical.data$age)
save(laml,file = "R_data/laml_clin.Rdata")
oncoplot(maf = laml,
         clinicalFeatures=c("group","gender","age","stage"),sortByAnnotation=T,
         top = 20)
library(export)


cp_list1=split(phe$bcr_patient_barcode,phe$group)
cg_tmp=snv_data@data
lapply(1:length(cp_list1), function(i){
  cp1=cp_list1[[i]]
  pro1=names(cp_list1)[i]
  cg_cp_tmp1=cg_tmp[substring(cg_tmp$Tumor_Sample_Barcode,1,12) %in% cp1,]
  laml1 = read.maf(maf = cg_cp_tmp1)
  save(laml1,file= paste0('maftools-',pro1,'.Rdata') ) 
})



#######High Low分开作图
library(maftools)
load(file = 'maftools-High.Rdata') 
oncoplot(maf = laml1, top = 10)
#加上临床信息
laml1@clinical.data
pos1=match(substring(as.character(laml1@clinical.data$Tumor_Sample_Barcode),1,12),
           phe$bcr_patient_barcode[
             laml1@clinical.data$Tumor_Sample_Barcode %in%
               phe$bcr_patient_barcode[phe$group=='High']])
laml1@clinical.data$stage=phe$new_stage[pos1]
table(laml1@clinical.data$stage)
#other s1 s2 s3 s4 
#24 173 590 220 19 
laml1@clinical.data$group=phe$group[pos1]
table(laml1@clinical.data$group)
laml1@clinical.data$gender=phe$gender[pos1]
table(laml1@clinical.data$gender)
laml1@clinical.data$age=phe$age[pos1]
table(laml1@clinical.data$age)
save(laml1,file = "R_data/laml_clin_High.Rdata")
load(file = "R_data/laml_clin_High.Rdata")
oncoplot(maf = laml1,
         #clinicalFeatures=c("age","gender","stage"),
         sortByAnnotation=T,
         top = 20)
plotmafSummary(maf = laml1, rmOutlier = TRUE, addStat = 'median')
library(export)


#############计算TMB
library(ggplot2)
library(ggpubr)
#load(file = 'maftools-HNSC.Rdata') 
#maf_tmb<-tmb(maf=laml,
             #captureSize = 50,
             #logScale = TRUE)
#save(maf_tmb,file = "output/0.7/maf_TMB.Rdata")
####箱线图
load(file = "R_data/metadata_tumor_riskgroup_0.4_bestcutpoint.Rdata")
#metadata$Tumor_Sample_Barcode<-substr(metadata$TCGA_id,1,12)

####
#BiocManager::install("PoisonAlien/TCGAmutations")
# TCGAmutations包整合了TCGA中全部样本的maf文件
#devtools::install_github(repo = "PoisonAlien/TCGAmutations")
library(TCGAmutations)
tmp=as.data.frame(tcga_available())

dt <- TCGAmutations::tcga_load(study = "HNSC")

dt <- dt@data

dt1 <- as.data.frame( table(dt$Tumor_Sample_Barcode))
dt1$ID<-substr(dt1$Var1,1,16)
metadata$ID<-substr(metadata$TCGA_id,1,16)
dt2<-merge(dt1,metadata,by="ID")
dt2<-dt2[,2:3]
names(dt2) <- c('Barcode', 'Freq')

dt2$tmb <- dt2$Freq/38

names(dt2)
write.csv(dt2, file = 'output/0.7/HNSC_TMB_TCGAmutations.csv')
dt2$ID<-substr(dt2$Barcode,1,16)

box_data<-merge(metadata,dt2,by="ID")
data_1<-box_data[,c(2,3,6)]
class(data_1$tmb)

data2 <- data_1 %>% 
  pivot_longer(cols=-c(1,2),
               names_to= "TMB",
               values_to = "score")

### 肿瘤对比正常作图

ggboxplot(data2,x="TMB",y="score",fill = "sample",
          ylab = "TMB",
          xlab = "",
          palette = c("#d7191c","#2b83ba"),
          outlier.shape = 19
          #label.select = list(top.up = 5000, top.down = -5000)
)+
  rotate_x_text(45)+
  stat_compare_means(mapping = aes(group=sample),
                     label = "p.signif",
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")),
                     method = "wilcox.test")



##################################################################
##不同算法免疫细胞表达量
rm(list = ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(reshape2)
load(file = "output/0.7/xCELL_hetadata.Rdata")
load("R_data/metadata_tumor_riskgroup_0.4_bestcutpoint.Rdata")
rownames(xCELL)<-xCELL[,1]
mer<-merge(xCELL,metadata,by.x = "rowname",by.y = "TCGA_id")
data_1<-mer[,c(1,2,12,20,21,24,3:11,25:28,23,22,16:19,13:15,29)]
rownames(data_1)<-data_1[,1]

data_1$ImmuneScore<-as.numeric(data_1$ImmuneScore)
data_1$StromalScore<-as.numeric(data_1$StromalScore)
data_1$ESTIMATEScore<-as.numeric(data_1$ESTIMATEScore)
data2 <- data_1 %>% 
  pivot_longer(cols=-c(1,29),
               names_to= "immunity",
               values_to = "score")

### 肿瘤对比正常作图

ggboxplot(data2,x="immunity",y="score",fill = "sample",
          ylab = "immunity score",
          xlab = "",
          palette = c("#d7191c","#2b83ba"),
          label.select = list(top.up = 5000, top.down = -5000))+
  rotate_x_text(45)+
  stat_compare_means(mapping = aes(group=sample),
                     label = "p.signif",
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")),
                     method = "wilcox.test")

rm(list = ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(reshape2)
load(file = "output/0.7/TIMER_hetadata.Rdata")
load("R_data/metadata_tumor_riskgroup_0.4_bestcutpoint.Rdata")
rownames(TIMER)<-TIMER[,1]
mer<-merge(TIMER,metadata,by.x = "rowname",by.y = "TCGA_id")
data_1<-mer
rownames(data_1)<-data_1[,1]

data_1$ImmuneScore<-as.numeric(data_1$ImmuneScore)
data_1$StromalScore<-as.numeric(data_1$StromalScore)
data_1$ESTIMATEScore<-as.numeric(data_1$ESTIMATEScore)
data2 <- data_1 %>% 
  pivot_longer(cols=-c(1,8),
               names_to= "immunity",
               values_to = "score")

### 肿瘤对比正常作图

ggboxplot(data2,x="immunity",y="score",fill = "sample",
          ylab = "immunity score",
          xlab = "",
          palette = c("#d7191c","#2b83ba"),
          label.select = list(top.up = 5000, top.down = -5000))+
  rotate_x_text(45)+
  stat_compare_means(mapping = aes(group=sample),
                     label = "p.signif",
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")),
                     method = "wilcox.test")

rm(list = ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(reshape2)
load(file = "output/0.7/mcp_counter_hetadata.Rdata")
load("R_data/metadata_tumor_riskgroup_0.4_bestcutpoint.Rdata")
rownames(mcp_counter)<-mcp_counter[,1]
mer<-merge(mcp_counter,metadata,by.x = "rowname",by.y = "TCGA_id")
data_1<-mer
rownames(data_1)<-data_1[,1]

data_1$ImmuneScore<-as.numeric(data_1$ImmuneScore)
data_1$StromalScore<-as.numeric(data_1$StromalScore)
data_1$ESTIMATEScore<-as.numeric(data_1$ESTIMATEScore)
data2 <- data_1 %>% 
  pivot_longer(cols=-c(1,13),
               names_to= "immunity",
               values_to = "score")

### 肿瘤对比正常作图

ggboxplot(data2,x="immunity",y="score",fill = "sample",
          ylab = "immunity score",
          xlab = "",
          palette = c("#d7191c","#2b83ba"),
          label.select = list(top.up = 5000, top.down = -5000))+
  rotate_x_text(45)+
  stat_compare_means(mapping = aes(group=sample),
                     label = "p.signif",
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")),
                     method = "wilcox.test")

rm(list = ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(reshape2)
load(file = "output/0.7/QUANTISEQ_hetadata.Rdata")
load("R_data/metadata_tumor_riskgroup_0.4_bestcutpoint.Rdata")
rownames(QUANTISEQ)<-QUANTISEQ[,1]
mer<-merge(QUANTISEQ,metadata,by.x = "rowname",by.y = "TCGA_id")
data_1<-mer
rownames(data_1)<-data_1[,1]

data_1$ImmuneScore<-as.numeric(data_1$ImmuneScore)
data_1$StromalScore<-as.numeric(data_1$StromalScore)
data_1$ESTIMATEScore<-as.numeric(data_1$ESTIMATEScore)
data2 <- data_1 %>% 
  pivot_longer(cols=-c(1,13),
               names_to= "immunity",
               values_to = "score")

### 肿瘤对比正常作图

ggboxplot(data2,x="immunity",y="score",fill = "sample",
          ylab = "immunity score",
          xlab = "",
          palette = c("#d7191c","#2b83ba"),
          label.select = list(top.up = 5000, top.down = -5000))+
  rotate_x_text(45)+
  stat_compare_means(mapping = aes(group=sample),
                     label = "p.signif",
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")),
                     method = "wilcox.test")


###X-cell score
rm(list = ls())
library(xCell)
library(dplyr)
library(tibble)
library(tidyverse)
library(tidyr)
###############制作xcell需要的数据
##加载tpm数据(9.527与之前一样)
load(file = "R_data/exprSet_protein_coding_gene_exprssion_tpm_499.Rdata")
sample_input_xcell<-as.data.frame(t(sample))
test<-sample_input_xcell[1:10,1:10]
save(sample_input_xcell,file = "R_data/sample_input_xcell.Rdata")

#####
load(file = "R_data/sample_input_xcell.Rdata")
exprSet<-sample_input_xcell
xCell_RNAseq_score<- xCellAnalysis(exprSet,rnaseq = T)
xCell_RNAseq_score<-as.data.frame(xCell_RNAseq_score)
save(xCell_RNAseq_score,file = "1/data/xCell_RNAseq_score.Rdata")
immune_xcell_score<-xCell_RNAseq_score[c(65,66,67),]
immune_xcell_score<-as.data.frame(t(immune_xcell_score))
save(immune_xcell_score,file = "1/data/immune_xcell_score.Rdata")
xCell_RNAseq_score1<-rownames_to_column(xCell_RNAseq_score)
xCELL_immunecell<-xCell_RNAseq_score[c(4,6:14,17,19,21,28,32:35,38,46,49,50,57,61:64),]
xCELL_immunecell<-xCell_RNAseq_score[-c(65:67),]
save(xCELL_immunecell,file = "1/data/xCELL_immunecells.Rdata")


library('dplyr')
library('ggplot2')
library('tidyr')
library('immunedeconv')
#devtools::install_github("GfellerLab/EPIC", build_vignettes=TRUE)
#install.packages("testit")
#install.packages("data.tree")
#install.packages("limSolve")
#install.packages("ComICS")
#install.packages("pracma")

#install_github("ebecht/MCPcounter",ref="master", subdir="Source")
#install.packages("remotes")
#library(remotes)
#remotes::install_github("omnideconv/immunedeconv")
library('tibble')
load(file = "R_data/exprSet_protein_coding_gene_exprssion_tpm_499.Rdata")
tpm_tumor<-as.data.frame(t(sample))
test<-tpm_tumor[1:10,1:10]
logtpm<-log(tpm_tumor)
deconvolution_methods
res_quantiseq = deconvolute(logtpm, "quantiseq", tumor = TRUE)
save(res_quantiseq,file = "1/data/quantiseq_immune_cell.Rdata")
#res_mcp_counter = deconvolute(logtpm, "mcp_counter",tumor = TRUE)
#save(res_mcp_counter,file = "output/0.7/res_mcp_counter_heatmap.Rdata")
res_epic = deconvolute(logtpm, "epic",tumor = TRUE)
save(res_epic,file = "1/data/epic_immune_cell.Rdata")
res_consensus_tme = deconvolute(logtpm, "consensus_tme",tumor = TRUE)
save(res_epic,file = "1/data/epic_immune_cell.Rdata")

library(MCPcounter)
genes <- data.table::fread("resource/MCPcounter-master/Signatures/genes.txt",data.table = F)
probesets <- data.table::fread("resource/MCPcounter-master/Signatures/probesets.txt",data.table = F,header = F)
results<- MCPcounter.estimate(logtpm,
                              featuresType= "HUGO_symbols", 
                              probesets=probesets,
                              genes=genes)
##清洗数据

save(results,file = "1/data/results_immune_cell.Rdata")

###作图对比
###immunity scores
rm(list = ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(tidyverse)

data<-as.data.frame(t(xCELL_immunecell))
identical(metadata$TCGA_id,rownames(data))
data_1<-cbind(metadata,data)

data2 <- data_1 %>% 
  pivot_longer(cols=-c(1,2),
               names_to= "immunity_cell",
               values_to = "score")

### 高低对比作图

ggboxplot(data2,x="immunity_cell",y="score",fill = "sample",
          ylab = "immunity score",
          xlab = "",
          palette = c("#d7191c","#2b83ba"),
          label.select = list(top.up = 5000, top.down = -5000))+
  rotate_x_text(45)+
  stat_compare_means(mapping = aes(group=sample),
                     label = "p.signif",
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")),
                     method = "t.test")

##quantiseq

data<-res_quantiseq %>%
  column_to_rownames("cell_type") 
data<-as.data.frame(t(data))
identical(metadata$TCGA_id,rownames(data))
data_1<-cbind(metadata,data)

data2 <- data_1 %>% 
  pivot_longer(cols=-c(1,2),
               names_to= "immunity_cell",
               values_to = "score")

### 高低对比作图

ggboxplot(data2,x="immunity_cell",y="score",fill = "sample",
          ylab = "quantiseq immunity score",
          xlab = "",
          palette = c("#d7191c","#2b83ba"),
          label.select = list(top.up = 5000, top.down = -5000))+
  rotate_x_text(45)+
  stat_compare_means(mapping = aes(group=sample),
                     label = "p.signif",
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")),
                     method = "t.test")


data<-res_epic %>%
  column_to_rownames("cell_type") 
data<-data[-8,]
data<-as.data.frame(t(data))
identical(metadata$TCGA_id,rownames(data))
data_1<-cbind(metadata,data)

data2 <- data_1 %>% 
  pivot_longer(cols=-c(1,2),
               names_to= "immunity_cell",
               values_to = "score")

### 高低对比作图

ggboxplot(data2,x="immunity_cell",y="score",fill = "sample",
          ylab = "epic immunity score",
          xlab = "",
          palette = c("#d7191c","#2b83ba"),
          label.select = list(top.up = 5000, top.down = -5000))+
  rotate_x_text(45)+
  stat_compare_means(mapping = aes(group=sample),
                     label = "p.signif",
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")),
                     method = "t.test")

#MCPcounter
load(file = "1/data/results_immune_cell.Rdata")
load(file = "1/metadata_tumor_4Gene_9.16.Rdata")
data<-as.data.frame(t(results))
identical(metadata$TCGA_id,rownames(data))
data_1<-cbind(metadata,data)

data2 <- data_1 %>% 
  pivot_longer(cols=-c(1,2),
               names_to= "immunity_cell",
               values_to = "score")

### 高低对比作图

ggboxplot(data2,x="immunity_cell",y="score",fill = "sample",
          ylab = "immunity score",
          xlab = "",
          palette = c("#d7191c","#2b83ba"),
          label.select = list(top.up = 5000, top.down = -5000))+
  rotate_x_text(45)+
  stat_compare_means(mapping = aes(group=sample),
                     label = "p.signif",
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")),
                     method = "t.test")


