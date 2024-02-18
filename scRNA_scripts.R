setwd("/home/lilixuan/GSE181919-revision/")
rm(list = ls())
library(Seurat)
library(AnnoProbe)
library(GEOquery)
library(infercnv)
library(tidyverse)
library(SingleR)
library(celldex)
library(RColorBrewer)
library(SingleCellExperiment)
library(ggsci)

##1.0 读取数据，并调节格式一致
GSE181919<-read.delim2("/share/share/wanglixuan/SRA/GSE181919/GSE181919_UMI_counts.txt.gz",header = TRUE,sep = "\t",quote = "\"", dec = ",", fill = TRUE, comment.char = "")
sce.meta <- read.table("/share/share/wanglixuan/SRA/GSE181919/GSE181919_Barcode_metadata.txt.gz",header = TRUE)
rownames(sce.meta) <- gsub("-",".",rownames(sce.meta))
dim(GSE181919)
##1.1 创建Seurat文件，供后续分析
scRNA = CreateSeuratObject(counts=GSE181919,
                           meta.data = sce.meta,
                           min.cells = 55, #过滤少于55个细胞（0.1%）检测出的基因（min.cells = 55）
                           min.features = 200) #过滤检测少于200个基因的细胞（min.features = 200）
#save(scRNA,file = "/share/share/wanglixuan/SRA/GSE181919/R_data/GSE181919seurat.Rdata")


scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
table(scRNA$sample.id)
table(scRNA$tissue.type)
VlnPlot(scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
scRNA <- subset(scRNA, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 10)
saveRDS(scRNA,"R_data/scRNA_after_QC.rds")

##去除LP的组织
#0.1 读取之前质控后的数据

scRNA<-readRDS("R_data/scRNA_after_QC.rds")
scRNA1<-subset(scRNA,tissue.type%in%c("CA","LN","NL"))
saveRDS(scRNA1,"R_data/scRNA_noLP.rds")
#0.2
#为了做进一步的分析，我们需要对数据进行归一化（Normalization）和标准化（Z-score）。
#这里我们介绍一下经典的Normalization方法，这个方法假设所有细胞均含有10,000个UMI。
#normalization.method有三种方法：LogNormalize、CLR、RC
#scale.factor可以设置，主要变化在细胞亚群的数量上
scRNA<-readRDS("R_data/scRNA_noLP.rds")
scRNA<-NormalizeData(scRNA,normalization.method = "LogNormalize",scale.factor = 10000)
###寻找高变基因
#寻找高变基因的方法：vst、mean.var.plot 和 dispersion（找到的高变基因各不相同）
scRNA<-FindVariableFeatures(scRNA,selection.method = "vst",nfeatures = 2000)
top10 <- head(VariableFeatures(scRNA), 10)
top10 
plot1 <- VariableFeaturePlot(scRNA)
LabelPoints(plot = plot1, 
            points = top10, 
            repel = TRUE, 
            xnudge = 0, 
            ynudge = 0)
###扩展数据（ScalingData）
#用于PCA降维，一般只需要对高变基因进行Scale（默认使用高变基因）
#all.genes <- rownames(scRNA)
scRNA<-ScaleData(scRNA,vars.to.regress = "percent.mt")#features	要标准化/居中的特征基因，默认是高可变基因;vars.to.regress	要回归的变量， 例如，nUMI 或percent.mito
#0.3 主成分分析
scRNA <- RunPCA(scRNA, features = VariableFeatures(object = scRNA))
## 可视化
VizDimLoadings(scRNA, dims = 1:9, 
               reduction = "pca") + 
  theme(axis.text=element_text(size=5), 
        axis.title=element_text(size=8,face="bold"))
DimPlot(scRNA, reduction = "pca", group.by = "sample.id")
DimPlot(scRNA, reduction = "pca", group.by = "tissue.type")

DimHeatmap(scRNA, dims = 1:8, nfeatures = 20, cells = 500, balanced = T)

#pca量化拐点
pct<-scRNA[["pca"]]@stdev/sum(scRNA[["pca"]]@stdev)*100
cumu <- cumsum(pct)
pct.use<-min(which(cumu > 90 & pct < 5)[1],
             sort(which((pct[1:length(pct) - 1]-pct[2:length(pct)])> 0.1),decreasing = T)[1] + 1)
#结果显示拐点是18

#0.4 聚类
#之前确定了10个PCs最合适，这里我们将PC1:18纳入，进行聚类。
#这一步非常重要，因为后面marker的选择要基于聚类的结果。
saveRDS(scRNA,file = "R_data/SCRNA_unUMAP_240103_wlx.rds")
scRNA<-readRDS("R_data/SCRNA_unUMAP_240103_wlx.rds")
library(reticulate)
use_condaenv(condaenv = "base", required = T)
use_python('/home/lilixuan/mambaforge/bin/python3.10',required = T)
py_config()
umap_learn<-import("umap")
py_module_available("umap-learn")
library(uwot)
scRNA2 <- FindNeighbors(scRNA, dims = 1:30)
scRNA2 <- FindClusters(scRNA2, resolution = 0.2)
scRNA3 <- RunUMAP(scRNA2,dims = 1:50,reduction = "pca",uwot.sgd = FALSE,metric = 'correlation',min.dist = 0.3,n.neighbors = 50L)
#umap.method = 'umap-learn',

saveRDS(scRNA3,"R_data/scRNA_nei30_clu0.2_umap50_0.3_50L_wlx240106.rds")
##以上是在ssh做
#0.5 加载分析后的scRNA，做dimplot，做findmarker找出标志基因对细胞进行分群
#scRNA<-readRDS("scRNA_umap_0.3_25_0.2.rds")

###结果可视化
DimPlot(scRNA3,raster=FALSE,label = T)
DimPlot(scRNA,raster=FALSE,label = T,group.by = "patient.id")

#0.5 寻找marker基因
###寻找marker基因
#setwd("/home/lilixuan/GSE181919-revision/")
library(tidyverse)
library(dplyr)
scRNA.markers <- FindAllMarkers(scRNA3, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
scRNA.markers_top30<-scRNA.markers %>%
  group_by(cluster) %>%
  #top_n(avg_log2FC,n=30)
  slice_max(n = 30, order_by = avg_log2FC)
seq<-rep(c(1:30),times=23)
scRNA.markers_top30<-data_frame(seq,scRNA.markers_top30$cluster,scRNA.markers_top30$gene)
scRNA.markers_top30_width<-pivot_wider(scRNA.markers_top30,names_from = "scRNA.markers_top30$cluster",values_from = "scRNA.markers_top30$gene")
write.csv(scRNA.markers_top30_width, "output/scRNA.markers_top30.csv", row.names = F)

markergene_list<-c("KRT14", "KRT17", "KRT6A", "KRT5", "KRT19", "KRT8","KRT16", "KRT18", "KRT6B", "KRT15", "KRT6C", "KRTCAP3","EPCAM", "SFN",#malignant
                   "FAP", "PDPN", "COL1A2", "DCN", "COL3A1", "COL6A1",#fibroblasts
                   "ACTA1", "ACTN2", "MYL2", "MYH2",#myocytes
                   "CD2", "CD3D", "CD3E", "CD3G", #T cells
                   "SLAMF7", "CD79A", "BLNK", "FCRL5", #B/Plasma cells
                   "CD14", "CD163", "CD68", "FCGR2A", "CSF1R", #macrophages
                   "CD40", "CD80", "CD83", "CCR7", #dendritic cells
                   "CMA1", "MS4A2", "TPSAB1", "TPSB2", #mast cells
                   "PECAM1", "VWF", "ENG" #endothelial cells
)
VlnPlot(scRNA3,features = markergene_list,pt.size = 0,stack = T)+NoLegend()
DimPlot(scRNA,label = T)

##0.6 注释细胞
scRNA3$seurat_clusters<- recode(scRNA3$seurat_clusters,
                               '0'='T_cells','1'='T_cells','12'='T_cells',
                               '3'='Fibroblasts','4'='Fibroblasts','8'='Fibroblasts','14'='Fibroblasts','16'='Fibroblasts',
                               '9'='Epithelial_cells','7'='Epithelial_cells','11'='Epithelial_cells','15'='Epithelial_cells','20'='Epithelial_cells','22'='Epithelial_cells',
                               '2'='Macrophages',
                               '10'='Endothelial_cells','17'='Endothelial_cells',
                               '6'='B/Plasma_cells','21'='B/Plasma_cells','5'='B/Plasma_cells',
                               '13'='Dendritic_cells',
                               '18'='Myocytes',
                               '19'='Mast_cells')
DimPlot(scRNA3,raster=FALSE,label = T,group.by = "seurat_clusters")

saveRDS(scRNA3,"R_data/scRNA_cellidentified_chubu240106.rds")


#0.4 聚类
#之前确定了10个PCs最合适，这里我们将PC1:18纳入，进行聚类。
#这一步非常重要，因为后面marker的选择要基于聚类的结果。
scRNA2 <- RunUMAP(scRNA,dims = 1:50,umap.method = 'umap-learn',reduction = "pca",uwot.sgd = FALSE,metric = 'correlation',min.dist = 0.2,n.neighbors = 30L)
scRNA2 <- FindNeighbors(scRNA2, dims = 1:25)
scRNA2 <- FindClusters(scRNA2, resolution = 0.2)

##以上是在ssh做
#0.5 加载分析后的scRNA，做dimplot，做findmarker找出标志基因对细胞进行分群
scRNA<-readRDS("/share/share/wanglixuan/SRA/GSE181919/Rdata_for submission/scRNA_umap_0.3_25_0.2.rds")

###结果可视化
DimPlot(scRNA,raster=FALSE,label = T)
DimPlot(scRNA,raster=FALSE,label = T,group.by = "patient.id")

#0.5 寻找marker基因
###寻找marker基因
scRNA.markers <- FindAllMarkers(scRNA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
scRNA.markers_top30<-scRNA.markers %>%
  group_by(cluster) %>%
  slice_max(n = 30, order_by = avg_log2FC)
seq<-rep(c(1:30),times=20)
scRNA.markers_top30<-data_frame(seq,scRNA.markers_top30$cluster,scRNA.markers_top30$gene)
scRNA.markers_top30_width<-pivot_wider(scRNA.markers_top30,names_from = "scRNA.markers_top30$cluster",values_from = "scRNA.markers_top30$gene")
write.csv(scRNA.markers_top30_width, "/share/share/wanglixuan/SRA/GSE181919/Rdata_for submission/scRNA.markers_top30.csv", row.names = F)

markergene_list<-c("KRT14", "KRT17", "KRT6A", "KRT5", "KRT19", "KRT8","KRT16", "KRT18", "KRT6B", "KRT15", "KRT6C", "KRTCAP3","EPCAM", "SFN",#malignant
                   "FAP", "PDPN", "COL1A2", "DCN", "COL3A1", "COL6A1",#fibroblasts
                   "ACTA1", "ACTN2", "MYL2", "MYH2",#myocytes
                   "CD2", "CD3D", "CD3E", "CD3G", #T cells
                   "SLAMF7", "CD79A", "BLNK", "FCRL5", #B/Plasma cells
                   "CD14", "CD163", "CD68", "FCGR2A", "CSF1R", #macrophages
                   "CD40", "CD80", "CD83", "CCR7", #dendritic cells
                   "CMA1", "MS4A2", "TPSAB1", "TPSB2", #mast cells
                   "PECAM1", "VWF", "ENG" #endothelial cells
)
VlnPlot(scRNA,features = markergene_list,pt.size = 0,stack = T)+NoLegend()
DimPlot(scRNA,label = T)

##0.6 注释细胞
scRNA<-readRDS("R_data/scRNA_umap_0.3_25_0.2.rds")
DimPlot(object = scRNA,group.by = "cell.type",label = T)
DimPlot(object = scRNA,label = T)
scRNA$seurat_clusters<- recode(scRNA$seurat_clusters,
                               '0'='T_cells','1'='T_cells','12'='T_cells',
                               '3'='Fibroblasts','4'='Fibroblasts','7'='Fibroblasts','15'='Fibroblasts','16'='Fibroblasts',
                               '9'='Epithelial_cells','10'='Epithelial_cells','11'='Epithelial_cells','13'='Epithelial_cells',
                               '2'='Macrophages',
                               '8'='Endothelial_cells',
                               '6'='B/Plasma_cells','19'='B/Plasma_cells','5'='B/Plasma_cells',
                               '14'='Dendritic_cells',
                               '17'='Myocytes',
                               '18'='Mast_cells')
DimPlot(scRNA,raster=FALSE,label = T,group.by = "seurat_clusters")

saveRDS(scRNA,"R_data/scRNA_cellidentified_chubu_240114.rds")


# 0.7 cnv鉴定上皮细胞恶性程度
rm(list = ls())
#setwd("/home/lilixuan/GSE181919-revision/")
scRNA<-readRDS("R_data/scRNA_cellidentified_chubu.rds")
#0.7.1 选出epi和参照细胞类型
levels(scRNA$seurat_clusters)
scRNA<-SetIdent(scRNA,value = "seurat_clusters")
DimPlot(object = scRNA,group.by = "cell.type",label = T)
DimPlot(object = scRNA,label = T)
DimPlot(object = scRNA,group.by = "RNA_snn_res.0.2",label = T)
scRNA_object <- subset(scRNA, seurat_clusters%in%c("Macrophages","Epithelial_cells"))
scRNA_object <- subset(scRNA, seurat_clusters=="Macrophages"&tissue.type%in%c("NL","LN")|seurat_clusters=="Epithelial_cells"&tissue.type%in%c("CA","LN"))
#3.2.2 制作counts信息
counts  <- GetAssayData(scRNA_object,assay = "RNA",slot = "counts")
counts <- as.data.frame(counts)
#3.2.3 细胞类型注释文件。第一列为细胞barcode，第二列为细胞类型名称。seurat对象提取如下
anno <- data.frame(scRNA_object@active.ident)
identical(colnames(counts),rownames(anno))
colnames(anno)<-"cluster"
anno$cluster<-ifelse(anno$scRNA_object.active.ident=="2","Macrophages","Epithelial_cells")
anno<-as.data.frame(anno[,2])

rownames(anno)<-colnames(counts)
colnames(anno)<-"cluster"
#3.2.4 制作基因染色体位置信息
library(AnnoProbe)
library(GEOquery)
geneInfor=annoGene(rownames(counts),"SYMBOL",'human')
geneInfor=geneInfor[with(geneInfor, order(chr, start)),c(1,4:6)]      
geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
length(unique(geneInfor[,1]))
head(geneInfor)
counts=counts[match(geneInfor[,1], rownames(counts)),]    #将表达矩阵的基因排列顺序与geneInfor的基因排列顺序弄成一致
rownames(geneInfor) <- geneInfor$SYMBOL   
geneInfor <- geneInfor[,-1]     #这样我们就制作好了染色体位置信息和排列顺序好的count表达矩阵

#验证1   表达矩阵的列名要与meta的行名一致
identical(colnames(counts),rownames(anno))  
#验证2   表达矩阵的行名要与geneInfor的行名一致
identical(rownames(counts),rownames(geneInfor))
#因此三个输入数据准备好了   dat-表达矩阵  meta-分组信息  geneInfor-基因染色体信息
#3.2.6 创建cnv project
library(infercnv)
infercnv_obj = CreateInfercnvObject(raw_counts_matrix = counts,
                                    annotations_file = anno,
                                    delim="\t",
                                    gene_order_file = geneInfor,
                                    #min_max_counts_per_cell = c(100, +Inf),
                                    ref_group_names = c("Macrophages"))
saveRDS(infercnv_obj,"output/infercnv_obj_240106.rds")

#3.2.7 计算cnv,工作站后台跑
#infercnv_obj1 = infercnv::run(infercnv_obj,cutoff = 0.1,out_dir = "cnv3/", cluster_by_groups = T,HMM = T,denoise = TRUE,write_expr_matrix = T,num_threads = 8) 


#3.2.8 聚类并定义肿瘤细胞
library(ComplexHeatmap)
library(circlize)
library("RColorBrewer")
infercnv_obj_F<-readRDS("/home/lilixuan/GSE181919-revision/output/cnv240114/run.final.infercnv_obj")

expr <- infercnv_obj_F@expr.data
normal_loc <- infercnv_obj_F@reference_grouped_cell_indices
normal_loc <- normal_loc$Macrophages
test_loc <- infercnv_obj_F@observation_grouped_cell_indices
test_loc <- test_loc$Epithelial_cells

anno.df=data.frame(
  CB=c(colnames(expr)[normal_loc],colnames(expr)[test_loc]),
  class=c(rep("normal",length(normal_loc)),rep("test",length(test_loc)))
)
head(anno.df)
#提取ca、ln的epi细胞
#scRNA<-readRDS("/share/share/wanglixuan/SRA/GSE181919/R_data/scRNA_after_QC.rds")
meta_epi<-scRNA@meta.data[,c("tissue.type","sample.id","seurat_clusters")]
meta_epi<-subset(meta_epi,meta_epi$seurat_clusters=="Epithelial_cells")
meta_epi<-meta_epi %>%
  filter(tissue.type !="NL")

meta_epi<-rownames_to_column(meta_epi)
anno.df1<-anno.df %>%
  filter(CB %in% meta_epi$rowname)
gn <- rownames(expr)

sub_geneFile <- geneInfor[intersect(gn,rownames(geneInfor)),]
expr=expr[intersect(gn,rownames(geneInfor)),]
head(sub_geneFile,4)
expr[1:4,1:4]
dim(expr)
#聚类，8类，提取结果
set.seed(20240114)
kmeans.result <- kmeans(t(expr), 8)
kmeans_df <- data.frame(kmeans_class=kmeans.result$cluster)
kmeans_df$CB=rownames(kmeans_df)
kmeans_df=kmeans_df%>%inner_join(anno.df,by="CB") #合并
kmeans_df_s=arrange(kmeans_df,kmeans_class) #排序
rownames(kmeans_df_s)=kmeans_df_s$CB
kmeans_df_s$CB=NULL
kmeans_df_s$kmeans_class=as.factor(kmeans_df_s$kmeans_class) #将kmeans_class转换为因子，作为热图的一个注释，最终在热图上就会按照1:7排序
head(kmeans_df_s)

#定义热图的注释，及配色
top_anno <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "NA",col="NA"), labels = 1:22,labels_gp = gpar(cex = 1.5)))
color_v=RColorBrewer::brewer.pal(8, "Set3")[1:8] #类别数
names(color_v)=as.character(1:8)
left_anno <- rowAnnotation(df = kmeans_df_s,col=list(class=c("test"="red","normal" = "blue"),kmeans_class=color_v))

#下面是绘图(只有-4是可以用的！)
Km_t<-t(expr)[rownames(kmeans_df_s),]
pdf("output/cnv240114/infercnv_wlx240115_1.pdf",width = 15,height = 10)
ht = Heatmap(Km_t, #绘图数据的CB顺序和注释CB顺序保持一致
             col = colorRamp2(c(0,1,2), c("#377EB8","#F0F0F0","#E41A1C")), #如果是10x的数据，这里的刻度会有所变化
             cluster_rows = F,cluster_columns = F,show_column_names = F,show_row_names = F,
             column_split = factor(sub_geneFile$chr, paste("chr",1:22,sep = "")), #这一步可以控制染色体顺序，即使你的基因排序文件顺序是错的
             column_gap = unit(2, "mm"),
             
             heatmap_legend_param = list(title = "Modified expression",direction = "vertical",title_position = "leftcenter-rot",at=c(0,1,2),legend_height = unit(3, "cm")),
             
             top_annotation = top_anno,left_annotation = left_anno, #添加注释
             row_title = NULL,column_title = NULL)
draw(ht, heatmap_legend_side = "right")
dev.off()

expr2=expr-1
expr2=expr2 ^ 2
CNV_score=as.data.frame(colMeans(expr2))
colnames(CNV_score)="CNV_score"
CNV_score$CB=rownames(CNV_score)

kmeans_df_s$CB=rownames(kmeans_df_s)
CNV_score=CNV_score%>%inner_join(kmeans_df_s,by="CB")

CNV_score%>%ggplot(aes(kmeans_class,CNV_score))+geom_violin(aes(fill=kmeans_class),color="NA")+
  scale_fill_manual(values = color_v)+
  theme_bw()

#将class2中的test定义为正常上皮，其他定义为肿瘤
#每一类对应的CB保存在kmeans_df_s数据框中
write.table(kmeans_df_s, file = "output/infercnv_kmeans_df_s.txt", quote = FALSE, sep = '\t', row.names = T, col.names = T)
kmeans_df_s<-read.table("output/infercnv_kmeans_df_s.txt")
anno1<-rownames_to_column(anno)
kmeans_df_s1<-rownames_to_column(kmeans_df_s)
kmeans_df_s1<-kmeans_df_s1[,-4]
kmeans_df_s1<-kmeans_df_s1 %>%
  filter(class=="test")
kmeans_df_s1$cell<-ifelse(kmeans_df_s1$kmeans_class=="2","Epithelial_cells","Malignant_cells")
ma_zhushi<-merge(kmeans_df_s1,meta_epi,by="rowname")

#取EP细胞和其他已定义非肿瘤细胞
ep<-scRNA@meta.data[,c("tissue.type","seurat_clusters")]
ep_zhushi<-subset(ep,ep$seurat_clusters=="Epithelial_cells"&ep$tissue.type=="NL")
other_zhushi<-ep %>%
  filter(!seurat_clusters=="Epithelial_cells")

##将各个细胞重新拼接
ma_zhushi1<-ma_zhushi[,c(1,4)]
colnames(ma_zhushi1)<-c("id","cell_rev")

ep_zhushi1<-ep_zhushi %>%
  rownames_to_column("id") %>%
  select(c("id","seurat_clusters"))
colnames(ep_zhushi1)[2]<-"cell_rev"

other_zhushi1<-other_zhushi %>%
  rownames_to_column("id") %>%
  select(c("id","seurat_clusters"))
colnames(other_zhushi1)[2]<-"cell_rev"

##纵向合并
cell_rev<-rbind(ep_zhushi1,ma_zhushi1,other_zhushi1)

#使行名与seurat一致
ep<-ep%>%
  rownames_to_column("id")

data<-inner_join(ep,cell_rev)
identical(ep$id,data$id)
##添加回seurat
scRNA <- AddMetaData(object = scRNA,                #seurat对象
                    metadata = data$cell_rev,               #需要添加的metadata
                    col.name = "cell_rev")        #给新添加的metadata命名为

saveRDS(scRNA,file = "R_data/scRNA_CELL_IDENTIFY_REV_240115.rds")


#1.0 提取肿瘤和上皮细胞，对比SGOC通路的表达差异
rm(list = ls())
scRNA<-readRDS("R_data/scRNA_CELL_IDENTIFY_REV_240115.rds")

markergene_list<-c("KRT14", "KRT17", "KRT6A", "KRT5", "KRT19", "KRT8","KRT16", "KRT18", "KRT6B", "KRT15", "KRT6C", "KRTCAP3","EPCAM", "SFN",#malignant
                   "FAP", "PDPN", "COL1A2", "DCN", "COL3A1", "COL6A1",#fibroblasts
                   "ACTA1", "ACTN2", "MYL2", "MYH2",#myocytes
                   "CD2", "CD3D", "CD3E", "CD3G", #T cells
                   "SLAMF7", "CD79A", "BLNK", "FCRL5", #B/Plasma cells
                   "CD14", "CD163", "CD68", "FCGR2A", "CSF1R", #macrophages
                   "CD40", "CD80", "CD83", "CCR7", #dendritic cells
                   "CMA1", "MS4A2", "TPSAB1", "TPSB2", #mast cells
                   "PECAM1", "VWF", "ENG" #endothelial cells
)
VlnPlot(scRNA,features = markergene_list,pt.size = 0,stack = T)+NoLegend()
DimPlot(scRNA,label = T)

epma<-subset(scRNA,cell_rev%in%c("Malignant_cells","Epithelial_cells"))
VlnPlot(epma,features =c("IMPDH1","IMPDH2") ,pt.size = 0)+NoLegend()
DotPlot(epma,features = c("IMPDH1","IMPDH2"))
#1.1重新进行归一化(用SCT做！！！)
DefaultAssay(epma) <- "RNA"
#epma<-NormalizeData(epma,normalization.method = "LogNormalize",scale.factor = 10000)
###寻找高变基因
#寻找高变基因的方法：vst、mean.var.plot 和 dispersion（找到的高变基因各不相同）
#epma<-FindVariableFeatures(epma,selection.method = "vst",nfeatures = 2000)

###扩展数据（ScalingData）
#用于PCA降维，一般只需要对高变基因进行Scale（默认使用高变基因）
#all.genes <- rownames(scRNA)
#epma<-ScaleData(epma,vars.to.regress = "percent.mt")#features	要标准化/居中的特征基因，默认是高可变基因;vars.to.regress	要回归的变量， 例如，nUMI 或percent.mito
#1.2 主成分分析
#epma <- RunPCA(epma, features = VariableFeatures(object = epma))
## 1.3可视化
#VizDimLoadings(epma, dims = 1:9, 
              # reduction = "pca") + 
 # theme(axis.text=element_text(size=5), 
  #      axis.title=element_text(size=8,face="bold"))
#DimPlot(epma, reduction = "pca", group.by = "sample.id")
#DimPlot(epma, reduction = "pca", group.by = "tissue.type")

#DimHeatmap(epma, dims = 1:8, nfeatures = 20, cells = 500, balanced = T)
#1.4 Runumap
epma<- SCTransform(epma)
epma<- RunPCA(epma, npcs=50, verbose=FALSE)
#epma <- FindNeighbors(epma, dims = 1:10)
#epma <- FindClusters(epma, resolution = 0.05)
#epma <- RunUMAP(epma,dims = 1:30,reduction = "pca",uwot.sgd = FALSE,metric = 'correlation',min.dist = 0.3,n.neighbors = 30L)

epma<-epma%>%RunUMAP(dims = 1:30,reduction = "pca",metric = 'correlation',min.dist = 0.3,n.neighbors = 50L)%>%FindNeighbors(dims =1:20)%>%FindClusters(resolution = 0.1)%>%identity()

###结果可视化
DimPlot(epma,raster=FALSE,label = T)
DimPlot(epma,raster=FALSE,label = T,group.by = "cell_rev")
DimPlot(scRNA,raster=FALSE,label = T,group.by = "cell_rev")

epma <- RunTSNE(object = epma, dims = 1:20)
TSNEPlot(object = epma, pt.size = 0.5, group.by = "cell_rev",label = F,cols=c("#B4EEB4","#E9967A"))
saveRDS(epma,"R_data/epma_240116_umap_tsne.rds")
epma<-SetIdent(epma,value = "cell_rev")
#1.5 GSVA分析

meta <- epma@meta.data[,c("cell_rev","tissue.type", "sample.id")]#分组信息，为了后续作图
epma_exp <- as.matrix(epma[["SCT"]]@data)
#epma_exp <- as.matrix(epma[["RNA"]]@data)
save(meta,file = "R_data/GSVA_fig2/epma_meta_forheatmap_240115.Rdata")
save(epma_exp,file = "R_data/GSVA_fig2/epma_exp_forrevision_240115.Rdata")

dim(epma_exp)
#加载SGOC基因集
gene_sets <- as.matrix(t(data.table::fread("resource/NUM_genelist.csv")))
#加载exp

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
epma_gsva<-gsva(allcell_exp,gs,kcdf="Gaussian")
#save(epma_gsva,file = "/share/share/wanglixuan/SRA/GSE181919/Rdata_for submission/epma_gsva.Rdata")

##后台运行代码
cd /home/lilixuan/GSE181919-revision/R_scripts
cat >Rscript_HALLMARK_WLX_gsva.R
setwd("/home/lilixuan/GSE181919-revision/")
library(tidyverse)
library(GSVA)
library(dplyr)
library(tidyr)
load(file = "R_data/GSVA_fig2/epma_exp_forrevision_240112.Rdata")
gene_sets <- as.matrix(t(data.table::fread("resource/HALLMARK_WLX.csv")))

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
epma_gsva_HALLMARK<-gsva(epma_exp,gs,kcdf="Gaussian")

save(epma_gsva_HALLMARK,file = "R_data/GSVA_fig2/epma_gsva_HALLMARK.Rdata")
Ctrl+D
#nohup Rscript Rscript_HALLMARK_WLX_gsva.R > nohup_HALLMARK.out 2>&1 &
  

#addmetadata
#epma<-readRDS("R_data/epma_240112_umap_tsne.rds")
epma<-SetIdent(epma,value = "cell_rev")

load(file = "R_data/GSVA_fig2/epma_gsva_SGOC_pathway.Rdata")
#load(file = "R_data/GSVA_fig2/epma_meta_forheatmap_240112.Rdata")
epma_gsva<-as.data.frame(t(epma_gsva))
identical(rownames(epma_gsva),rownames(meta))
epma_SGOCmeta<-cbind(meta,epma_gsva)
epma_SGOCmeta2<-epma_SGOCmeta[,c(1,4:9)]
epma_SGOCmeta2<-rownames_to_column(epma_SGOCmeta2)
colnames(epma_SGOCmeta2)[8]<-"ASGOC"
colnames(epma_SGOCmeta2)[1]<-"ID"
data2 <- epma_SGOCmeta2 %>% 
  pivot_longer(cols=-c(1,2),
               names_to= "pathways",
               values_to = "GSVA_socre")

library(ggplot2)
library(ggpubr)
my_comparisons <- list(
  c("Epithelial_cells","Malignant_cells" ))
ggplot(data = data2,aes(x=cell_rev,y=GSVA_socre))+ #fill=sample填充图形，color=sample填充边线
  geom_violin(aes(fill=cell_rev))+
  geom_boxplot(width = 0.2)+ # Combine with box plot to add median and quartiles
  scale_fill_manual(values=c('Epithelial_cells' = "skyblue", 'Malignant_cells' = "brown1"))+ #使用scale_fill_manual函数自定义颜色,scale_fill_brewer(palette = "Blues")定义单一颜色
  #scale_fill_manual(breaks = 'cell_type_me',values = c("blue", "red"))+
  theme_bw()+
  facet_grid(.~pathways)+
  theme(strip.text.x = element_text(
    size = 10,  face = "bold"), 
    axis.title.x = element_blank(),#去x轴标签
    axis.title.y=element_text(face = "bold",size = 10,colour = "black"),
    axis.text.y = element_text(face = "bold",size = 10,colour = "black",angle=0),
    axis.text.x =  element_text(face = "bold",size = 8,colour = "black" ),# 这里设置x轴方向的字体类型，
    legend.title=element_blank(),#图例名称为空
    legend.direction = "vertical",
    legend.text = element_text(size = 10, face = "bold",colour = "black"), #图例文字大小
    legend.key.size=unit(0.7,'cm') #图例大小
  )+
  ggtitle("GSE181919 cohort") +
  theme(plot.title = element_text(size = 15,face = "bold", vjust = 0.5, hjust = 0.5))+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")),method = "t.test")
#1.6 添加ASGOC分组信息
SGOC<-epma_SGOCmeta2[,c(1,3:8)]
rownames(SGOC)<-SGOC[,1]
SGOC<-data.frame(SGOC[,-1])

epma<-AddMetaData(object = epma,     #seurat对象
                  metadata = SGOC,    #需要添加的metadata
                  col.name = c("serine","folate","methionine","purine","pyrimidine","SGOC"))
head(epma@meta.data,5)
Allmarker=c("SGOC","serine","folate","methionine","purine","pyrimidine")
VlnPlot(epma, features = Allmarker,pt.size=0,stack = T)+ NoLegend() 

##1.7 114_pathway画热图
library(dplyr)
library(tidyverse)
library(pheatmap)
library(tibble)
library(GSEABase)
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(TCGAbiolinks)
library(limma)
library(edgeR)
library(readr)

load(file = "R_data/GSVA_fig2/epma_gsva_114_pathway.Rdata")
epma_gsva_114pat<-as.data.frame(t(epma_gsva_114pat))
identical(rownames(epma_gsva_114pat),rownames(meta))
epma_114meta<-cbind(meta,epma_gsva_114pat)

gsva_Malignant_NC2018<-epma_114meta %>%
  filter(cell_rev=="Malignant_cells") 
gsva_Malignant_NC2018<-gsva_Malignant_NC2018[,-c(1:3)]
gsva_Malignant_NC2018<-t(gsva_Malignant_NC2018)

gsva_Epithelial_NC2018<-epma_114meta %>%
  filter(cell_rev=="Epithelial_cells") 
gsva_Epithelial_NC2018<-gsva_Epithelial_NC2018[,-c(1:3)]
gsva_Epithelial_NC2018<-t(gsva_Epithelial_NC2018)

library(TCGAbiolinks)
DEA.gs <- TCGAanalyze_DEA(
  mat1 = gsva_Epithelial_NC2018,
  mat2 = gsva_Malignant_NC2018,
  metadata = FALSE,
  pipeline = "limma",
  Cond1type = "Epithelial_cells",
  Cond2type = "Malignant_cells",
  fdr.cut = 0.05,
  logFC.cut = F,
)
DEA.gsva_NC2018_scrna_EP_MA<-DEA.gs %>% 
  rownames_to_column("Pathways")
write.csv(DEA.gsva_NC2018_scrna_EP_MA,file = "R_data/GSVA_fig2/DEA.ssgsea_NC2018_scRNA_ep_ma_ssgsea.csv")


pathway114_t0<-DEA.gsva_NC2018_scrna_EP_MA %>%
  filter(t>0)
epma_114meta$rank<-ifelse(epma_114meta$cell_rev=="Epithelial_cells",1,2)
epma_114meta1<-epma_114meta %>%
  rownames_to_column("ID") #%>%
epma_114meta1<-epma_114meta1[order(epma_114meta1$rank),]
rownames(epma_114meta1)<-epma_114meta1[,1]

annotation_col<-as.data.frame(epma_114meta1$cell_rev)
rownames(annotation_col)<-rownames(epma_114meta1)
colnames(annotation_col)[1]<-"group"
heatdata<-epma_114meta1[,-c(1:4,119)]
heatdata<-as.data.frame(t(heatdata))
heatdata<-heatdata%>%
  rownames_to_column("Pathways")
heatdata<-merge(pathway114_t0,heatdata,by="Pathways")
heatdata<-heatdata[order(heatdata$t,decreasing = T),]
rownames(heatdata)<-heatdata[,1]
heatdata<-heatdata[,-c(1:7)]

identical(rownames(annotation_col),colnames(heatdata))
annotation_col$group<-as.factor(annotation_col$group)
ann_colors = list( group = c(Epithelial_cells="#24708B", Malignant_cells="#CDC673"))
pheatmap(heatdata, #热图的数据
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
         cellwidth = 0.05, cellheight = 6,# 格子比例
         fontsize =6)


#1.8 对比得出T与N的差异基因
DefaultAssay(epma) <- "SCT"
epimalignant.markers <- FindMarkers(epma,only.pos =F,ident.1 ="Malignant_cells" ,ident.2 = "Epithelial_cells",threshold=0.05, min.pct=0.1)
epimalignant.markers_1<-rownames_to_column(epimalignant.markers)
epimalignant.markers2<-subset(epimalignant.markers,epimalignant.markers$p_val_adj<0.05)#&abs(epimalignant.markers$avg_log2FC)>0.10)
epimalignant.markers2_1<-rownames_to_column(epimalignant.markers2)
save(epimalignant.markers2,file = "/home/lilixuan/GSE181919-revision/output/epimalignant.markers_logFC0.2.Rdata")
library(ReactomePA)
library(reactome.db)
library(dplyr)


##1.9 ggplot画评分tsne图

tsne = epma@reductions$tsne@cell.embeddings %>%  #坐标信息
  as.data.frame() %>% 
  cbind(score = epma@meta.data$folate)


gb <- ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, color = score)) +
  geom_point() +
  labs(x = "tSNE_1", y = "tSNE_2", color = "score") + 
  ggtitle("folate")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(breaks = seq(1,7,0.5))+
  scale_y_continuous(breaks = seq(4,8,0.5))+
  theme(
    legend.title =element_blank(), #去掉legend.title 
    legend.key=element_rect(fill='white'), #
    legend.text = element_text(size=20), #设置legend标签的大小
    legend.key.size=unit(1,'cm') )  # 设置legend标签之间的大小


gb1<-gb+ theme(panel.grid.major = element_blank(), #主网格线
               panel.grid.minor = element_blank(), #次网格线
               panel.border = element_blank(), #边框
               axis.title = element_blank(),  #轴标题
               axis.text = element_blank(), # 文本
               axis.ticks = element_blank(),
               panel.background = element_rect(fill = 'white'), #背景色
               plot.background=element_rect(fill="white"))
#gb + scale_color_continuous()
mid <- mean(tsne$score)  ## midpoint

gb1 + scale_color_gradient2(midpoint = mid, low = "#00B2EE",
                            mid = "white", high ="#EE2C2C")


##ep_mad的tsne图
celltsne= epma@reductions$tsne@cell.embeddings %>%  #坐标信息
  as.data.frame() %>% 
  cbind(cell_type = epma@meta.data$cell)
ga <- ggplot(celltsne, aes(x = tSNE_1, y = tSNE_2, color = cell_type)) +
  geom_point() +
  labs(x = "tSNE_1", y = "tSNE_2", color = NULL)
ga
ga1<-ga +
  theme(panel.grid.major = element_blank(), #主网格线
        panel.grid.minor = element_blank(), #次网格线
        panel.border = element_blank(), #边框
        axis.title = element_blank(),  #轴标题
        axis.text = element_blank(), # 文本
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), #背景色
        plot.background=element_rect(fill="white"))
ga2<-ga1+         
  theme(
    legend.title = element_blank(), #去掉legend.title 
    legend.key=element_rect(fill='white'), #
    legend.text = element_text(size=20), #设置legend标签的大小
    legend.key.size=unit(1,'cm') ) +  # 设置legend标签之间的大小
  guides(color = guide_legend(override.aes = list(size=5))) #设置legend中 点的大小 

ga2 + scale_color_manual(values = c("#B4EEB4","#E9967A"))

saveRDS(epma,file = "R_data/scripts1_ending_epma.rds")


rm(list = ls())
library(harmony)
epma<-readRDS(file = "/home/lilixuan/GSE181919-revision/R_data/scripts1_ending_epma.rds")
##2.0单独提取肿瘤细胞，进行分群分析
ca_malig<-subset(epma,tissue.type=="CA")
ca_malig<- SCTransform(ca_malig)
ca_malig<- RunPCA(ca_malig, npcs=50, verbose=FALSE)
ca_malig<- RunHarmony(ca_malig,group.by.vars ="sample.id",assay.use="SCT", lambda=1,theta=0, plot_convergence = TRUE,project.dim = FALSE)
ca_malig<-ca_malig%>%RunUMAP(dims = 1:20,reduction = "harmony",metric = 'correlation',min.dist = 0.3,n.neighbors = 30L)%>%FindNeighbors(reduction = "harmony",dims =1:10)%>%FindClusters(resolution = 0.02)%>%identity()
DimPlot(ca_malig,label = T,)

#2.1 绘制t-SNE
ca_malig <- RunTSNE(object = ca_malig, dims = 1:25)
#saveRDS(ca_malig,"NO_LP_AND_LN/ca_malig_DONE_CLUSTER_AND_TSNE.rds")
TSNEPlot(object = ca_malig, pt.size = 1, label = F)#,cols=c("#B4EEB4","#E9967A")group.by = "SCT_snn_res.0.01",
DimPlot(ca_malig,group.by ="SCT_snn_res.0.02" )



#2.2 创建反卷积的数据格式
identical(colnames(ca_malig@assays$RNA@counts),rownames(ca_malig@meta.data))
sce_malig<-SingleCellExperiment(assays = list(counts = as.matrix(ca_malig@assays$RNA@counts)),
                                colData = ca_malig@meta.data)
save(sce_malig,file = "/home/lilixuan/GSE181919-revision/output/MuSiC-240121/sce_malig_only_CA_forMuSiC.Rdata")

load(file = "Rdata_for_submission/MuSiC/sce_malig_cluster_more.Rdata")

##2.3 跑MuSiC的代码 必须加载这些R包！！！！
library(SummarizedExperiment)
library(Biobase)
library(MuSiC)
library(tidyverse)
load(file ="resource/MuSic_tcga_data.Rdata")

c.genes <- sort(intersect(rownames(bulk.mtx), rownames(sce_malig)))

sce_malig <- sce_malig[c.genes,]
bulk.mtx <- as.matrix(bulk.mtx[c.genes,])

dim(sce_malig)
dim(bulk.mtx)
# dim(case.mtx)

Est.prop.TCGA <- music_prop(bulk.mtx = bulk.mtx,sc.sce = sce_malig,clusters = "SCT_snn_res.0.02",samples="patient.id",verbose = F)

NNLS <- Est.prop.TCGA$Est.prop.allgene###MuSiC的结果有两种，源自于两种方法，NNLS和MuSiC
MuSiC <- Est.prop.TCGA$Est.prop.weighted###这两种结果都可以inner_join到临床数据中进行预后分析。

save(NNLS,file = "output/MuSiC-240121/NNLS_ONLY_CA_SCT_snn_res.0.02.Rdata")
save(MuSiC,file = "output/MuSiC-240121/MuSiC_ONLY_CA_SCT_snn_res.0.02.Rdata")
save(Est.prop.TCGA,file = "output/MuSiC-240121/MuSiC_SCT_snn_res.0.02_results.Rdata")
#反卷积发现SCT_snn_res.0.02投射到TCGA有3个亚组且0 2有生存意义

#2.4 将SCT_snn_res.0.02设置成默认的ident
ca_malig<-SetIdent(ca_malig,value = "SCT_snn_res.0.02")
#2.5 画dimplot
DimPlot(ca_malig,raster=FALSE,label = F)
DimPlot(ca_malig,raster=FALSE,label = F,group.by = "patient.id")
#2.6 gglpot2画umap图
umap = ca_malig@reductions$umap@cell.embeddings %>%  #坐标信息
  as.data.frame() %>% 
  cbind(score = ca_malig@meta.data$SGOC)


gb <- ggplot(umap, aes(x = UMAP_1, y = UMAP_2, color = score)) +
  geom_point() +
  labs(x = "UMAP_1", y = "UMAP_2", color = "score") +         
  theme(
    legend.title =element_blank(), #去掉legend.title 
    legend.key=element_rect(fill='white'), #
    legend.text = element_text(size=20), #设置legend标签的大小
    legend.key.size=unit(1,'cm') )  # 设置legend标签之间的大小

gb1<-gb+ theme(panel.grid.major = element_blank(), #主网格线
               panel.grid.minor = element_blank(), #次网格线
               panel.border = element_blank(), #边框
               axis.title = element_blank(),  #轴标题
               axis.text = element_blank(), # 文本
               axis.ticks = element_blank(),
               panel.background = element_rect(fill = 'white'), #背景色
               plot.background=element_rect(fill="white"))
#gb + scale_color_continuous()
mid <- mean(umap$score)  ## midpoint

#gb + scale_color_gradient2(midpoint = mid)


gb1 + scale_color_gradient2(midpoint = mid, low = "#00B2EE",
                            mid = "white", high ="#EE2C2C")

#2.7 GSVA分析emt cell.cycle等
library(GSVA)
meta <- ca_malig@meta.data[,c("cell_rev","tissue.type", "sample.id","SCT_snn_res.0.02")]#分组信息，为了后续作图
maexp <- as.matrix(ca_malig[["SCT"]]@data)
gene_sets <- as.matrix(t(data.table::fread("resource/HALLMARK_WLX.csv")))
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
epma_gsva_HALLMARK<-gsva(maexp,gs,kcdf="Gaussian")

load(file = "output/riskscore/ssgsea_cellcycle_scRNA_CLUSTER_ssgsea.Rdata")
gsva_hallmark<-as.data.frame(t(gsva_hallmark))
identical(rownames(gsva_hallmark),rownames(meta))


ca_malig<-AddMetaData(object = ca_malig,     #seurat对象
                  metadata = gsva_hallmark,    #需要添加的metadata
                  col.name = c("cell_cycle","EMT","Glycolysis","MITOTIC_SPINDLE","OXPHOS"))
head(ca_malig@meta.data,5)
marker=c("EMT","cell_cycle","Glycolysis","MITOTIC_SPINDLE","OXPHOS","SGOC","serine","folate","methionine","purine","pyrimidine")
DotPlot(object = ca_malig,features =marker,cols =c("#C7C7C7",  "#000080"),scale = T)+
  scale_size_area(max_size = 10) + coord_flip()
Allmarker=c("pyrimidine","purine","serine","folate","methionine","SGOC","cell_cycle","EMT")
DotPlot(object = ca_malig,features =Allmarker,cols =c("#C7C7C7",  "#000080"),scale = F)+
  scale_size_area(max_size = 8) + coord_flip()


#2.8 画相关性的图
COR<-ca_malig@meta.data[,c("cell_cycle","SGOC","serine","folate","methionine","purine","pyrimidine")]


scores_cor<-cor(COR,method = 'spearman')
library(corrplot)
col<- colorRampPalette(c("#1C86EE", "white", "#CD5555"))(20)
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
p.mat <- cor.mtest(COR)
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

#2.9 画比例图 提取病人和分组信息

table(ca_malig$patient.id)#查看各组细胞数
prop.table(table(ca_malig$patient.id))
table(ca_malig$patient.id, ca_malig$SCT_snn_res.0.02)#各组不同细胞群细胞数
Cellratio <- prop.table(table(ca_malig$SCT_snn_res.0.02, ca_malig$patient.id), margin = 2)#计算各组样本不同细胞群比例

Cellratio <- as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1))
library(ggplot2)
library(ggsci)
library("scales")
#查看选中配色方案的帮助文档 （以Nature的配色NPG为例）
#?scale_color_jco()+scale_fill_jco()

ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5)+ #,colour = '#222222'
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  #coord_flip()+
  theme(panel.border = element_rect(fill=NA, linetype="solid"))+#color="black",
  scale_color_npg(palette = c("nrc"), alpha = 0.6)+
  scale_fill_npg(palette = c("nrc"), alpha = 0.6)

#2.10 对比各亚组差异基因，并做富集分析
cluster0.02_markergenes<-FindAllMarkers(object =ca_malig,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cluster0.02_markergenes_1<-cluster0.02_markergenes %>%
  filter(p_val_adj<0.05)
save(cluster0.02_markergenes_1,file = "NO_LP_AND_LN/cluster0.02_markergenes_1_only_CA.Rdata")

##2.11 对0.02各组进行GO和KEGG通路分析
library(gplots)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
ids=bitr(cluster0.02_markergenes_1$gene,'SYMBOL','ENTREZID','org.Hs.eg.db')
sce.markers=merge(cluster0.02_markergenes_1,ids,by.x='gene',by.y='SYMBOL')
#接下来进行kegg注释
## 函数split()可以按照分组因子，把向量，矩阵和数据框进行适当的分组。
## 它的返回值是一个列表，代表分组变量每个水平的观测。
gcSample=split(sce.markers$ENTREZID, sce.markers$cluster) 
## KEGG
xx <- compareCluster(gcSample,
                     fun = "enrichKEGG",
                     organism = "hsa", pvalueCutoff = 0.05
)
xx_cluster<-as.data.frame(xx)
write.csv(xx_cluster,"output/Fig3_plot/xx_cluster_KEGG_c0-4.csv")
p <- dotplot(xx)
p + theme(axis.text.x = element_text(
  angle = 45,
  vjust = 0.5, hjust = 0.5
))
library(ggplot2)
cluster_0.02<-read.csv("output/Fig3_plot/cluster_KEGG_c0-4_forplot.csv")

cluster_0.02<-cluster_0.02[order(cluster_0.02$Cluster),]
cluster_0.02_1<-cluster_0.02[c(12,13,10,11,18,16,17,19,9,20:23,4:6,24:26,2,3,15,14,1,8,7),]
cluster_0.02$order<-c(1:28)
rownames(cluster_0.02_1)<-c(1:26)
cluster_0.02$order=factor(rev(as.integer(rownames(cluster_0.02))),labels =  rev(cluster_0.02$Description))

ggplot(cluster_0.02, aes(Cluster, Description)) +
  geom_point(aes(y=order,color=p.adjust, size=Count))+
  theme_bw()+
  theme(axis.text.x=element_text(angle=0,hjust = 1,vjust=0.5))+
  scale_color_gradient(low = '#CD2626', high = '#C1CDCD')+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=10))+
  scale_size_area(max_size = 5)

#GO
xx_GO <- compareCluster(gcSample,
                        fun = "enrichGO",
                        OrgDb = "org.Hs.eg.db",
                        ont = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.01,
                        qvalueCutoff = 0.05
)
GOBP_CA_malig<-as.data.frame(xx_GO)
write.csv(GOBP_CA_malig,"投稿补充/only_CA/GOBP_CA_malig.csv")
p <- dotplot(xx_GO)
p + theme(axis.text.x = element_text(
  angle = 45,
  vjust = 0.5, hjust = 0.5
))

cluster_0.02<-read.csv("投稿补充/only_CA/GOBP_sc_choose.csv")
#cluster_0.02<-cluster_0.02[order(cluster_0.02$Cluster),]
cluster_0.02_1<-cluster_0.02[c(4,9,10,13,14,16,17:21,36,1:3,5:8,22,24:35,12,15,23,37,11),]

rownames(cluster_0.02_1)<-c(1:37)
cluster_0.02_1$order=factor(rev(as.integer(rownames(cluster_0.02_1))),labels =  rev(cluster_0.02_1$Description))

ggplot(cluster_0.02_1, aes(Cluster, Description)) +
  geom_point(aes(y=order,color=p.adjust, size=Count))+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  scale_color_gradient(low = '#CD5555', high = '#1C86EE')+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=10))+
  scale_size_area(max_size = 8)

saveRDS(ca_malig,file = "R_data/Fig3/ca_malig.rds")


#2.12 画比例图 提取病人和分组信息
ca_malig<-readRDS("/home/lilixuan/GSE181919-revision/R_data/Fig3/ca_malig.rds")
DimPlot(ca_malig,raster=FALSE,group.by = "hpv")

table(ca_malig$orig.ident)#查看各组细胞数
prop.table(table(ca_malig$orig.ident))
table(ca_malig$orig.ident, ca_malig$SCT_snn_res.0.02)#各组不同细胞群细胞数
Cellratio <- prop.table(table(ca_malig$SCT_snn_res.0.02, ca_malig$orig.ident), margin = 2)#计算各组样本不同细胞群比例

Cellratio <- as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1))
library(ggplot2)
library(ggsci)
library("scales")
#查看选中配色方案的帮助文档 （以Nature的配色NPG为例）
#?scale_color_jco()+scale_fill_jco()

ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5)+ #,colour = '#222222'
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  #coord_flip()+
  theme(panel.border = element_rect(fill=NA, linetype="solid"))+#color="black",
  scale_color_npg(palette = c("nrc"), alpha = 0.6)+
  scale_fill_npg(palette = c("nrc"), alpha = 0.6)



##3 signature-riskscore
##3.1 加载epma，计算差异基因
epimalignant.markers <- FindMarkers(epma,only.pos =F,ident.1 ="Malignant_cells" ,ident.2 = "Epithelial_cells",threshold=0.25, min.pct=0.1)
epimalignant.markers_1<-rownames_to_column(epimalignant.markers)
epimalignant.markers2<-subset(epimalignant.markers,epimalignant.markers$p_val_adj<0.05&abs(epimalignant.markers$avg_log2FC)>0.3)
epimalignant.markers2_1<-rownames_to_column(epimalignant.markers2)
save(epimalignant.markers2,file = "output/epimalignant.markers_logFC0.3.Rdata")

##3.2 计算riskscore并添加到Seurat文件中
maexp<-as.data.frame(t((ca_malig[["SCT"]]@data)))
maexp[1:5,1:5]
meta <- ca_malig@meta.data[,c("cell_rev", "patient.id")]
identical(rownames(maexp),rownames(meta))
da<-cbind(meta,maexp)
da1<-as.data.frame(t(da))
aaa2<-as.data.frame(data.table::fread("/home/lilixuan/GSE181919-revision/resource/lasso_cox_folum_4GENE.txt"))

rownames(aaa2)<-aaa2[,1]
aaa2<-aaa2[,-1]

tcgaexp <- da1[c(which(rownames(da1)%in%rownames(aaa2))),]
tcgaexp <- tcgaexp[order(rownames(tcgaexp)),]
aaa2<- aaa2[order(rownames(aaa2)),]

#将tcgaexp中数字由“chr”变成“num”
tcgaexp<-as.data.frame(lapply(tcgaexp, as.numeric))
rownames(tcgaexp)<-c(rownames(aaa2))
colnames(tcgaexp)<-c(colnames(da1))
tcgascore <- aaa2$coef*tcgaexp
tcgascore <- t(colSums(tcgascore))
tcgascore <- t(tcgascore)
colnames(tcgascore) <- c('score')

identical(rownames(meta),rownames(tcgascore))
risk<-cbind(meta,tcgascore)
risk<-risk %>% 
  rownames_to_column("ID")
#risk1<-risk[,c(1,2,5)]
risk<-risk[order(risk$score),]
rownames(risk)<-c(1:5284)
data<-tcgascore[,"score"]
quantile(data,c(0.3,0.7))
quantile(data,c(0.5))
# quantile(data,c(0.3,0.7))
#30%       70% 
#  0.2814705 0.6333408   
tcgascore<-as.data.frame(tcgascore)
tcgascore$group<-ifelse(tcgascore$score>0.6333408,"High",
                        ifelse(tcgascore$score<0.2814705,"Low","Median"))
#tcgascore$group2<-ifelse(tcgascore$score>0.4554323,"High","Low")
#将riskscore加入seurat文件中
identical(rownames(tcgascore),rownames(maexp))
ca_malig<-AddMetaData(object = ca_malig,     #seurat对象
                      metadata = tcgascore,    #需要添加的metadata
                      col.name = c("riskscore","group") )#给新添加的metadata命名

ca_malig<-SetIdent(ca_malig,value = "group")
saveRDS(ca_malig,"/home/lilixuan/GSE181919-revision/R_data/Fig3/CA_MALIG_add_riskscore.rds")
H_L_markergene_g2<-FindMarkers(ca_malig,ident.1 = "High",ident.2 = "Low")
H_L_markergene_g2_1<-rownames_to_column(H_L_markergene_g2)
H_L_markergene1<-H_L_markergene_g2_1 %>%
  filter(abs(avg_log2FC)>0.25) %>%
  filter(p_val_adj<0.05)
save(H_L_markergene1,file = "output/riskscore/H_L_markergene_0.25.Rdata")

metadata<-ca_malig@meta.data[,c("SGOC", "group")]  #"riskscore","cell_cycle"
save(metadata,file = "output/Fig4/metadata_cellcycle_sgoc_forcor.Rdata")
library("ggpubr")

ggscatter(metadata, x = "riskscore", y = "cell_cycle", 
          add = "reg.line", conf.int = TRUE,color = "#1874CD",
          add.params = list(color = "#104E8B", fill = "lightgray"),
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "riskscore", ylab = "Cell cycle")

##
metadata<-metadata %>%
  filter(!group=="Median")%>%
  rownames_to_column("ID")
data2 <- metadata %>% 
  pivot_longer(cols=-c(1,3),
               names_to= "pathways",
               values_to = "GSVA_socre")
colnames(data2)[2]<-"sample"
### 肿瘤对比正常作图
my_comparisons <- list(
  c("Low","High" ))
ggplot(data = data2,aes(x=sample,y=GSVA_socre))+ #fill=sample填充图形，color=sample填充边线
  geom_violin(aes(fill=sample))+
  geom_boxplot(width = 0.2)+ # Combine with box plot to add median and quartiles
  #scale_fill_manual(values=c(Normal = "white", OSCC = "white"))+ #使用scale_fill_manual函数自定义颜色,scale_fill_brewer(palette = "Blues")定义单一颜色
  scale_fill_manual(breaks = sample,values = c("brown1", "skyblue"))+
  theme_bw()+
  facet_grid(.~pathways)+
  theme(strip.text.x = element_text(
    size = 15,  face = "bold"), 
    axis.title.x = element_blank(),#去x轴标签
    axis.title.y=element_text(face = "bold",size = 15,colour = "black"),
    axis.text.y = element_text(face = "bold",size = 15,colour = "black",angle=0),
    axis.text.x =  element_text(face = "bold",size = 12,colour = "black" ),# 这里设置x轴方向的字体类型，
    legend.title=element_blank(),#图例名称为空
    legend.direction = "vertical",
    legend.text = element_text(size = 15, face = "bold",colour = "black"), #图例文字大小
    legend.key.size=unit(0.7,'cm') #图例大小
  )+
  ggtitle("sc cohort") +
  theme(plot.title = element_text(size = 15,face = "bold", vjust = 0.5, hjust = 0.5))+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")),method = "t.test")


#3.3 做GO分析
colnames(H_L_markergene1)[1]<-"gene"
diffgene_H <- H_L_markergene1 %>% 
filter(p_val_adj < 0.05) %>% 
  filter(avg_log2FC> 0)

diffgene_L <- H_L_markergene1 %>% 
  filter(p_val_adj < 0.05) %>% 
  filter(avg_log2FC< 0)

gene_L <- diffgene_L$gene
gene_H <- diffgene_H$gene
#基因名称转换，返回的是数据框
gene_H = bitr(gene_H, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
gene_L = bitr(gene_L, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

head(gene_H)
ego_BP_H <- enrichGO(gene = gene_H$ENTREZID,
                     OrgDb= org.Hs.eg.db,
                     ont = "BP",
                     pAdjustMethod = "BH",
                     minGSSize = 1,
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.01,
                     readable = TRUE)
bp_H<-as.data.frame(ego_BP_H)
bp_L<-as.data.frame(ego_BP_L)
write.csv(bp_H,"output/riskscore/bp_H_forxianbiaotu.csv")
write.csv(bp_L,"output/riskscore/bp_L_forxianbiaotu.csv")

library(GOplot)
GOplotIn_BP<-ego_BP_L[c(9,26,28,37,50,83,88,100,150),c(1,2,6,8)]
library(stringr)
GOplotIn_BP$geneID <-str_replace_all(GOplotIn_BP$geneID,'/',',') #把GeneID列中的’/’替换成‘,’
names(GOplotIn_BP)<-c('ID','Term','adj_pval','Genes')#修改列名,后面弦图绘制的时候需要这样的格式
GOplotIn_BP$Category = "BP"#分类信息
GOplotIn_BP<-rownames_to_column(GOplotIn_BP)
GOplotIn_BP<-GOplotIn_BP[,-1]
genedata<-data.frame(ID=diffgene_L$gene,logFC=-diffgene_L$avg_log2FC)
#colnames(genedata)<-c("ID","logFC")
circ_BP<-GOplot::circle_dat(GOplotIn_BP,genedata) 
#circ_BP <- circle_dat(GOplotIn_BP,genedata)#GOplot导入数据格式整理
GOCircle(circ_BP, 
         nsub = 9, # 也可以指定数字
         rad1 = 2, rad2 = 3, # 内径、外径设置
         #zsc.col = c('firebrick3', 'white', 'royalblue3'),# z-score颜色设置
         #lfc.col = c('firebrick3', 'royalblue3'),# 上调下调颜色设置
         lfc.col = "red",
         label.size = 4,
         label.fontface='bold',# 字体大小格式设置
         table.legend = T 
         # 右侧表格设置,为TRUE则无法设置theme
         # 如果解除表格加theme会释放出网格线和坐标轴
         )

#high-risk group
GOplotIn_BP<-ego_BP_H[c(3,4,6,120,146,181),c(1,2,6,8)]
library(stringr)
GOplotIn_BP$geneID <-str_replace_all(GOplotIn_BP$geneID,'/',',') #把GeneID列中的’/’替换成‘,’
names(GOplotIn_BP)<-c('ID','term','adj_pval','genes')#修改列名,后面弦图绘制的时候需要这样的格式
GOplotIn_BP$category = "BP"#分类信息
GOplotIn_BP<-rownames_to_column(GOplotIn_BP)
GOplotIn_BP<-GOplotIn_BP[,-1]
genedata<-data.frame(ID=diffgene_H$gene,logFC=diffgene_H$avg_log2FC)
#colnames(genedata)<-c("ID","logFC")
circ_BP<-GOplot::circle_dat(GOplotIn_BP,genedata) 
#circ_BP <- circle_dat(GOplotIn_BP,genedata)#GOplot导入数据格式整理
GOCircle(circ_BP, 
         nsub = 6, # 也可以指定数字
         rad1 = 2, rad2 = 3, # 内径、外径设置
         #zsc.col = c('firebrick3', 'white', 'royalblue3'),# z-score颜色设置
         #lfc.col = c('firebrick3', 'royalblue3'),# 上调下调颜色设置
         lfc.col = "red",
         label.size = 4,
         label.fontface='bold',# 字体大小格式设置
         table.legend = T 
         # 右侧表格设置,为TRUE则无法设置theme
         # 如果解除表格加theme会释放出网格线和坐标轴
) #弦表图

#3.4 H_L hallmark做GSVA分析
save(maexp,file = "R_data/riskscore/HL_exp_forhallmark.Rdata")
#ssh做
cd /home/lilixuan/scRNA/GSE164690/R_scripts
cat >Rscript_HvsL_hallmark_gsva.R
setwd("/home/lilixuan/scRNA/GSE164690/")
library(tidyverse)
library(GSVA)
library(dplyr)
library(tidyr)
load(file = "R_data/riskscore/HL_exp_forhallmark.Rdata")
pathway <- read_delim("/home/lilixuan/GSE181919-revision/resource/h.all.v7.5.1.symbols_hallmark.gmt", "\t", 
                      escape_double = FALSE, trim_ws = TRUE,col_names = F,show_col_types = FALSE)
pathway <- as.data.frame(pathway)

rownames(pathway)<-pathway$X1
#colnames(pathway)<-pathway[,1]
pathway <- data.frame(t(pathway))
pathway <- pathway[3:202,]
pathway_list <- vector("list",length(pathway))
pathway_list <- lapply(pathway, function(x) {
  unique(na.omit(x)) })
gsva_pathway <- gsva(as.matrix(t(maexp)), pathway_list,method='ssgsea',
                     kcdf='Gaussian',abs.ranking=TRUE)

gsva_hallmark <-as.data.frame(t(gsva_pathway))
save(gsva_hallmark,file = "output/riskscore/ssgsea_hallmark_scRNA_high_low_ssgsea.Rdata")
#nohup Rscript Rscript_HvsL_hallmark_gsva.R > nohup_HL_hallmark.out 2>&1 &
### TCGAanalyze_DEA 函数用于执行差异表达分析（DEA）来识别差异表达基因（DEG）

tcgascore<-tcgascore %>%
  rownames_to_column("ID")

gsva_hallmark_high<-gsva_hallmark %>% 
  rownames_to_column("ID") %>% 
  inner_join(tcgascore,"ID") %>% 
  select(ID,group,everything()) %>% 
  filter(group=="High") %>% 
  column_to_rownames("ID") %>% 
  select(-c("group","score"))
gsva_hallmark_high<-t(gsva_hallmark_high)

gsva_hallmark_low<-gsva_hallmark %>% 
  rownames_to_column("ID") %>% 
  inner_join(tcgascore,"ID") %>% 
  select(ID,group,everything()) %>% 
  filter(group=="Low") %>% 
  column_to_rownames("ID") %>% 
  select(-c("group","score"))
gsva_hallmark_low<-t(gsva_hallmark_low)
###
DEA.gs <- TCGAanalyze_DEA(
  mat1 = gsva_hallmark_low,
  mat2 = gsva_hallmark_high,
  metadata = FALSE,
  pipeline = "limma",
  Cond1type = "Low",
  Cond2type = "High",
  fdr.cut = 0.05,
  logFC.cut = F,
)

DEA.gsva_hallmark_scRNA_high_low<-DEA.gs %>% 
  rownames_to_column("Pathways")
write.csv(DEA.gsva_hallmark_scRNA_high_low,"output/Fig4/DEA.gsva_hallmark_scRNA_high_low.csv")

hall<-read.csv("output/Fig4/ALL_dataset_hallmark_gsva.csv")

#hall<-hall[order(hall$cluster,decreasing = T),]
hall$order<-c(1:68)
rownames(hall)<-c(1:68)
hall1<-hall
hall1$order=factor(rev(as.integer(rownames(hall1))),labels =  rev(hall1$Pathways))

hall1$cluster <- factor(hall1$cluster, level=c("scRNA", "TCGA", "GSE41613","GSE65858","GSE42743"))
ggplot(hall1, aes(cluster, Pathways)) +
  geom_point(aes(x=cluster,y=order,color=adj.P.Val, size=t))+
  theme_bw()+
  theme(axis.text.x=element_text(angle=0,hjust = 1,vjust=0.5))+
  scale_color_gradient(low = '#CD2626', high = '#C1CDCD')+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=10))+
  scale_size_area(max_size = 6)

##3.5 小提琴图SGOC
My_levels <- c('Low','Median','High')
Idents(ca_malig) <- factor(Idents(ca_malig), levels= My_levels)
VlnPlot(ca_malig,features = "riskscore",pt.size=0)

#3.6 加载全部细胞进行重新命名
###提取肿瘤细胞的metadata
#DimPlot(object = ca_malig,group.by = "patient.id",pt.size = 2, label = T)
meta_patient<-ca_malig@meta.data[,c("riskscore","patient.id")]
table(meta_patient$patient.id)
meta_patient<-meta_patient %>%
  rownames_to_column("ID")
patient_risk<-aggregate(meta_patient$riskscore, by=list(type=meta_patient$patient.id),mean)
pait<-patient_risk$type

scRNA<-readRDS("/home/lilixuan/GSE181919-revision/R_data/scRNA_CELL_IDENTIFY_REV_240115.rds")
patient_allcell<-subset(scRNA, tissue.type=="CA") #patient.id %in% pait &
patient_allcell<-subset(patient_allcell, cell_rev=="Epithelial_cells",invert = TRUE)
patient_allcell$cell<-ifelse(patient_allcell$cell_rev %in% c("Endothelial_cells","Fibroblasts","Myocytes"),"Matrix_cell",
                             ifelse(patient_allcell$cell_rev=="Malignant_cells","Malignant_cells","Immune_cell"))

patient_allcell<-SetIdent(patient_allcell,value = "cell")
patient_allcell<-RunTSNE(patient_allcell)
TSNEPlot(patient_allcell, pt.size = 1,label=F)
table(patient_allcell$patient.id)#查看各组细胞数
prop.table(table(patient_allcell$patient.id))
table(patient_allcell$patient.id, patient_allcell$cell)#各组不同细胞群细胞数
Cellratio <- prop.table(table(patient_allcell$cell, patient_allcell$patient.id), margin = 2)#计算各组样本不同细胞群比例

Cellratio <- as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1))

Cellratio1 <- Cellratio %>%
  pivot_wider(names_from = "Var1",
              values_from = "Freq")
colnames(Cellratio1)[1]<-"type"
pa_data<-merge(patient_risk,Cellratio1,by="type")
my_data<-pa_data[,c(2,3)]
ggscatter(my_data, x = "x", y = "Immune_cell", 
          add = "reg.line", conf.int = TRUE,color = "#1874CD",
          add.params = list(color = "#104E8B", fill = "lightgray"),
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "riskscore", ylab = "Relativelevel of immune cell infiltration")

#T cell与riskscore cor
table(patient_allcell$patient.id)#查看各组细胞数
prop.table(table(patient_allcell$patient.id))
table(patient_allcell$patient.id, patient_allcell$cell_rev)#各组不同细胞群细胞数
Cellratio <- prop.table(table(patient_allcell$cell_rev, patient_allcell$patient.id), margin = 2)#计算各组样本不同细胞群比例

Cellratio <- as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1))

Cellratio1 <- Cellratio %>%
  pivot_wider(names_from = "Var1",
              values_from = "Freq")
colnames(Cellratio1)[1]<-"type"
pa_data<-merge(patient_risk,Cellratio1,by="type")
my_data<-pa_data[,c(2,3)]
ggscatter(my_data, x = "x", y = "T_cells", 
          add = "reg.line", conf.int = TRUE,color = "#1874CD",
          add.params = list(color = "#104E8B", fill = "lightgray"),
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "riskscore", ylab = "Relativelevel of T cell infiltration")

saveRDS(epma,file="R_data/scripts3_ending_epma.rds")


library(Seurat)
library(SeuratData)
library(CellChat)
library(NMF)
library(ggalluvial)
library(patchwork)
library(svglite)
library(tidyverse)
options(stringsAsFactors = FALSE)
setwd("/home/lilixuan/GSE181919-revision/")
scRNA_181919<-readRDS("/home/lilixuan/GSE181919-revision/R_data/scRNA_181919_forcellchart.rds")
##4.0 提取T细胞并重新注释
Tcell_181919<-subset(scRNA_181919,tissue.type=="CA"&cell_rev=="T_cells")
#标准化流程

Tcell_181919<- SCTransform(Tcell_181919)
Tcell_181919<- RunPCA(Tcell_181919, npcs=50, verbose=FALSE)
Tcell_181919<- RunHarmony(Tcell_181919,group.by.vars ="patient.id",assay.use="pca",  plot_convergence = TRUE,project.dim = FALSE) 
Tcell_181919<-Tcell_181919%>%FindNeighbors(reduction = "harmony",dims =1:30)%>%FindClusters(resolution = 0.5)%>%RunUMAP(dims = 1:30,reduction = "harmony",metric = 'correlation',min.dist = 0.3,n.neighbors = 30L)%>%identity()

DimPlot(Tcell_181919,raster=FALSE,label = T)
#找各组markergene
library(tidyverse)
library(dplyr)
tcell.markers <- FindAllMarkers(Tcell_181919, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.3)
tcell.markers_top30<-tcell.markers %>%
  group_by(cluster) %>%
  #top_n(avg_log2FC,n=30)
  slice_max(n = 30, order_by = avg_log2FC)
seq<-rep(c(1:30),times=12)
tcell.markers_top30<-data_frame(seq,tcell.markers_top30$cluster,tcell.markers_top30$gene)
tcell.markers_top30_width<-pivot_wider(tcell.markers_top30,names_from = "tcell.markers_top30$cluster",values_from = "tcell.markers_top30$gene")
write.csv(tcell.markers_top30_width, "/home/lilixuan/GSE181919-revision/output/Immune/tcell.markers_top30_12cluster240216.csv", row.names = F)
#T细胞亚型注释
marker.list<-c("CD8A","GZMA","GZMB","GZMH","GZMK","PRFI","CCL5",#cd8+
               "PDCD1","LAG3","TIGIT","CTLA4","CXCR3","HAVCR2",#CD8+ exhausted
               "CD4","FOXP3","IL2RA",#Treg,CD25
               "CCR7","TCF7","CD154"#CD4+
               )
marker.list2<-c("CRTAM",#activated CD8+
                "IFNG","ISG15",#CD8+ effector CD8
                "KLRG1","HSPA1A",#effector stress-state CD8+
               "MKI67",#transitory exhausted
               "GZMK","CXCR3",#precursor exhausted
               "IL7R",#memory CD8+
               "ENTPD1","HAVCR2" #terminally exhausted    wangzhi TOD2
)
VlnPlot(Tcell_181919,features = marker.list2,pt.size = 0,stack = T)+NoLegend()
Tcell_181919$seurat_clusters<- recode(Tcell_181919$seurat_clusters,
                                '0'='CD4+ T','2'='CD4+ T','7'='CD4+ T','9'='CD4+ T','10'='CD4+ T',
                                '1'='CD8+ effect T','5'='Cycling',
                                '4'='CD8+ exhausted T','11'='CD8+ exhausted T','8'='CD8+ exhausted T',
                                '3'='NK','6'='NK')
Tcell_181919<-SetIdent(Tcell_181919,value = "CELL_240217")
Tcell_181919$CELL_240217<-Tcell_181919$seurat_clusters
DimPlot(Tcell_181919,raster=FALSE,label = T)
saveRDS(Tcell_181919,"/home/lilixuan/GSE181919-revision/R_data/Tcell_181919_240217.rds")

#提取CD8细胞再次进行分群
CD8_T<-subset(Tcell_181919,seurat_clusters%in%c("CD4+ T","NK"),invert = TRUE)
#标准化流程
CD8_T<- SCTransform(CD8_T)
CD8_T<- RunPCA(CD8_T, npcs=50, verbose=FALSE)
CD8_T<- RunHarmony(CD8_T,group.by.vars ="patient.id",assay.use="pca",  plot_convergence = TRUE,project.dim = FALSE) 
CD8_T<-CD8_T%>%FindNeighbors(reduction = "harmony",dims =1:30)%>%FindClusters(resolution = 1)%>%RunUMAP(dims = 1:30,reduction = "harmony",metric = 'correlation',min.dist = 0.3,n.neighbors = 30L)%>%identity()

DimPlot(CD8_T,raster=FALSE,label = T)

cd8tcell.markers <- FindAllMarkers(CD8_T, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.3)
cd8tcell.markers_top30<-cd8tcell.markers %>%
  group_by(cluster) %>%
  #top_n(avg_log2FC,n=30)
  slice_max(n = 30, order_by = avg_log2FC)
seq<-rep(c(1:30),times=14)
cd8tcell.markers_top30<-data_frame(seq,cd8tcell.markers_top30$cluster,cd8tcell.markers_top30$gene)
cd8tcell.markers_top30_width<-pivot_wider(cd8tcell.markers_top30,names_from = "cd8tcell.markers_top30$cluster",values_from = "cd8tcell.markers_top30$gene")
write.csv(cd8tcell.markers_top30_width, "/home/lilixuan/GSE181919-revision/output/Immune/cd8tcell.markers_top30_width_14cluster240216.csv", row.names = F)
CD8_T$seurat_clusters<-CD8_T$SCT_snn_res.1
CD8_T$seurat_clusters<- recode(CD8_T$seurat_clusters,
                               '4'='CD4+ T','11'='CD4+ T',
                               '3'='CD8+ effect T','6'='CD8+ effect T','10'='CD8+ effect T','13'='CD8+ effect T',
                               '2'='Cycling',
                               '0'="Naive",
                               '1'='CD8+ exhausted T','5'='CD8+ exhausted T','7'='CD8+ exhausted T',
                               '8'='CD8+ exhausted T','9'='CD8+ exhausted T','12'='CD8+ exhausted T')
CD8_T<-SetIdent(CD8_T,value = "seurat_clusters")
DimPlot(CD8_T,raster=FALSE,label = T)
Idents(Tcell_181919,cells=colnames(CD8_T))<-Idents(CD8_T)
Tcell_181919$CELL_240217<-Tcell_181919@active.ident
Tcell_181919<-SetIdent(Tcell_181919,value = "CELL_240217")
DimPlot(Tcell_181919,raster=FALSE,label = T)
saveRDS(Tcell_181919,file = "R_data/Tcell_181919_200217_identifytcell.rds")
#4.0 把注释后的肿瘤细胞放回到最开始的大群里
#加载ca_malig
ca_malig<-readRDS("/home/lilixuan/GSE181919-revision/R_data/Fig3/ca_malig.rds")
dim(ca_malig)
scRNA_181919$cell_2040217<-scRNA_181919@active.ident
scRNA_181919<-SetIdent(scRNA_181919,value = "cell_2040217")
Idents(scRNA_181919,cells=colnames(ca_malig))<-Idents(ca_malig)
Idents(scRNA_181919,cells=colnames(Tcell_181919))<-Idents(Tcell_181919)
scRNA_181919$cell_2040217<-scRNA_181919@active.ident
DimPlot(scRNA_181919,label = T)
scRNA_181919$CELL_forimmune<-scRNA_181919@active.ident
saveRDS(scRNA_181919,"/home/lilixuan/GSE181919-revision/R_data/scRNA_181919_forcellchart_240217.rds")
#4.1 提取CA的肿瘤细胞和免疫细胞
patient_allcell<-subset(scRNA_181919, tissue.type=="CA") #patient.id %in% pait &
patient_allcell$cell<-patient_allcell@active.ident
patient_allcell<-subset(patient_allcell, cell%in%c('CD8+ effect T','CD8+ exhausted T','High',"Low","Median"))#反选,invert = TRUE
DimPlot(patient_allcell,raster=FALSE,label = F)
patient_allcell<-RunTSNE(patient_allcell)
#标准化流程

patient_allcell<- SCTransform(patient_allcell)
patient_allcell<- RunPCA(patient_allcell, npcs=50, verbose=FALSE)

patient_allcell<-patient_allcell%>%RunUMAP(dims = 1:30,reduction = "pca",metric = 'correlation',min.dist = 0.3,n.neighbors = 30L)%>%FindNeighbors(dims =1:20)%>%FindClusters(resolution = 0.1)%>%identity()
DimPlot(patient_allcell,raster=FALSE,label = T,group.by = "cell")
patient_allcell<-SetIdent(patient_allcell,value = "cell")

saveRDS(patient_allcell,"/home/lilixuan/GSE181919-revision/output/cellchat/patient_allcell_240216.rds")
#4.2 准备cellchat的数据  从Seurat对象直接创建输入的是log后的数据
#分别提取eff和exh的cd8
patient_allcell<-readRDS("/home/lilixuan/GSE181919-revision/output/cellchat/patient_allcell_240216.rds")
eff_patient<-subset(patient_allcell, cell%in%c('CD8+ effect T','High',"Low","Median"))
exh_patient<-subset(patient_allcell, cell%in%c('CD8+ exhausted T','High',"Low","Median"))


##CYY的代码
###创建Cellchat对象
Idents(eff_patient)<-eff_patient@active.ident
eff_patient$sub_cell<-eff_patient@active.ident
DefaultAssay(eff_patient)<-'SCT'
cellchat <- createCellChat(object = eff_patient@assays$SCT@data, meta =eff_patient@meta.data,group.by = "sub_cell")#创建cellchat对象，加上meta注释与细胞类型
cellchat <- setIdent(cellchat, ident.use = "sub_cell")#设置默认的ident
#4.3 导入配体受体数据库
CellChatDB <- CellChatDB.human
str(CellChatDB) #查看数据库信息
#包含interaction、complex、cofactor和geneInfo这4个dataframe
colnames(CellChatDB$interaction)
CellChatDB$interaction[1:4,1:4]
head(CellChatDB$cofactor)
head(CellChatDB$complex)
head(CellChatDB$geneInfo)
#dev.new() #下一步不出图的时候运行
showDatabaseCategory(CellChatDB)

unique(CellChatDB$interaction$annotation)#查看可以选择的侧面，也就是上图左中的三种
#选择"Secreted Signaling"进行后续分析
CellChatDB.use <- CellChatDB


cellchat@DB <- CellChatDB.use# set the used database in the object

# 4.4 预处理
#对表达数据进行预处理,用于细胞间的通信分析。首先在一个细胞组中识别过表达的配体或受
#体，然后将基因表达数据投射到蛋白-蛋白相互作用(PPI)网络上。如果配体或受体过表达，
#则识别过表达配体和受体之间的相互作用。
## 在矩阵的所有的基因中提取signaling gene,结果保存在data.signaling
cellchat <- subsetData(cellchat)
library(future)
future::plan(multisession, workers = 12)
#相当于Seurat的FindMarkers，找每个细胞群中高表达的配体受体
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)#Identify over-expressed ligand-receptor interactions (pairs) within the used CellChatDB
#上一步运行的结果储存在cellchat@LR$LRsig
cellchat <- projectData(cellchat, PPI.human)
#找到配体受体关系后，projectData将配体受体对的表达值投射到PPI上，
#来对@data.signaling中的表达值进行校正。结果保存在@data.project

#4.5 推断细胞通讯网络
#通过为每个相互作用分配一个概率值并进行置换检验来推断生物意义上的细胞-细胞通信。
#推断配体-受体水平细胞通讯网络结果储存在@net下面，有一个概率值和对应的pval）
#这一步也是CellChat比CellPhoneDB多的一步！！！！
#通过计算与每个信号通路相关的所有配体-受体相互作用的通信概率来推断信号通路水平上的
#通信概率。
##根据表达值推测细胞互作的概率（cellphonedb是用平均表达值代表互作强度）。
cellchat <- computeCommunProb(cellchat, raw.use = FALSE,population.size = TRUE)#如果不想用上一步PPI矫正的结果，raw.use = TRUE即可。
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <-filterCommunication(cellchat, min.cells = 10)
df.net1 <- subsetCommunication(cellchat)
write.csv(df.net1,"/home/lilixuan/GSE181919-revision/output/cellchat/net_lr_onlyeffT.csv")
#推断信号通路水平的细胞通讯网络，结果储存在@netP下面，有一个概率值和对应的pval
cellchat <- computeCommunProbPathway(cellchat)
df.netp1 <- subsetCommunication(cellchat,slot.name = "netP")
write.csv(df.netp1,"/home/lilixuan/GSE181919-revision/output/cellchat/df.netp_effT.csv")

#4.6 细胞互作关系展示
#统计细胞和细胞之间通信的数量（有多少个配体-受体对）和强度（概率）
cellchat <-aggregateNet(cellchat)
#计算每种细胞各有多少个
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count,vertex.weight = groupSize,weight.scale = T, 
                 label.edge= F,title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight,vertex.weight = groupSize,weight.scale = T,
                 label.edge= F, title.name="Interaction weights/strength")

#检查每种细胞发出的信号
mat <- cellchat@net$count
par(mfrow = c(1,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0,nrow = nrow(mat), ncol = ncol(mat),dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2,vertex.weight = groupSize, weight.scale = T,
                   arrow.width = 0.2, arrow.size = 0.1,
                   edge.weight.max = max(mat),
                   title.name =rownames(mat)[i])
}


#互作强度和概率
mat <- cellchat@net$weight
par(mfrow = c(1,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0,nrow = nrow(mat), ncol = ncol(mat),dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2,vertex.weight = groupSize, weight.scale = T,
                   arrow.width = 0.2, arrow.size = 0.1,
                   edge.weight.max = max(mat),
                   title.name =rownames(mat)[i])
}
write.csv(mat,"/home/lilixuan/GSE181919-revision/output/cellchat/mat_1qiangdu.csv")
write.csv(mat2,"/home/lilixuan/GSE181919-revision/output/cellchat/mat_2gailv.csv")

saveRDS(cellchat1,file = "/home/lilixuan/GSE181919-revision/output/cellchat/cellchat_effT.rds")

netVisual_circle(cellchat@net$weight,weight.scale = T,
                 sources.use = c("CD8+ effect T"),
                 targets.use = c("Low","Median","High"),
                 label.edge = F,title.name = "interaction weight/strength",remove.isolate=F)

#单个信号通路或配体-受体介导的细胞互作可视化
cellchat@netP$pathways #查看都有哪些信号通路
pathways.show <- c("CXCL") 

#层次图
levels(cellchat@idents)
vertex.receiver = c(1,2,3)
par(mfrow=c(1,1))
netVisual_aggregate(cellchat1, signaling= pathways.show,vertex.receiver=vertex.receiver,layout ="circle") #"circle", "chord" or "spatial"
netVisual_heatmap(cellchat1, signaling= pathways.show,color.heatmap = "Reds") #layout ="chord""circle", "chord" or "spatial"

options(repr.plot.width=10,repr.plot.height=10)
netVisual_heatmap(cellchat,signaling = "CXCL",color.heatmap = "Reds")

# Chord diagram
group.cellType <- c(rep("T", 941), rep("Median", 2073), rep("High", 1594),rep("Low",1617)) # grouping cell clusters into fibroblast, DC and TC cells
names(group.cellType) <- levels(cellchat@idents)
netVisual_chord_cell(cellchat1, signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network"))
#> Plot the aggregated cell-cell communication network at the signaling pathway level
#> Note: The first link end is drawn out of sector 'Inflam. FIB'.
#>
#气泡图
# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(cellchat1, sources.use = 4, targets.use = c(5:11), remove.isolate = FALSE)
#> Comparing communications on a single object
# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_bubble(cellchat1, sources.use = 4, targets.use = c(5:11), signaling = c("CCL","CXCL"), remove.isolate = FALSE)
#> Comparing communications on a single object


pathways.show <- c("CXCL") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat1, signaling = pathways.show,  vertex.receiver = vertex.receiver)

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat1, signaling = pathways.show, color.heatmap = "Reds")
#> Do heatmap based on a single object
#> 
netVisual_bubble(cellchat1, sources.use = 1, targets.use = c(1:4), remove.isolate = FALSE,angle.x =0)


##提取CD8和肿瘤细胞
patient_allcell<-subset(scRNA_181919, tissue.type=="CA")
scRNA_cellphone<-subset(patient_allcell,cell_2040217%in%c("Median","High","Low","Naive","CD8+ exhausted T","Cycling","CD8+ effect T"))
##准备cellphonedb数据
meta_data <- scRNA_cellphone@meta.data %>% 
  rownames_to_column("ID") %>% 
  select(ID,cell_2040217)
write.table(meta_data, 'output/Immune/test_meta_cellphone.txt', sep='\t', quote=F, row.names=F)
#矩阵文件
scRNA_cellphone@assays$RNA@data[1:4,1:4]
write.table(as.matrix(scRNA_cellphone@assays$RNA@data), 'output/Immune/test_count_cellphonedb.txt', sep='\t', quote=F)


#### 2.批量输出Output
library(Seurat)
library(tidyverse)
library(data.table)
library(dplyr)
library(readr)
scRNA_cellphone$celltype<-scRNA_cellphone@active.ident
out.dir = "output/Step8.CellphoneDB/"
for(sample in unique(scRNA_cellphone$patient.id)){
  sp1 <- subset(scRNA_cellphone,patient.id == sample)
  sp1_counts <- as.data.frame(sp1@assays$RNA@data) %>% 
    rownames_to_column(var = "Gene")
  
  sp1_meta <- data.frame(Cell=rownames(sp1@meta.data), 
                         cell_type=sp1@meta.data$celltype)
  
  fwrite(sp1_counts,file = paste0(out.dir,sample,"_counts.txt"), row.names=F, sep='\t')
  fwrite(sp1_meta, file = paste0(out.dir,sample,"_meta.txt"), row.names=F, sep='\t')
}

##进入目的文件夹
#cd Step8.CellphoneDB/
  
  ##先写个脚本
  cat >gc_cellphoneDB.bash
file_count=./$1_counts.txt
file_anno=./$1_meta.txt
outdir=./$1_Output

if [ ! -d ${outdir} ]; then
mkdir ${outdir}
fi
cellphonedb method statistical_analysis ${file_anno} ${file_count} --counts-data hgnc_symbol --output-path ${outdir} --threshold 0.01 --threads 10 
##如果细胞数太多，可以添加下采样参数，默认只分析1/3的细胞
#--subsampling
#--subsampling-log true #对于没有log转化的数据，还要加这个参数



##
#library(devtools)
#devtools::install_github("ShellyCoder/cellcall")
library(cellcall)
##分亚群进行分析
cd8_all<-subset(scRNA_cellphone,cell_2040217%in%c("Median","High","Low","CD8+ effect T","CD8+ exhausted T"))

cd8_all$cellcall_all<-cd8_all@active.ident
cd8_all$cellcall_all<-recode(cd8_all$cellcall_all,'Median'='Median','High'='High','Low'='Low',
                             'CD8+ effect T'='CD8T','CD8+ exhausted T'="CD8T")
cd8_all<-SetIdent(cd8_all,value ="cellcall_all")
##从Seurat对象中创建数据格式
cd8_all_run <- CreateObject_fromSeurat(Seurat.object = cd8_all,
                                       slot = "counts",
                                       cell_type = "cellcall_all",  #细胞类型
                                       data_source = "UMI",
                                       scale.factor = 10^6,
                                       Org = "Homo sapiens")  #物种信息

##使用TransCommuProfile函数计算细胞间通讯的相关性————时间很长
cd8all_com <- TransCommuProfile(object = cd8_all_run,
                                pValueCor = 0.05,
                                CorValue = 0.1,
                                topTargetCor=1,
                                p.adjust = 0.05,
                                use.type="median",
                                probs = 0.9,
                                method="weighted",
                                IS_core = TRUE,   #IS_core参数设置是否只考虑核心基因
                                Org = 'Homo sapiens')
#CTL_cmuni是一个大的CellInter（几G）
saveRDS(cd8all_com,file = "output/Immune/cellcall_cd8_all_comm_cor.rds")

##圈图可视化
cell_color <- data.frame(color=c("#00F5FF","#BC8F8F","#20B2AA","#FF34B3"), stringsAsFactors = FALSE)  #多少种细胞设置多少种颜色,
#"#FFA500", "#FFD700", "#7FFF00", "#DC143C", "#8A2BE2", 
#"#FF4500", "#6495ED", "#FF69B4", "#DAA520", "#FF0000"
rownames(cell_color) <- c("Median","High","Low","CD8_T")
par(mfrow = c(1,2), xpd=TRUE)
ViewInterCircos(object = cd8all_com, font = 2, cellColor = cell_color, 
                lrColor = c("#F16B6F", "#84B1ED"),
                arr.type = "big.arrow", arr.length = 0.04,
                trackhight1 = 0.05, slot="expr_l_r_log2_scale",
                linkcolor.from.sender = TRUE,
                linkcolor = NULL, gap.degree = 0.5,
                trackhight2 = 0.032, track.margin2 = c(0.1, 0.1), DIY = FALSE)


###细胞间通路分析————亮点——时间比较长
n <- cd8_com@data$expr_l_r_log2_scale   #提取表达数据
write.csv(n,"output/Immune/n_expr_Cd8effect.csv")
pathway.hyper.list <- lapply(colnames(n), function(i){
  print(i)
  tmp <- getHyperPathway(data = n, object = cd8_com, cella_cellb = i, Org="Homo sapiens")
  return(tmp)
})


#绘制气泡图
a1<-as.data.frame(pathway.hyper.list[[14]][["NES"]])
colnames(a1)<-"NES"
a2<-as.data.frame(pathway.hyper.list[[15]][["NES"]])
colnames(a2)<-"NES"
a<-rbind(a1,a2)
mid <- mean(a2$NES)
names(pathway.hyper.list) <- colnames(n)
myPub.df <- getForBubble(pathway.hyper.list[13:15], cella_cellb=colnames(n)[13:15])
plotBubble(myPub.df)

viewPheatmap(object = cd8_com, slot="expr_l_r_log2_scale", show_rownames = T,
             show_colnames = T,treeheight_row=0, treeheight_col=10,
             cluster_rows = T,cluster_cols = F,fontsize = 4,angle_col = "45",  
             main="score")


#将list变成dataframe
pathway.hyper.list[[7]]                  #提取第一个list对象
unlist(pathway.hyper.list[[7]])          #将第一个list对象展开成普通向量
length(pathway.hyper.list[[7]])          #第一个list对象的长度
names(unlist(pathway.hyper.list[[7]]))
data_frame <- as.data.frame(do.call(rbind, pathway.hyper.list))
