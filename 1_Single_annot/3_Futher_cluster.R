options(stringsAsFactors = FALSE)
library(Seurat)
library(dplyr)
library(DoubletFinder)
library(limma)
library(Cairo)
library(car)
library(SingleCellExperiment)
library(ggplot2)
library(data.table)
options(bitmapType = "cairo")
library(MuDataSeurat)
library(harmony)

#工作路径
indir = '/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC'
topdir = '/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC'
workdir = '/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC1'
setwd(workdir)
# if(!dir.exists(paste0(workdir, "/Results/1.1.4_seurat_notharmony"))){dir.create(paste0(workdir, "/Results/1.1.4_seurat_notharmony"), recursive = TRUE)}
if(!dir.exists(paste0(workdir, "/Results/1.1.4_seurat_harmony"))){dir.create(paste0(workdir, "/Results/1.1.4_seurat_harmony"), recursive = TRUE)}
if(!dir.exists(paste0(workdir, "/Results/Futher_Legend/"))){dir.create(paste0(workdir, "/Results/Futher_Legend/"), recursive = TRUE)}

#多线程运行
plan()
plan("multisession", workers = as.numeric(2))
# Enable parallelization
plan()
options(future.globals.maxSize= 429496729600)

#====================================================================================
#Read Data ####正常来说 读进来以后就已经走到了聚类后的结果
seurat_merge <- ReadH5AD("/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/BC/R23BC/seurat1.h5ad")
load('/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/BC/R28BC/1.1.4_seurat_merge.RData')
##Normalization——Variable Gene——Scale
seurat_merge <- NormalizeData(seurat_merge)
seurat_merge <- FindVariableFeatures(seurat_merge, selection.method = "vst", nfeatures = 2000)
seurat_merge <- ScaleData(seurat_merge)
#seurat_merge <- RunPCA(seurat_merge)  #因为数据是根据scVI进行降维的，因此不需要跑pca

#ElbowPlot 肘图
pdf(file = "Results/1.1.4_seurat_notharmony/1.1.4_ElbowPlot.pdf", width = 18, height = 6)
print(ElbowPlot(seurat_merge))
dev.off()

#可视化 注意:reduction
seurat_merge <- RunUMAP(seurat_merge, dims = 1:10,reduction="scVI")
seurat_merge <- RunTSNE(seurat_merge, dims = 1:10,check_duplicates = FALSE,reduction="scVI")
seurat_merge <- FindNeighbors(seurat_merge, dims = 1:10,reduction = "scVI")
seurat_merge <- FindClusters(seurat_merge, resolution = 0.9,algorithm = 1)
#===================================================================================================

load('/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC/Data/seurat_harmony/1.1.4_seurat_merge_harmony.RData')
seurat_merge <- FindClusters(seurat_merge, resolution =1.3)
#用subset进行选择聚类信息（自我判断几大类）
#将聚类标题转变为字符串  注意：resolution
cluster <- unique(as.character(seurat_merge@meta.data$RNA_snn_res....))
seurat_merge1 <- subset(seurat_merge,subset = RNA_snn_res.0.5 %in% c("11","12","13","15","16","22","24","28"))
seurat_merge2 <- subset(seurat_merge,subset = RNA_snn_res.0.5 %in% c("21","6","3","4","9"))
seurat_merge3 <- subset(seurat_merge,subset = RNA_snn_res.0.5 %in% c("1","2","7","5","8","14","19","10","17","25","26","27"))
seurat_merge4 <- subset(seurat_merge,subset = RNA_snn_res.0.5 %in% c("0","18","20","23"))
dim(seurat_merge1)
dim(seurat_merge2)
dim(seurat_merge3) #查看三者加起来是不是总的cell
#标准化
seurat_merge1 <- NormalizeData(seurat_merge1)
seurat_merge2 <- NormalizeData(seurat_merge2)
seurat_merge3 <- NormalizeData(seurat_merge3)
save(seurat_merge1, file = paste0(workdir, '/Data/seurat_merge_Immune.raw.RData'))
save(seurat_merge2, file = paste0(workdir, '/Data/seurat_merge_Mesenchymal.raw.RData'))
save(seurat_merge3, file = paste0(workdir, '/Data/seurat_merge_Epithelial.raw.RData'))
save(seurat_merge4, file= paste0(workdir,"/Data/seurat_merge_Unsure.raw.RData"))
load(paste0(workdir, '/Data/seurat_merge_Immune.raw.RData'))
# load(paste0(workdir, '/Data/seurat_merge_Mesenchymal.raw.RData'))
# load(paste0(workdir, '/Data/seurat_merge_Epithelial.raw.RData'))
#查看QC结果
pdf(file = "Results/Futher_Legend/1.1.4_qc_post.pdf", width = 30, height = 6)
print(VlnPlot(seurat_merge1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "RNA_snn_res.0.5"))
print(VlnPlot(seurat_merge1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "Patients"))

print(VlnPlot(seurat_merge2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "RNA_snn_res.0.5"))
print(VlnPlot(seurat_merge2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "Patients"))

print(VlnPlot(seurat_merge3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "RNA_snn_res.0.5"))
print(VlnPlot(seurat_merge3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "Patients"))
dev.off()


#基于cluster类型进一步重新重新进行（标准化）——找高变——缩放——降维（PCA，scVI）——寻找临近距离——聚类
seurat_merge1 <- FindVariableFeatures(seurat_merge1, selection.method = "vst", nfeatures = 1000)
seurat_merge1 <- ScaleData(seurat_merge1)
seurat_merge1 <- RunPCA(seurat_merge1)

seurat_merge2 <- FindVariableFeatures(seurat_merge2, selection.method = "vst", nfeatures = 1000)
seurat_merge2 <- ScaleData(seurat_merge2)
seurat_merge2 <- RunPCA(seurat_merge2)

seurat_merge3 <- FindVariableFeatures(seurat_merge3, selection.method = "vst", nfeatures = 1000)
seurat_merge3 <- ScaleData(seurat_merge3)
seurat_merge3 <- RunPCA(seurat_merge3,features = VariableFeatures(object = seurat_merge3))

#降维以后查看肘图 后续的可视化维度设置
pdf(file = "Results/Futher_Legend/ElbowPlot.pdf", width = 18, height = 6)
print(ElbowPlot(seurat_merge1))
print(ElbowPlot(seurat_merge2))
print(ElbowPlot(seurat_merge3))
dev.off()

#可视化降维  umap和tsne 默认的reduction是pca，   

seurat_merge3 <- RunUMAP(seurat_merge3, dims = 1:10)
seurat_merge2 <- RunUMAP(seurat_merge2, dims = 1:10)
seurat_merge1 <- RunUMAP(seurat_merge1, dims = 1:10)
# seurat_merge1 <- RunTSNE(seurat_merge1, dims = 1:10,check_duplicates = FALSE)
# seurat_merge2 <- RunTSNE(seurat_merge2, dims = 1:10,check_duplicates = FALSE)
# seurat_merge3 <- RunTSNE(seurat_merge3, dims = 1:10,check_duplicates = FALSE)
save(seurat_merge1, file=paste0(workdir,'/Data/seurat_merge_Immune.RData'))
save(seurat_merge2, file=paste0(workdir,'/Data/seurat_merge_Mesenchymal.RData'))
save(seurat_merge3, file=paste0(workdir,'/Data/seurat_merge_Epithelial.RData'))

#=====================去批次=========================
# # #去批次 然后后面所有的reduction的pca全部换成harmony，除了后面的umap或tsne不换以外。而且下面两行的run~也加上reduction+harmony。
# load('/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC/Data/seurat_merge_Immune.RData')
load('/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC1/Data/seurat_merge_Mesenchymal.RData')
seurat_merge1 <- seurat_merge1 %>% RunHarmony("Patients", plot_convergence = FALSE,reduction = "pca") 
seurat_merge1 <- RunUMAP(seurat_merge1, dims = 1:10,reduction = "harmony")
seurat_merge2 <- seurat_merge2 %>% RunHarmony("Patients", plot_convergence = FALSE,reduction = "pca") 
seurat_merge2 <- RunUMAP(seurat_merge2, dims = 1:10,reduction = "harmony")
seurat_merge3 <- seurat_merge3 %>% RunHarmony("Patients", plot_convergence = FALSE,reduction = "pca") 
seurat_merge3 <- RunUMAP(seurat_merge3, dims = 1:10,reduction = "harmony")
# seurat_merge3 <- RunTSNE(seurat_merge1, dims = 1:10,check_duplicates = FALSE,reduction="harmony")
save(seurat_merge2, file=paste0(workdir,'/Data/seurat_merge_Mesenchymal_harmony.RData'))
save(seurat_merge1, file=paste0(workdir,'/Data/seurat_merge_Immune_harmony.RData'))
save(seurat_merge3, file=paste0(workdir,'/Data/seurat_merge_Epithelial_harmony.RData'))


##使用seurat_merge做umap图
# load('seurat_merge1_notharmony.RData')
# load('seurat_merge2_notharmony.RData')
# load('seurat_merge3_notharmony.RData')

#画高变gene
top20genes_1 <- head(VariableFeatures(seurat_merge1),20)
top20genes_2 <- head(VariableFeatures(seurat_merge2),20)
top20genes_3 <- head(VariableFeatures(seurat_merge3),20)
pdf(file = "Results/Futher_Legend/vstgenes.pdf", width = 18, height = 6)
plot1 <- VariableFeaturePlot(seurat_merge1)
plot2 <- LabelPoints(plot = plot1, points = top20genes_1, repel = TRUE)
print(plot1 + plot2)
plot3 <- VariableFeaturePlot(seurat_merge2)
plot4 <- LabelPoints(plot = plot3, points = top20genes_2, repel = TRUE)
print(plot3 + plot4)
plot5 <- VariableFeaturePlot(seurat_merge3)
plot6 <- LabelPoints(plot = plot5, points = top20genes_3, repel = TRUE)
print(plot5 + plot6)
dev.off()

#Plot：pca
pdf(file = "Results/Futher_Legend/pca_reduce.pdf", width = 8, height = 6)
print(DimPlot(seurat_merge1, reduction="pca", group.by = "orig.ident"))
print(DimPlot(seurat_merge1, reduction="pca", group.by = "Patients"))
print(DimPlot(seurat_merge2, reduction="pca", group.by = "orig.ident"))
print(DimPlot(seurat_merge2, reduction="pca", group.by = "Patients"))
print(DimPlot(seurat_merge3, reduction="pca", group.by = "orig.ident"))
print(DimPlot(seurat_merge3, reduction="pca", group.by = "Patients"))
print(DimHeatmap(seurat_merge1,dims = 1:10, cells = 500, balanced = TRUE))
print(DimHeatmap(seurat_merge2,dims = 1:10, cells = 500, balanced = TRUE))
print(DimHeatmap(seurat_merge3,dims = 1:10, cells = 500, balanced = TRUE))
dev.off()


#Plot：umap_tsne_patients
load('/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC_1/Data/seurat_merge_Immune.RData')
pdf(file = "Results/Futher_Legend/umap_tsne_patients_Immune_harmony.pdf", width = 9, height = 8)
# pdf(file = "Results/Futher_Legend/umap_tsne_patients_Immune.pdf", width = 9, height = 8)
print(DimPlot(seurat_merge1, reduction = "umap", pt.size = 0.1, group.by = "Patients", label = TRUE))
print(DimPlot(seurat_merge1, reduction = "umap", pt.size = 0.1, group.by = "group", label = TRUE))
#print(DimPlot(seurat_merge, reduction = "tsne", pt.size = 0.1, group.by = "Patients", label = TRUE))
#print(DimPlot(seurat_merge, reduction = "tsne", pt.size = 0.1, group.by = "group", label = TRUE))
dev.off()
pdf(file = "Results/Futher_Legend/umap_tsne_patients_Mesenchymal_harmony.pdf", width = 9, height = 8)
print(DimPlot(seurat_merge2, reduction = "umap", pt.size = 0.1, group.by = "Patients", label = TRUE))
print(DimPlot(seurat_merge2, reduction = "umap", pt.size = 0.1, group.by = "group", label = TRUE))
#print(DimPlot(seurat_merge, reduction = "tsne", pt.size = 0.1, group.by = "Patients", label = TRUE))
#print(DimPlot(seurat_merge, reduction = "tsne", pt.size = 0.1, group.by = "group", label = TRUE))
dev.off()
pdf(file = "Results/Futher_Legend/umap_tsne_patients_Epithelial_harmony.pdf", width = 9, height = 8)
print(DimPlot(seurat_merge3, reduction = "umap", pt.size = 0.1, group.by = "Patients", label = TRUE))
print(DimPlot(seurat_merge3, reduction = "umap", pt.size = 0.1, group.by = "group", label = TRUE))
#print(DimPlot(seurat_merge, reduction = "tsne", pt.size = 0.1, group.by = "Patients", label = TRUE))
#print(DimPlot(seurat_merge, reduction = "tsne", pt.size = 0.1, group.by = "group", label = TRUE))
dev.off()


pdf(file = "Results/Futher_Legend/RNA_featureplots.pdf", width = 12, height = 5)
print(FeaturePlot(seurat_merge1, features = c("nFeature_RNA", "nCount_RNA")))
print(FeaturePlot(seurat_merge2, features = c("nFeature_RNA", "nCount_RNA")))
print(FeaturePlot(seurat_merge3, features = c("nFeature_RNA", "nCount_RNA")))
dev.off()
#pdf(file = "Results/1.1.4_seurat_notharmony/1.1.4_umap_tsne_patients_p4.pdf", width = 9, height = 8)
#print(DimPlot(seurat_merge, reduction = "umap", pt.size = 0.1, group.by = "Patients", label = TRUE, cols=c('P1T'='red', 'P22T'='green', 'P35T'='blue', 'P38T'='purple')))
#dev.off()

##resolution
seurat_merge1 <- FindNeighbors(seurat_merge1, dims = 1:10,reduction = "harmony")
seurat_merge1 <- FindNeighbors(seurat_merge1, dims = 1:10,reduction = "pca")
seurat_merge1 <- FindClusters(seurat_merge1, resolution = seq(0.5,1.4,0.1),algorithm = 1)
seurat_merge2 <- FindNeighbors(seurat_merge2, dims = 1:10,reduction = "harmony")
seurat_merge2 <- FindNeighbors(seurat_merge2, dims = 1:10,reduction = "pca")
seurat_merge2 <- FindClusters(seurat_merge2, resolution = seq(0.5,1.4,0.1),algorithm = 1)
seurat_merge3 <- FindNeighbors(seurat_merge3, dims = 1:10,reduction = "harmony")
seurat_merge3 <- FindNeighbors(seurat_merge3, dims = 1:10,reduction = "pca")
seurat_merge3 <- FindClusters(seurat_merge3, resolution = seq(0.5,1.4,0.1),algorithm = 1)

save(seurat_merge1, file=paste0(workdir,'/Data/seurat_merge_Immune_harmony.RData'))
save(seurat_merge1, file=paste0(workdir,'/Data/seurat_merge_Immune.RData'))
save(seurat_merge2, file=paste0(workdir,'/Data/seurat_merge_Mesenchymal.RData'))
save(seurat_merge3, file=paste0(workdir,'/Data/seurat_merge_Epithelial.RData'))
save(seurat_merge, file = "seurat_scVI_cluster.RData")

#上皮1.1，其他两0.7
#查看一下聚类后的质量
pdf(file = "Results/Futher_Legend/cluster_qc_post_Immune.pdf", width = 30, height = 6)
plot1 <- FeatureScatter(seurat_merge1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_merge1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(seurat_merge1, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "Patients")
plot4 <- FeatureScatter(seurat_merge1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "Patients")
print(VlnPlot(seurat_merge1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "Patients"))
print(VlnPlot(seurat_merge1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3),group.by ="orig.ident")
print(plot1 + plot2)
print(plot3 + plot4)
dev.off()

pdf(file = "Results/Futher_Legend/cluster_qc_post_Mesenchymal.pdf", width = 30, height = 6)
plot1 <- FeatureScatter(seurat_merge2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_merge2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(seurat_merge2, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "Patients")
plot4 <- FeatureScatter(seurat_merge2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "Patients")
print(VlnPlot(seurat_merge2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "Patients"))
print(VlnPlot(seurat_merge2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3),group.by ="orig.ident")
print(plot1 + plot2)
print(plot3 + plot4)
dev.off()

pdf(file = "Results/Futher_Legend/cluster_qc_post_Epithelial.pdf", width = 30, height = 6)
plot1 <- FeatureScatter(seurat_merge3, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_merge3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(seurat_merge3, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "Patients")
plot4 <- FeatureScatter(seurat_merge3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "Patients")
print(VlnPlot(seurat_merge3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "Patients"))
print(VlnPlot(seurat_merge3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3),group.by ="orig.ident")
print(plot1 + plot2)
print(plot3 + plot4)
dev.off()

#######画聚类图
load(paste0(workdir,'/Data/seurat_merge_Immune.RData'))
load(paste0(workdir,'/Data/seurat_merge_Mesenchymal.RData'))
load(paste0(workdir,'/Data/seurat_merge_Epithelial.RData'))

pdf(file = "Results/Futher_Legend/umap_cluster_Epithelial_harmony.pdf",width = 7, height = 6)
#print(DimPlot(seurat_merge3, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.4", label = TRUE))
print(DimPlot(seurat_merge3, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.5", label = TRUE))
print(DimPlot(seurat_merge3, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.6", label = TRUE))
print(DimPlot(seurat_merge3, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.7", label = TRUE))
print(DimPlot(seurat_merge3, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.8", label = TRUE))
print(DimPlot(seurat_merge3, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.9", label = TRUE))
print(DimPlot(seurat_merge3, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.1",   label = TRUE))
print(DimPlot(seurat_merge3, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.1.1", label = TRUE))
print(DimPlot(seurat_merge3, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.1.2", label = TRUE))
print(DimPlot(seurat_merge3, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.1.3", label = TRUE))
print(DimPlot(seurat_merge3, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.1.4", label = TRUE))
#print(DimPlot(seurat_merge3, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.1.5", label = TRUE))
dev.off()

pdf(file = "Results/Futher_Legend/umap_cluster_Immune_harmony.pdf",width = 7, height = 6)
#print(DimPlot(seurat_merge3, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.4", label = TRUE))
print(DimPlot(seurat_merge1, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.5", label = TRUE))
print(DimPlot(seurat_merge1, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.6", label = TRUE))
print(DimPlot(seurat_merge1, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.7", label = TRUE))
print(DimPlot(seurat_merge1, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.8", label = TRUE))
print(DimPlot(seurat_merge1, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.9", label = TRUE))
print(DimPlot(seurat_merge1, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.1",   label = TRUE))
print(DimPlot(seurat_merge1, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.1.1", label = TRUE))
print(DimPlot(seurat_merge1, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.1.2", label = TRUE))
print(DimPlot(seurat_merge1, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.1.3", label = TRUE))
print(DimPlot(seurat_merge1, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.1.4", label = TRUE))
# print(DimPlot(seurat_merge1, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.1.5", label = TRUE))
dev.off()

pdf(file = "Results/Futher_Legend/umap_cluster_Mesenchymal_harmony.pdf",width = 7, height = 6)
#print(DimPlot(seurat_merge3, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.4", label = TRUE))
print(DimPlot(seurat_merge2, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.5", label = TRUE))
print(DimPlot(seurat_merge2, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.6", label = TRUE))
print(DimPlot(seurat_merge2, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.7", label = TRUE))
print(DimPlot(seurat_merge2, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.8", label = TRUE))
print(DimPlot(seurat_merge2, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.9", label = TRUE))
print(DimPlot(seurat_merge2, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.1",   label = TRUE))
print(DimPlot(seurat_merge2, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.1.1", label = TRUE))
print(DimPlot(seurat_merge2, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.1.2", label = TRUE))
print(DimPlot(seurat_merge2, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.1.3", label = TRUE))
print(DimPlot(seurat_merge2, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.1.4", label = TRUE))
#print(DimPlot(seurat_merge3, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.1.5", label = TRUE))
dev.off()

# save(seurat_merge, file = "seurat_merge_scVI_cluster.RData")
save(seurat_merge1, file=paste0(workdir,'/Data/seurat_merge_Immune_harmony.RData'))
save(seurat_merge2, file=paste0(workdir,'/Data/seurat_merge_Mesenchymal_harmony.RData'))
save(seurat_merge3, file=paste0(workdir,'/Data/seurat_merge_Epithelial_harmony.RData'))











#Feature Plot
length(names(Immune_cell))
n <- 1
setwd("/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/BC/R28BC/Results/SignatureGenes/seurat_merge_FeaturePlot1.1.1")

while( as.numeric(n) < 31 ){
  a <- unlist(Immune_cell[n])
  b <- FeaturePlot(seurat_merge1, features =a)
  c <- paste0(names(Immune_cell)[n],".png")
  ggsave(b,filename = c,width = 25,height = 25,device = "png")
  n= n+1
}


#之后寻找每个cluster的差异基因，之后根据差异基因进行确定亚型，还有featurePlot进行画图
#group_by就是按照变量进行分组，ungroup就是取消分组信息，回到默认
setwd("/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/BC/R23BC")
# seurat_merge2 <- FindNeighbors(seurat_merge2, dims = 1:10,reduction = "pca")
seurat_merge1 <- FindClusters(seurat_merge1, resolution = 0.9,algorithm = 1)
seurat_merge2 <- FindClusters(seurat_merge2, resolution = 0.8,algorithm = 1)
seurat_merge3 <- FindClusters(seurat_merge3, resolution = 0.7,algorithm = 1)

seurat_merge1_marker <- FindAllMarkers(seurat_merge1,resolution = 0.6 ,only.pos = TRUE)
seurat_merge1_marker %>%
            group_by(cluster) %>%
            dplyr::filter(avg_log2FC > 1) %>%
            slice_head(n = 20) %>%
            ungroup() -> top20
            
pdf(file = "Results/Futher_Legend/cluster_Marker_heatmap.pdf",width = 15, height = 16)
print(DoHeatmap(seurat_merge1, features = top20$gene))
dev.off()
write.table(seurat_merge1_marker, file = "/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/BC/R28BC/Results/1.1.4_seurat_futher_cluser/seurat_merge1_marker.txt", sep = "\t", row.names = FALSE) #file后面写入保存路径


#聚类后的feature PLot看不出来的话 就单独看簇的marker gene
seurat_merge1.23_marker <- FindMarkers(seurat_merge1,ident.1 =2)
head(seurat_merge1.23_marker, n = 10)


#细胞注释
seurat_merge3 <- FindNeighbors(seurat_merge3, dims = 1:10)
seurat_merge3 <- FindClusters(seurat_merge3, resolution =0.5)

# # cellTypes_new <- c("Macrophage","Tcell","Macrophage","Tcell","Plasmocyte","Macrophage",
#                   "Macrophage","Tcell","Macrophage","DCS","Bcell","Macrophage","Macrophage",
#                   "DCS","Mastcell","Bcell","Plasmocyte")
# cellTypes_new <- c("Myocyte","Fibroblast","Fibroblast","Fibroblast","Myocyte","Myocyte","Endothelial","Endothelial","Myocyte","Fibroblast","Undefined","Fibroblast",
                  #  "Fibroblast","Pericyte","Myocyte","Fibroblast","Pericyte","Myocyte","Endothelial")
# cellTypes_new <- c("Epithelial","Epithelial","Epithelial","Epithelial","Epithelial","Epithelial","Epithelial","Epithelial","Epithelial","Epithelial","Epithelial","Epithelial","Epithelial","Epithelial","Epithelial","Epithelial","Epithelial","Epithelial","Epithelial","Epithelial","Epithelial","Epithelial","Epithelial")

cellTypes_new <- c("tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","Normal","tumor")
length(cellTypes_new)
names(cellTypes_new) <- levels(seurat_merge3)   #names是将函数名称提取出来
seurat_merge3 <- RenameIdents(seurat_merge3, cellTypes_new)  #通过 Idents 函数查看或修改。而 RenameIdents 函数则用于重命名聚类标识符
#细胞注释图
pdf(file = "Results/Futher_Legend/umap_cluster3.pdf",width = 7, height = 6)
print(DimPlot(seurat_merge3, reduction = "umap",pt.size = 0.2, label = TRUE))#,cols = c("grey","red")
dev.off()
#新增一列，将细胞类型放入数据框新的一列中
seurat_merge1$cellTypes_new <- Idents(seurat_merge1)  #Ident() 函数是用于获取或设置每个细胞的聚类（cluster）标识符的函数
seurat_merge2$cellTypes_new <- Idents(seurat_merge2)
seurat_merge3$cellTypes_new <- Idents(seurat_merge3)
save(seurat_merge1, file = "seurat_merge_Immune.RData")
save(seurat_merge2, file = "seurat_merge_fibroblast.RData")
save(seurat_merge3, file = "seurat_merge_Epithelial.RData")

load("/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/BC/R28BC/seurat_merge_Epithelial.RData")
#csv文件 统计细胞类型在每个病人的数量
seurat_merge =merge(seurat_merge1,y=seurat_merge2)
seurat_merge = merge(seurat_merge,y=seurat_merge3)
table(seurat_merge$cellTypes_new)
table(seurat_merge$Patients,seurat_merge$cellTypes_new)
write.csv(a,"/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/BC/R28BC/Results/Futher_Legend/Patients_1.csv",row.names = TRUE)
##
a <- table(seurat_merge3@meta.data$RNA_snn_res.1.3,seurat_merge3@meta.data$Patients) #看每个cluster里面的Patients的分布
write.table(a, file=paste0(workdir,"/shangpi_cluster.txt"), sep = "\t", col.names = TRUE, quote = TRUE)
#图片修改，标题字体等
plot <- DimPlot(seurat_merge, reduction = "umap", label = TRUE, label.size = 4.5) + xlab("UMAP 1") + ylab("UMAP 2") +
    theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10)))
ggsave(filename = "../output/images/pbmc3k_umap.jpg", height = 7, width = 12, plot = plot, quality = 50)


save(seurat_merge, file = "../output/seurat_merge.RData")

