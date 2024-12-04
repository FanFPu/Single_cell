#qsub -cwd -l vf=200g,p=40 -P P21H28400N0232 -q st.q,st_supermem.q work.sh
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
library(sctransform)

#Work directory
indir = '/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC/Sample'
topdir = '/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC'
workdir = topdir
setwd(workdir)


#多线程分析
plan()
plan("multisession", workers = as.numeric(2))
# Enable parallelization
plan()
options(future.globals.maxSize= 429496729600)

##Merge
pids = read.delim('SID-PID.txt', col.names=1)
pid_0 = rownames(pids)[1]

#load(paste0(indir, '/', '1.1.4_seurat_merge_raw.RData'))
#首先获得一个数据作为第一个seurat_merge
load(paste0(indir, '/', pid_0, '/seurat_', pid_0, '_single.RData'))  #里面的对象是seurat_comb_singlet
seurat_merge = seurat_comb_singlet
seurat_merge$Patients = pids[pid_0,1]


for (sample_id in rownames(pids)[-1]) {
	print(sample_id)
	load(paste0(indir, '/', sample_id, '/seurat_', sample_id, '_single.RData'))
  seurat_comb_singlet$Patients = pids[sample_id,1]
	print('begin merge...')
	seurat_merge = merge(seurat_merge, seurat_comb_singlet)
	print(c(table(seurat_merge$orig.ident), table(seurat_merge$Patients), dim(seurat_merge)))
}
save(seurat_merge,file = "/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC/Data/rawData/1.1.4_seurat_merge_raw_NO_QC.RData")

#QC after Merger，合并后的质控
maxGenes <- as.numeric(8000) #8000
minUMIs <- as.numeric(1000) #1000
minGenes <- as.numeric(500) #500
maxPercent.mt <- as.numeric(10) #10
seurat_merge <- subset(seurat_merge, subset = nCount_RNA>minUMIs & nFeature_RNA > minGenes & nFeature_RNA < maxGenes & percent.mt < maxPercent.mt)
save(seurat_merge, file = paste0(workdir, '/1.1.4_seurat_merge_raw.RData'))

load('1.1.4_seurat_merge_raw.RData')
if(!dir.exists(paste0(workdir, "/Results/1.1.4_seurat_notharmony"))){dir.create(paste0(workdir, "/Results/1.1.4_seurat_notharmony"), recursive = TRUE)}

#质控图
pdf(file = "Results/1.1.4_seurat_notharmony/1.1.4_qc_post.pdf", width = 30, height = 6)
plot1 <- FeatureScatter(seurat_merge, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_merge, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(seurat_merge, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "Patients")
plot4 <- FeatureScatter(seurat_merge, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "Patients")
print(VlnPlot(seurat_merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "Patients"))
print(VlnPlot(seurat_merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3),group.by ="orig.ident")
print(plot1 + plot2)
print(plot3 + plot4)
dev.off()

##normalize,标准化
# load('1.1.4_seurat_merge_raw.RData')
seurat_merge <- NormalizeData(seurat_merge)

##seurat_merge，去基因
print(dim(seurat_merge))
removeBiasGenes <- function(genes){
  RPgenes <- genes[intersect(grep("^RP", genes), grep("-", genes))]
  RPgenes2 <- genes[grep("^RP[SL]", genes)]
  MTgenes <- genes[grep("^MT-", genes)]
  CTCgenes <- genes[intersect(grep("^CTC", genes), grep("-", genes))]
  MIRgenes <- genes[grep("^MIR", genes)]
  ACgenes <- genes[intersect(grep("^AC[0-9]", genes), grep(".", genes))]
  CTgenes <- genes[intersect(grep("^CT", genes), grep("-", genes))]
  LINCgenes <- genes[grep("^LINC[0-9]", genes)]
  ALgenes <- genes[intersect(grep("^AL", genes), grep(".", genes))]
  HEMOgene <- c("HBB", "HBA1", "HBA2")
  rmgenes <- c(RPgenes, RPgenes2, MTgenes, CTCgenes, MIRgenes, ACgenes, CTgenes, LINCgenes, ALgenes, HEMOgene)
  return(rmgenes)
}
rmgenes = removeBiasGenes(rownames(seurat_merge))
seurat_merge = seurat_merge[setdiff(rownames(seurat_merge), rmgenes),]
save(seurat_merge, file = '1.1.4_seurat_merge_normalize.RData')



#magic处理
library(Rmagic)
library(reticulate)
# use_python("/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/software/miniconda/envs/miniconda/bin/python")
use_python("/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/miniconda/envs/gch_base/bin/python")
seurat_merge <- magic(seurat_merge)
seurat_merge@assays$RNA@data <- seurat_merge@assays$MAGIC_RNA@data



###Find Varaiable Gene —— Scale（缩放）—— reduction（PCA，scVI）
load('1.1.4_seurat_merge_normalize.RData')
#=====================================  scVI ===============================
# scvi <- read.csv('/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/BC/R28BC/X_scVI.csv',row.names=1)
# colnames(scvi) <- paste0('scVI_',c(0:9))
# seurat_merge[['scVI']] <- CreateDimReducObject(embeddings = as.matrix(scvi), key = 'scVI_', assay = 'RNA')
#======================================  PCA  =====================
seurat_merge <- FindVariableFeatures(seurat_merge, selection.method = "vst", nfeatures = 1000)     #找高变基因，让后进行scale，有了sacle才能进行pca，umap、tsne可视化之类的后续
seurat_merge <- ScaleData(seurat_merge)   ##normalization一般是对文库处理，目的消除一些技术差异；scale一般对基因表达量处理（典型的z-score：表达量减均值再除以标准差），目的是后续分析不受极值影响。
seurat_merge <- RunPCA(seurat_merge)
#======================================  sctransform  =====================
# store mitochondrial percentage in object meta data
seurat_merge <- PercentageFeatureSet(seurat_merge, pattern = "^MT-", col.name = "percent.mt")
# run sctransform
seurat_merge <- SCTransform(seurat_merge, vars.to.regress = "percent.mt", verbose = FALSE)
seurat_merge <- RunPCA(seurat_merge)
#======================================  CCA  =====================
seurat_merge1 <- subset(seurat_merge, cells = colnames(seurat_merge)[1:50000])
seurat_merge2 <- subset(seurat_merge, cells = colnames(seurat_merge)[50001:102028 ])
seurat_merge1[["group"]] <- "group1"
seurat_merge2[["group"]] <- "group2"
seurat_merge_cca <- RunCCA(object1 = seurat_merge1, object2 = seurat_merge2)
print(x = pbmc_cca[["cca"]])




#ElbowPlot 肘图  (只针对pca)
pdf(file = "Results/1.1.4_seurat_notharmony/1.1.4_ElbowPlot.pdf", width = 18, height = 6)
# print(ElbowPlot(seurat_merge,reduction = 'scVI'))
pdf(file = "1.1.4_ElbowPlot.pdf", width = 18, height = 6)
print(ElbowPlot(seurat_merge,reduction = 'pca'))
dev.off()


#Seurat对象瘦身,默认会把降维矩阵（pca）删掉，添加dimreduces参数,使pca保存
# seurat_merge <- DietSeurat(seurat_merge,dimreducs = "scVI")
seurat_merge <- DietSeurat(seurat_merge,dimreducs = "pca")
#可视化PCA（umap和tsne 默认的reduction是pca），根据肘图进行维度确认
# seurat_merge <- RunUMAP(seurat_merge, dims = 1:10,reduction = 'scVI')
#seurat_merge <- RunTSNE(seurat_merge, dims = 1:10,check_duplicates = FALSE,reduction = 'scVI')
seurat_merge <- RunUMAP(seurat_merge, dims = 1:10,reduction = 'pca')
#seurat_merge <- RunTSNE(seurat_merge, dims = 1:10,check_duplicates = FALSE,reduction = 'pca')
save(seurat_merge, file='1.1.4_seurat_merge.RData')
print(dim(seurat_merge))


##使用seurat_merge做umap图
load('1.1.4_seurat_merge.RData')
#找高变基因前20
top20genes <- head(VariableFeatures(seurat_merge),20)
#可视化基因表达差异程度（前20）
pdf(file = "Results/1.1.4_seurat_notharmony/1.1.4_vstgenes.pdf", width = 18, height = 6)
plot1 <- VariableFeaturePlot(seurat_merge)
plot2 <- LabelPoints(plot = plot1, points = top20genes, repel = TRUE)
print(plot1 + plot2)
dev.off()
#=========================
# #查看差异效果是否显著
pdf(file = "Results/1.1.4_seurat_notharmony/1.1.4_pca_reduce.pdf", width = 8, height = 6)
# # print(DimPlot(seurat_merge, reduction="scVI", group.by = "orig.ident"))
# # print(DimPlot(seurat_merge, reduction="scVI", group.by = "Patients"))
print(DimPlot(seurat_merge, reduction="pca", group.by = "orig.ident"))
print(DimPlot(seurat_merge, reduction="pca", group.by = "Patients"))
print(DimHeatmap(seurat_merge,dims = 1:10, cells = 500, balanced = TRUE)) # #热图，500个细胞，画图尺寸为10个，基因分数一样（balance = TRUE）
dev.off()
# ==========================
#(额外)因为查看后分了两部分，因此取出来一些查看是不是分别布置在两边
# seurat_merge4 <- subset(seurat_merge,subset = Patients %in% c("HBCP1","HBCP2","HBCP3A","HBCP12","HBCP13","HBCP14"))
#pdf(file = "Results/1.1.4_seurat_notharmony/1.1.4_umap_tsne_patients_p4.pdf", width = 9, height = 8)
#print(DimPlot(seurat_merge, reduction = "umap", pt.size = 0.1, group.by = "Patients", label = TRUE, cols=c('P1T'='red', 'P22T'='green', 'P35T'='blue', 'P38T'='purple')))
#dev.off()

#可视化经过pca降维后二维空间的情况
pdf(file = "Results/1.1.4_seurat_notharmony/1.1.4_umap_tsne_patients.pdf", width = 9, height = 8)
pdf(file = "1.1.4_umap_tsne_patients1.pdf", width = 9, height = 8)
print(DimPlot(seurat_merge, reduction = "umap", pt.size = 0.1, group.by = "Patients",raster=FALSE,  label = TRUE))
print(DimPlot(seurat_merge, reduction = "umap", pt.size = 0.1, group.by = "group", raster=FALSE, label = TRUE))
#print(DimPlot(seurat_merge, reduction = "tsne", pt.size = 0.1, group.by = "Patients", label = TRUE))
#print(DimPlot(seurat_merge, reduction = "tsne", pt.size = 0.1, group.by = "group", label = TRUE))
dev.off()
#
pdf(file = "Results/1.1.4_seurat_notharmony/1.1.4_RNA_featureplots.pdf", width = 12, height = 5)
pdf(file = "1.1.4_RNA_featureplots.pdf", width = 12, height = 5)
FeaturePlot(seurat_merge, features = c("nFeature_RNA", "nCount_RNA"), raster=FALSE)
dev.off()


##resolution(分辨率越低，聚类的数量越少)
load('1.1.4_seurat_merge.RData')
seurat_merge <- FindNeighbors(seurat_merge, dims = 1:10,reduction = "pca")
# seurat_merge <- FindNeighbors(seurat_merge, dims = 1:10,reduction = "scVI")
seurat_merge <- FindClusters(seurat_merge, resolution = seq(0.4,1.5,0.1),algorithm = 1)


#cluster_umap_Plot，聚类后的umap图
setwd('/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC')
pdf(file = "Results/1.1.4_seurat_notharmony/1.1.4_umap_cluster2.pdf",width = 7, height = 6)
print(DimPlot(seurat_merge, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.2",raster=FALSE, label = TRUE))
print(DimPlot(seurat_merge, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.4", label = TRUE))
print(DimPlot(seurat_merge, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.5", label = TRUE))
print(DimPlot(seurat_merge, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.6", label = TRUE))
print(DimPlot(seurat_merge, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.7", label = TRUE))
print(DimPlot(seurat_merge, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.8", label = TRUE))
print(DimPlot(seurat_merge, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.9", label = TRUE))
print(DimPlot(seurat_merge, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.1",   label = TRUE))
print(DimPlot(seurat_merge, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.1.1", label = TRUE))
print(DimPlot(seurat_merge, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.1.2", label = TRUE))
print(DimPlot(seurat_merge, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.1.3", label = TRUE))
print(DimPlot(seurat_merge, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.1.4", label = TRUE))
print(DimPlot(seurat_merge, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.1.5", label = TRUE))
dev.off()


save(seurat_merge, file = "1.1.4_seurat_merge_notharmony.RData")


#之后寻找每个cluster的差异基因，之后根据差异基因进行确定亚型，还有featurePlot进行画图
#group_by就是按照变量进行分组，ungroup就是取消分组信息，回到默认
setwd("/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC")
seurat_merge1 <- FindClusters(seurat_merge1, resolution = 0.4,algorithm = 1)
seurat_merge_marker <- FindAllMarkers(seurat_merge1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
seurat_merge_marker %>%
            group_by(cluster) %>%
            dplyr::filter(avg_log2FC > 1) %>%
            arrange(desc(avg_log2FC)) %>%
            slice_head(n = 20) %>%
            ungroup() -> top20
write.table(seurat_merge_marker, file = "seurat_merge_marker.txt", sep = "\t", row.names = TRUE)  #file后面写入保存路径
save(seurat_merge_marker,file= 'seurat_merge_marker.RData')
#根据Marker gene 画热图，看那几个cluster是不是一类，根据表达情况
seurat_merge_marker <- read.table('/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC/Data/seurat_harmony/seurat_merge_marker.txt')
pdf(file = "Results/1.1.4_seurat_notharmony/cluster_Marker_heatmap.pdf",width = 14, height = 18)
print(DoHeatmap(seurat_merge, features = top20$gene))
dev.off()



##=====相关性

table(seurat_merge$SCT_snn_res.0.7)
av <- AverageExpression(seurat_merge, group.by = 'SCT_snn_res.0.7', assays = 'SCT')
av <- av[[1]]
head(av)
#选出标准差最大的1000个基因
cg <- names(tail(sort(apply(av,1,sd)),1000)) #这个部分是在数据集av的每一行上应用sd函数，即计算每一行的标准差。
#查看这1000个基因在各细胞群中的表达矩阵
View(av[cg,])
#查看细胞群的相关性矩阵
View(cor(av[cg,],method='spearman'))
# heatmap
pdf("cluster_correlation.pdf", width = 10, height = 10)
print(pheatmap::pheatmap(cor(av[cg,],method='spearman')))
dev.off()

#聚类后的feature PLot看不出来的话 就单独看簇的marker gene
# seurat_merge0_marker <- FindMarkers(seurat_merge,ident.1 =0)
# write.table(seurat_merge0_marker, file = "/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC/Data/seurat_harmony/seurat_merge_marker_0.txt", sep = "\t", row.names = TRUE) 
# head(seurat_merge0_marker, n = 20)
# ###改变gene类型
# gene_n.df <- bitr(seurat_merge_marker[,7][1],
#         fromType = "ENSEMBL",
#         toType = "SYMBOL",
#         OrgDb = org.Hs.eg.db
#     )








