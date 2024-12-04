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
library(harmony)
# ##运行harmony的时候删掉我的Rlib 好像不太对
# current_lib_paths <- .libPaths()
# new_lib_paths <- c(current_lib_paths[1], current_lib_paths[-1])
# .libPaths(new_lib_paths)


#Work directory
indir = '/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC/Sample'
topdir = '/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC/Results/Futher_Legend/Epithelial'
workdir = topdir
if(!dir.exists(workdir)){dir.create(workdir, recursive = TRUE)}
if(!dir.exists(paste0(workdir, "/Results/1.1.4_seurat_harmony"))){dir.create(paste0(workdir, "/Results/1.1.4_seurat_harmony"), recursive = TRUE)}
setwd(workdir)

#加载数据
load('/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC/Data/seurat_harmony/1.1.4_seurat_merge_harmony.RData')


# seurat_merge <- ReadH5AD("/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/BC/R28BC/adata_scvi_raw.h5ad")
###Find Varaiable Gene —— Scale（缩放）—— reduction（PCA，scVI）
#=====================================  scVI ===============================
# scvi <- read.csv('/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/BC/R28BC/X_scVI.csv',row.names=1)
# colnames(scvi) <- paste0('scVI_',c(0:9))
# seurat_merge[['scVI']] <- CreateDimReducObject(embeddings = as.matrix(scvi), key = 'scVI_', assay = 'RNA')
#===================================================================
seurat_merge <- FindVariableFeatures(seurat_merge, selection.method = "vst", nfeatures = 1000)   #找高变基因，让后进行scale，有了sacle才能进行pca，umap、tsne可视化之类的后续
seurat_merge <- ScaleData(seurat_merge)   ##normalization一般是对文库处理，目的消除一些技术差异；scale一般对基因表达量处理（典型的z-score：表达量减均值再除以标准差），目的是后续分析不受极值影响。
seurat_merge <- RunPCA(seurat_merge)

#ElbowPlot 肘图  (只针对pca)
pdf(file = "Results/1.1.4_seurat_harmony/1.1.4_ElbowPlot.pdf", width = 18, height = 6)
# print(ElbowPlot(seurat_merge,reduction = 'scVI'))
print(ElbowPlot(seurat_merge,reduction = 'pca'))
dev.off()


#瘦身,默认会把降维矩阵（pca）删掉，添加dimreduces参数,使pca保存
# seurat_merge <- DietSeurat(seurat_merge,dimreducs = "scVI")
# seurat_merge <- DietSeurat(seurat_merge,dimreducs = "pca")
#可视化PCA（umap和tsne 默认的reduction是pca），根据肘图进行维度确认
# seurat_merge <- RunUMAP(seurat_merge, dims = 1:10,reduction = 'scVI')
#seurat_merge <- RunTSNE(seurat_merge, dims = 1:10,check_duplicates = FALSE,reduction = 'scVI')
seurat_merge <- RunUMAP(seurat_merge, dims = 1:10,reduction = 'pca')
seurat_merge <- RunTSNE(seurat_merge, dims = 1:10,check_duplicates = FALSE,reduction = 'pca')
save(seurat_merge, file='1.1.4_seurat_merge.RData')
print(dim(seurat_merge))
seurat_merge <- Tumor_Epi

##如果结果不是很理想，可以查看一下去批次前的情况，或者用Python的scVI进行降维，然后标准化后，再进行umap和tsne可视化，找临近距离，然后聚类画图

#去批次
seurat_merge <- seurat_merge %>% RunHarmony("Patients", plot_convergence = FALSE)  #这一行我把reduction= pca给删掉了 需要查看查看一些harmony这个包的参数情况
seurat_merge <- RunHarmony(object = seurat_merge, group.by.vars = "Patients", plot_convergence = FALSE, reduction = "pca")
#然后后面所有的reduction的pca全部换成harmony，除了后面的umap或tsne不换以外。而且下面两行的run~也加上reduction+harmony。
seurat_merge <- RunUMAP(seurat_merge, dims = 1:10,reduction = "harmony")
seurat_merge <- RunTSNE(seurat_merge, dims = 1:10,check_duplicates = FALSE,reduction="harmony")
save(seurat_merge, file='Tumor_Epi_harmony.RData')

pdf(file = "umap_tsne_patients.pdf", width = 9, height = 8)
# pdf(file = "1.1.4_umap_tsne_patients.pdf", width = 9, height = 8)
print(DimPlot(seurat_merge, reduction = "umap", pt.size = 0.1, group.by = "Patients", label = TRUE))
print(DimPlot(seurat_merge, reduction = "umap", pt.size = 0.1, group.by = "cluster_Types", label = TRUE))
dev.off()
##使用seurat_merge做umap图
# load('1.1.4_seurat_merge_harmony.RData')
#找高变基因前20
top20genes <- head(VariableFeatures(seurat_merge),20)
#可视化基因表达差异程度（前20）
pdf(file = "Results/1.1.4_seurat_harmony/1.1.4_vstgenes.pdf", width = 18, height = 6)
plot1 <- VariableFeaturePlot(seurat_merge)
plot2 <- LabelPoints(plot = plot1, points = top20genes, repel = TRUE)
print(plot1 + plot2)
dev.off()
#=========================
# #查看差异效果是否显著
pdf(file = "Results/1.1.4_seurat_harmony/1.1.4_pca_reduce.pdf", width = 8, height = 6)
# print(DimPlot(seurat_merge, reduction="scVI", group.by = "orig.ident"))
# print(DimPlot(seurat_merge, reduction="scVI", group.by = "Patients"))
print(DimPlot(seurat_merge, reduction="harmony", group.by = "orig.ident"))
print(DimPlot(seurat_merge, reduction="harmony", group.by = "Patients"))
print(DimHeatmap(seurat_merge,dims = 1:10, cells = 500, balanced = TRUE)) # #热图，500个细胞，画图尺寸为10个，基因分数一样（balance = TRUE）
dev.off()
#==========================
#(额外)因为查看后分了两部分，因此取出来一些查看是不是分别布置在两边
# seurat_merge4 <- subset(seurat_merge,subset = Patients %in% c("HBCP1","HBCP2","HBCP3A","HBCP12","HBCP13","HBCP14"))

#可视化经过pca降维后二维空间的情况
pdf(file = "Results/1.1.4_seurat_harmony/1.1.4_umap_tsne_patients.pdf", width = 9, height = 8)
# pdf(file = "1.1.4_umap_tsne_patients.pdf", width = 9, height = 8)
print(DimPlot(seurat_merge, reduction = "umap", pt.size = 0.1, group.by = "Patients", label = TRUE))
print(DimPlot(seurat_merge, reduction = "umap", pt.size = 0.1, group.by = "group", label = TRUE))
print(DimPlot(seurat_merge, reduction = "tsne", pt.size = 0.1, group.by = "Patients", label = TRUE))
print(DimPlot(seurat_merge, reduction = "tsne", pt.size = 0.1, group.by = "group", label = TRUE))
dev.off()
#
pdf(file = "Results/1.1.4_seurat_harmony/1.1.4_RNA_featureplots.pdf", width = 12, height = 5)
FeaturePlot(seurat_merge, features = c("nFeature_RNA", "nCount_RNA"))
dev.off()
#pdf(file = "Results/1.1.4_seurat_notharmony/1.1.4_umap_tsne_patients_p4.pdf", width = 9, height = 8)
#print(DimPlot(seurat_merge, reduction = "umap", pt.size = 0.1, group.by = "Patients", label = TRUE, cols=c('P1T'='red', 'P22T'='green', 'P35T'='blue', 'P38T'='purple')))
#dev.off()

##resolution(分辨率越低，聚类的数量越少)
# load('/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC/Data/seurat_harmony/1.1.4_seurat_merge_harmony.RData')
seurat_merge <- FindNeighbors(seurat_merge, dims = 1:10,reduction = "harmony")
seurat_merge <- FindClusters(seurat_merge, resolution = seq(0.4,1.5,0.1),algorithm = 1)


#cluster_umap_Plot，聚类后的umap图
pdf(file = "Results/1.1.4_seurat_harmony/1.1.4_umap_cluster.pdf",width = 7, height = 6)
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

save(seurat_merge, file = "1.1.4_seurat_merge_Epithelial_harmony.RData")

#之后寻找每个cluster的差异基因，之后根据差异基因进行确定亚型，还有featurePlot进行画图
#group_by就是按照变量进行分组，ungroup就是取消分组信息，回到默认
# setwd("/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC")
seurat_merge <- FindClusters(seurat_merge, resolution = 0.5,algorithm = 1)
seurat_merge_marker <- FindAllMarkers(seurat_merge, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
seurat_merge_marker %>%
            group_by(cluster) %>%
            dplyr::filter(avg_log2FC > 1) %>%
            arrange(desc(avg_log2FC)) %>%
            slice_head(n = 20) %>%
            ungroup() -> top20
write.table(seurat_merge_marker, file = "seurat_merge_marker.txt", sep = "\t", row.names = TRUE)  #file后面写入保存路径
save(seurat_merge_marker,file= 'seurat_merge_marker.RData')
write.table(top20, file = "top20.txt", sep = "\t", row.names = TRUE) 
#根据Marker gene 画热图，看那几个cluster是不是一类，根据表达情况
# seurat_merge_marker <- read.table('/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC/Data/seurat_harmony/seurat_merge_marker.txt')
pdf(file = "Results/1.1.4_seurat_harmony/cluster_Marker_heatmap1.pdf",width = 14, height = 18)
print(DoHeatmap(seurat_merge, features = top20$gene))
dev.off()








