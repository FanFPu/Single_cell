#================= SoupX ==========================
library(Seurat)
library(dplyr)
library(patchwork)
library(DoubletFinder)
library(ggplot2)
library(future)
library(Matrix)
library(Cairo)
options(bitmapType='cairo')
library(SoupX)
library(Seurat)
library(DropletUtils)

##参数简介
#toc是分析矩阵，即有过滤的矩阵
#tod是全矩阵，即没有任何过滤的矩阵
#rho是污染比例系数，可自行设置，如果不设置则会自动计算
#加载数据
# run_soupx <- function(toc,tod,rho=NULL) {
# toc <- Read10X('/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC/Test/SoupX/toc1_matrix',gene.column=1)
# tod <- Read10X('/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC/Test/SoupX/tod1_matrix',gene.column=1)

indir <- '/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC/SoupX/'
Patients <- "HPCP6"
workdir <- paste0(indir,Patients)
setwd(workdir)
getwd()

#读取tod全矩阵
load(paste0(workdir,'/RawMat.RData')) #里面有三个对象，但是只要mat_comb
# SampleID <- 'R231113001' #"HPCP1"
# seurat_comb <- CreateSeuratObject(counts = as.matrix(mat_comb), meta.data = metadata_comb, names.field = FALSE,project = SampleID, min.cells = 1, min.features = 1)
tod <- mat_comb
#读取toc分析矩阵
load('/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC/Sample/HPCP6/seurat_HPCP6_single.RData') #对象为seurat_comb_singlet
toc <- seurat_comb_singlet@assays$RNA@data

#保证基因名一致
tod <- tod[rownames(toc),]

#提取聚类后的meta.data信息
all <- seurat_comb_singlet
matx <- all@meta.data
#
sc = SoupChannel(tod, toc)
sc = setClusters(sc, setNames(matx$seurat_clusters, rownames(matx)))
#寻找污染比例系数 1.自动寻找 2.根据基因表达细胞得情况进行计算污染比例系数
sc = autoEstCont(sc,forceAccept =T,doPlot =F)
rho <- unique(sc$metaData$rho)
# method1:auto sc
pdf(file = "auto_rho.pdf")
print(autoEstCont(sc))
dev.off()
# # method2: Define with genes，此处的基因列表选择Pomc基因，可替换为任何基因
# pomc_gene=list(pomc=c('EPACM'))
# ute2 = estimateNonExpressingCells(sc,pomc_gene)
# sc_pomc2 = calculateContaminationFraction(sc,pomc_gene,ute2,forceAccept=TRUE)
# #此处设置为上述建议得污染系数
# sc_filter = setContaminationFraction(sc,0.16)

###
sc = setContaminationFraction(sc, rho,forceAccept=TRUE)
#保存取出污染后得
save(sc,file = "seurat_SoupX.RData")
#校正矩阵
out = adjustCounts(sc)
# #保存校正后的矩阵
save(out,file= "Matrix.Rdata")

#auto
auto <- CreateSeuratObject(out, min.cells = 1, min.features=1)
auto <- NormalizeData(auto,verbose = FALSE)
auto <- FindVariableFeatures(auto, selection.method = "vst", nfeatures = 2000)
auto <- ScaleData(auto)
auto <- RunPCA(auto, features = VariableFeatures(object = auto))
auto <- RunUMAP(auto, dims = 1:10)
# auto <- RunTSNE(auto, dims = 1:10)
auto@meta.data <- seurat_comb_singlet@meta.data
# auto@reductions <- seurat_comb_filter@reductions
seurat <- auto 
save(seurat,file = paste0(Patients,"_soupX_auto.RData"))



pdf(file ="2_pca_reduce.pdf", width = 8, height = 6)
print(DimPlot(auto, reduction = "pca", group.by = "orig.ident"))
print(DimPlot(auto, reduction = "pca", group.by = "group"))
print(DimHeatmap(auto,dims = 1:10, cells = 500, balanced = TRUE))
dev.off()

pdf(file = "2_ElbowPlot.pdf", width = 18, height = 6)
print(ElbowPlot(auto))
dev.off()

## umap and tsne reduction
pdf(file ="3_umap_tsne.pdf",width = 5, height = 4)
print(DimPlot(auto, reduction = "umap", pt.size = 0.2, group.by = "orig.ident", label = TRUE))
print(DimPlot(auto, reduction = "umap", pt.size = 0.2, group.by = "group", label = TRUE))
# print(DimPlot(seurat_comb_filter, reduction = "tsne", pt.size = 0.2, group.by = "orig.ident", label = TRUE))
# print(DimPlot(seurat_comb_filter, reduction = "tsne", pt.size = 0.2, group.by = "group", label = TRUE))
dev.off()

# ## Clusters
# seurat_comb_singlet <- FindNeighbors(seurat_comb_singlet, dims = 1:10)
# seurat_comb_singlet <- FindClusters(seurat_comb_singlet, resolution = 0.4)
# seurat_comb_singlet <- FindClusters(seurat_comb_singlet, resolution = 0.5)
# seurat_comb_singlet <- FindClusters(seurat_comb_singlet, resolution = 0.8)
# seurat_comb_singlet <- FindClusters(seurat_comb_singlet, resolution = 1)
# seurat_comb_singlet <- FindClusters(seurat_comb_singlet, resolution = 1.2)

pdf(file = "FeaturePlot.pdf",width = 19, height =14)
FeaturePlot(auto, features = c("EPCAM", "KRT3", "KRT14", "MUC1","ACTA2", "MYH11","CNN1","TAGLN","TP63","CDH1","IGLC2","IGLC3","IGHG3","IGHG1","IGHG4","MYH11","NR2F2","CRYAB","LMOD1","TPPP3","TAGLN","ACTA2"))
dev.off()




#所有样本均已经去除污染，进行合并
## Merge
workdir <- '/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC/SoupX'
sampList <- list.dirs(workdir, recursive = FALSE, full.names = FALSE)
pids = read.delim('SID-PID.txt', col.names=1)
pid_0 = rownames(pids)[1]

#首先获得一个数据作为第一个seurat_merge
load(paste0(workdir, '/', pid_0, '/', pid_0, '_soupX_auto.RData'))  #里面的对象是seurat_comb_singlet
seurat_merge = seurat
seurat_merge$Patients = pids[pid_0,1]

for (sample_id in rownames(pids)[-1]) {
	print(sample_id)
	load(paste0(workdir, '/', sample_id, '/', sample_id, '_soupX_auto.RData'))
    seurat$Patients = pids[sample_id,1]
	print('begin merge...')
	seurat_merge = merge(seurat_merge, seurat)
	print(c(table(seurat_merge$orig.ident), table(seurat_merge$Patients), dim(seurat_merge)))
}
save(seurat_merge,file = "/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC/SoupX/Seurat_merge/Data/seurat_merge_SoupX_raw.RData")


workdir <- '/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC/SoupX/Seurat_merge'
if(!dir.exists(paste0(workdir, "/Results/1.1.4_seurat_notharmony"))){dir.create(paste0(workdir, "/Results/1.1.4_seurat_notharmony"), recursive = TRUE)}
if(!dir.exists(paste0(workdir, "/Results/sc_transform"))){dir.create(paste0(workdir, "/Results/sc_transform"), recursive = TRUE)}
#QC after Merger，合并后的质控
maxGenes <- as.numeric(8000) #8000
minUMIs <- as.numeric(1000) #1000
minGenes <- as.numeric(500) #500
maxPercent.mt <- as.numeric(10) #10
seurat_merge <- subset(seurat_merge, subset = nCount_RNA>minUMIs & nFeature_RNA > minGenes & nFeature_RNA < maxGenes & percent.mt < maxPercent.mt)
save(seurat_merge, file = paste0(workdir, '/Data/seurat_merge_SoupX_QC.RData'))

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
save(seurat_merge, file = paste0(workdir,'/Data/seurat_merge_SoupX_normalize.RData'))


#magic处理
library(Rmagic)
library(reticulate)
# use_python("/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/software/miniconda/envs/miniconda/bin/python")
use_python("/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/miniconda/envs/gch_base/bin/python")
seurat_merge <- magic(seurat_merge)
seurat_merge@assays$RNA@data <- seurat_merge@assays$MAGIC_RNA@data


###Find Varaiable Gene —— Scale（缩放）—— reduction（PCA，scVI）
#=====================================  scVI ===============================
# scvi <- read.csv('/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/BC/R28BC/X_scVI.csv',row.names=1)
# colnames(scvi) <- paste0('scVI_',c(0:9))
# seurat_merge[['scVI']] <- CreateDimReducObject(embeddings = as.matrix(scvi), key = 'scVI_', assay = 'RNA')
#======================================  PCA  =====================
seurat_merge <- FindVariableFeatures(seurat_merge, selection.method = "vst", nfeatures = 1000)     #找高变基因，让后进行scale，有了sacle才能进行pca，umap、tsne可视化之类的后续
all.genes <- rownames(seurat_merge)
seurat_merge <- ScaleData(seurat_merge,features = all.genes)   ##normalization一般是对文库处理，目的消除一些技术差异；scale一般对基因表达量处理（典型的z-score：表达量减均值再除以标准差），目的是后续分析不受极值影响。
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
save(seurat_merge, file='1.1.4_seurat_merge_UMAP.RData')
print(dim(seurat_merge))


##使用seurat_merge做umap图
load('1.1.4_seurat_merge.RData')
#找高变基因前20
top20genes <- head(VariableFeatures(seurat_merge),20)
#可视化基因表达差异程度（前20）
pdf(file = "Results/1.1.4_seurat_notharmony/1.1.4_vstgenes.pdf", width = 18, height = 6)
pdf(file = "1.1.4_vstgenes.pdf", width = 18, height = 6)
plot1 <- VariableFeaturePlot(seurat_merge)
plot2 <- LabelPoints(plot = plot1, points = top20genes, repel = TRUE)
print(plot1 + plot2)
dev.off()
#=========================
# #查看差异效果是否显著
pdf(file = "Results/1.1.4_seurat_notharmony/1.1.4_pca_reduce.pdf", width = 8, height = 6)
pdf(file = "1.1.4_pca_reduce.pdf", width = 8, height = 6)
# # print(DimPlot(seurat_merge, reduction="scVI", group.by = "orig.ident"))
# # print(DimPlot(seurat_merge, reduction="scVI", group.by = "Patients"))
print(DimPlot(seurat_merge, reduction="pca", raster=FALSE, group.by = "orig.ident"))
print(DimPlot(seurat_merge, reduction="pca",  raster=FALSE,group.by = "Patients"))
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
pdf(file = "1.1.4_umap_tsne_patients.pdf", width = 9, height = 8)
# pdf(file = "1.1.4_umap_tsne_patients1.pdf", width = 9, height = 8)
print(DimPlot(seurat_merge, reduction = "umap", pt.size = 0.1, group.by = "Patients",raster=FALSE,  label = TRUE))
print(DimPlot(seurat_merge, reduction = "umap", pt.size = 0.1, group.by = "group", raster=FALSE, label = TRUE))
#print(DimPlot(seurat_merge, reduction = "tsne", pt.size = 0.1, group.by = "Patients", label = TRUE))
#print(DimPlot(seurat_merge, reduction = "tsne", pt.size = 0.1, group.by = "group", label = TRUE))
dev.off()
#
pdf(file = "Results/1.1.4_seurat_notharmony/1.1.4_RNA_featureplots.pdf", width = 12, height = 5)
pdf(file = "1.1.4_RNA_featureplots.pdf", width = 12, height = 5)
# pdf(file = "1.1.4_RNA_featureplots.pdf", width = 12, height = 5)
FeaturePlot(seurat_merge, features = c("nFeature_RNA", "nCount_RNA"), raster=FALSE)
dev.off()


##resolution(分辨率越低，聚类的数量越少)
# load('1.1.4_seurat_merge.RData')
seurat_merge <- FindNeighbors(seurat_merge, dims = 1:10,reduction = "pca")
# seurat_merge <- FindNeighbors(seurat_merge, dims = 1:10,reduction = "scVI")
seurat_merge <- FindClusters(seurat_merge, resolution = seq(0.4,1.5,0.1),algorithm = 1)


#cluster_umap_Plot，聚类后的umap图

save(seurat_merge, file = "1.1.4_seurat_merge_notharmony.RData")

pdf(file = "umap_cluster.pdf",width = 7, height = 6)
print(DimPlot(seurat_merge, reduction = "umap", pt.size = 0.2, group.by = "SCT_snn_res.0.4", raster=FALSE,label = TRUE))
print(DimPlot(seurat_merge, reduction = "umap", pt.size = 0.2, group.by = "SCT_snn_res.0.5",raster=FALSE, label = TRUE))
print(DimPlot(seurat_merge, reduction = "umap", pt.size = 0.2, group.by = "SCT_snn_res.0.6",raster=FALSE, label = TRUE))
print(DimPlot(seurat_merge, reduction = "umap", pt.size = 0.2, group.by = "SCT_snn_res.0.7",raster=FALSE, label = TRUE))
print(DimPlot(seurat_merge, reduction = "umap", pt.size = 0.2, group.by = "SCT_snn_res.0.8",raster=FALSE, label = TRUE))
print(DimPlot(seurat_merge, reduction = "umap", pt.size = 0.2, group.by = "SCT_snn_res.0.9", raster=FALSE,label = TRUE))
print(DimPlot(seurat_merge, reduction = "umap", pt.size = 0.2, group.by = "SCT_snn_res.1", raster=FALSE,  label = TRUE))
print(DimPlot(seurat_merge, reduction = "umap", pt.size = 0.2, group.by = "SCT_snn_res.1.1",raster=FALSE, label = TRUE))
print(DimPlot(seurat_merge, reduction = "umap", pt.size = 0.2, group.by = "SCT_snn_res.1.2",raster=FALSE, label = TRUE))
print(DimPlot(seurat_merge, reduction = "umap", pt.size = 0.2, group.by = "SCT_snn_res.1.3",raster=FALSE, label = TRUE))
print(DimPlot(seurat_merge, reduction = "umap", pt.size = 0.2, group.by = "SCT_snn_res.1.4",raster=FALSE, label = TRUE))
# print(DimPlot(seurat_merge1, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.1.5", label = TRUE))
dev.off()


#之后寻找每个cluster的差异基因，之后根据差异基因进行确定亚型，还有featurePlot进行画图
#group_by就是按照变量进行分组，ungroup就是取消分组信息，回到默认
setwd("/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC/SoupX/Seurat_merge/Results/1.1.4_seurat_notharmony/Results")
seurat_merge <- FindClusters(seurat_merge, resolution = 0.8,algorithm = 1)
seurat_merge_marker <- FindAllMarkers(seurat_merge, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
seurat_merge_marker %>%
            group_by(cluster) %>%
            dplyr::filter(avg_log2FC > 0.6) %>%
            arrange(desc(avg_log2FC)) %>%
            slice_head(n = 10) %>%
            ungroup() -> top20
write.table(seurat_merge_marker, file = "seurat_merge_marker.txt", sep = "\t", row.names = TRUE)  #file后面写入保存路径
save(seurat_merge_marker,file= 'seurat_merge_marker.RData')
#根据Marker gene 画热图，看那几个cluster是不是一类，根据表达情况
gene <- c(top20$gene,'EPCAM','MUC1') 
seurat_merge <- ScaleData(seurat_merge,features = gene) 

pdf(file = "cluster_Marker_heatmap.pdf",width = 14, height = 18)
print(DoHeatmap(seurat_merge, features = gene))
dev.off()


#细胞注释
new.cluster.ids <- c("Nervous","Epithelial","Mesenchymal","Epithelial","Epithelial","Mesenchymal","Mesenchymal","Epithelial","Immune","Mesenchymal","Mesenchymal","Immune","Epithelial","Epithelial","Epithelial","Epithelial","Immune","Immune","Mesenchymal","Mesenchymal","Epithelial","Epithelial","Epithelial","Melanocytes","Epithelial")


#seurat_merge <- FindNeighbors(seurat_merge, dims = 1:30)
seurat_merge <- FindClusters(seurat_merge, resolution =0.9,algorithm = 1)
names(new.cluster.ids) <- levels(seurat_merge)
seurat_merge <- RenameIdents(seurat_merge, new.cluster.ids)
seurat_merge$cellTypes_new <- Idents(seurat_merge) 
pdf(file = "/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC/SoupX/Seurat_merge/Results/sc_transform/Results/umap_bigType.pdf",width = 7, height = 6)
# print(DimPlot(seurat_merge, reduction = "umap",  pt.size = 0.2, raster = FALSE,label = TRUE))
print(DimPlot(seurat_merge, reduction = "umap", pt.size = 0.2, raster = FALSE,group.by = "cellTypes_new", label = T,label.size = 4,cols=c('#a2cffe', '#87a922', '#ffa62b', '#f8481c', '#cffdbc'))) #cols=c(Bcell='red')'#a6814c', '#a484ac', '#fc86aa', '#952e8f', '#02ccfe','#2000b1', '#009337', '#ad0afd'
dev.off()

### FeaturePlot
PC <- list(immuneCells= c("PTPRC"),
                      macrophage = c("CD163","ADGRE1","CD163", "SLC11A1", "APOC1", "CD86", "CSF1R", "SLCO2B1", "CD68", "F13A1", "CD14", "AIF1", "CD80","FCER1G", "FCGR3A", "TYROBP"),
                       Tcell = c("CD3E","CD2", "CD3D",  "CD3G","CD4", "CD8A", "CD8B", "GZMK", "FOXP3","TRDC","NKG7","CD79A"),
                       Bcell = c("KIF4A","FANCI","CHAF1A", "GTSE1", "ASPM", "SPC25","NCAPG2","POLA2","NCAPD3", "CD19", "MS4A1", "CD79A", "CD79B", "BLNK"), #"CD21", 
                       plasmocyte = c("IGLC2","IGLC3","IGHG3","IGHG1","IGHG4","JCHAIN","IGLC1"),
                       DCs = c("HLA-DRB1","ITGAX","CD83", "HLA-DMA","HLA-DQB1", "HLA-DPA1","HLA-DPB1","AIF1","LST1","FTL","HLA-C"),
                       NKcell = c("GNLY", "FCER1G", "KLRB1", "KLRC1", "AREG", "XCL1"),
                       mastCells = c("RHOH","BTK","FER","GATA2","IL4R","KIT","LCP2","ENPP3","RAC2","LAT2","CD84","LAT","ADGRE2","UNC13D","NDEL1", 
                                     "CPA3","TPSB2","TPSAB1","MS4A2","SLC18A2","IL1RL1","PTGS1","HPGDS"),
                       neutrophils = c("CSF3R","S100A8", "S100A9", "LTF", 'CEBPE', 'CEBPA', 'ETV6', 'FOXP1', 'GFI1', 'SPI1', 'STAT3', 'ELANE', 'GSTM1', 'LCN2',  'MPO', 'PRTN3'),
                       fibroblast = c("PDGFRA", "ACTA2", "SULF1", "TAGLN","TPM1", "GINS1" ,"THY1", "RBP1", "COL1A2", "COL1A1", "C1R", "IGFBP7", "SFRP2", "MGP","C1S", "DCN", "CXCL14", "COL3A1","LUM"), #"CTGF",
                       endothelial = c("VWF","PECAM1", "CD34", "HSPG2", "LDB2",  "PTPRB", "DOCK9", "CDH5", "SELE"), #"CD31","GPR116",
                       muscle = c("ACTA2","MYH11","NR2F2","CRYAB","TPM2","GUCA2B","GUCA2A","LMOD1","TPPP3","TAGLN"),
                       Pericyte = c("KCNJ8","ACTA2","PDGFRB","MCAM","HIGD1B","ANGPT2","VIM","MFGE8","MYO1B","NOTCH3","COX4I2"),
                       Cyclegene = c("MKI67","CDC20", "CENPF","PTTG1","TOP2A","CCNB1","PCNA"),
                       Epithelial = c("EPCAM", "KRT3", "KRT14", "MUC1", "TP63","CDH1"),
                       Epithe_Lum = c("HIF1A","CEBPD","BTG1","KRT8","CD9","AQP3","MUC1","KRT18","SLPI","AGR2","LCN2","CLDN4","ANXA1","CD74"),
                       Epithe_Mam = c("PRLR","CLDN4","CSN3","CSN1S1","KRT19","KRT7","KRT8","KRT18"),
                       Adipocytes = c("ADIPOQ","PNPLA2","CFD","CIDEC","APOC1","LRP1","PLIN1","PCK1","CDO1", "STAT6", "TFE3"),#"CIDEA"
                       tuftCells = c("AZGP1", "PLCG2", "HPGDS", "AVIL", "PAEP", "SH2D6", "BMX", "LRMP") #'LYZ2'
                       )
setwd("/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC/SoupX/Seurat_merge/Results/sc_transform/FeaturesPlot")
n <- 1

marker <- PC
seurat <- seurat_merge
while( as.numeric(n) <= length(names(marker))){
    name <- paste0(names(marker)[n],".png")
    genes <- unlist(marker[n])
    genes_count <- length(genes)
    png(name, width = genes_count * 100, height = genes_count * 20, units = "px", res = 300)
    combined_plot <- ggplot() + theme_void()
    for (gene in genes) {
        p <- FeaturePlot(object = seurat, features = gene,raster = TRUE)
        # p <- p + geom_text(aes(x = 1, y = 1,label = ""), vjust = -0.5, size = 3, nudge_y = 0.1)
        combined_plot <- combined_plot + p 
    }
    combined_plot <- combined_plot + facet_wrap(~ genes, scales = "free", nrow  = 5) #combined_plot +
    ggsave(name, plot = combined_plot, width = 20, height = 14)
    dev.off()
    n <- n + 1
}


#####亚型注释
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