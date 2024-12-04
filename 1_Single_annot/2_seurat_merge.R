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

#Work directory
indir = '/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC/Sample'
topdir = '/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC_magic/Futher_Legend/Epithelial'
workdir = topdir
if(!dir.exists(paste0(workdir, "/Data/seurat_notharmony"))){dir.create(paste0(workdir, "/Data/seurat_notharmony"), recursive = TRUE)}
if(!dir.exists(paste0(workdir, "/Results/1.1.4_seurat_notharmony"))){dir.create(paste0(workdir, "/Results/1.1.4_seurat_notharmony"), recursive = TRUE)}
#if(!dir.exists(paste0(workdir, "/Results/1.1.4_seurat_harmony"))){dir.create(paste0(workdir, "/Results/1.1.4_seurat_harmony"), recursive = TRUE)}
outdir  = paste0(topdir, "/Data/seurat_notharmony")
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

# #load(paste0(indir, '/', '1.1.4_seurat_merge_raw.RData'))
# # load('/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/BC/Sample/1.1.4_seurat_merge_raw.RData')
# # seurat_merge1 =seurat_merge
# #
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

#QC after Merger，合并后的质控
maxGenes <- as.numeric(8000) #8000
minUMIs <- as.numeric(1000) #1000
minGenes <- as.numeric(500) #500
maxPercent.mt <- as.numeric(10) #10
seurat_merge <- subset(seurat_merge, subset = nCount_RNA>minUMIs & nFeature_RNA > minGenes & nFeature_RNA < maxGenes & percent.mt < maxPercent.mt)
save(seurat_merge, file = paste0(outdir, '/1.1.4_seurat_merge_raw.RData'))

# load(paste0(outdir, '/1.1.4_seurat_merge_raw.RData'))

load("/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC_magic/Futher_Legend/Mesenchymal/Data/seurat_notharmony/1.1.4_seurat_merge.RData")
seurat_merge <- seurat_merge1
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

#normalize,标准化
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
rmgenes=removeBiasGenes(rownames(seurat_merge))
seurat_merge = seurat_merge[setdiff(rownames(seurat_merge), rmgenes),]
save(seurat_merge, file = paste0(outdir, '/1.1.4_seurat_merge_normalize.RData'))


#magic处理
library(Rmagic)
library(reticulate)
use_python("/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/software/miniconda/envs/miniconda/bin/python")
seurat_merge <- magic(seurat_merge)



# seurat_merge <- ReadH5AD("/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/BC/R28BC/adata_scvi_raw.h5ad")
###Find Varaiable Gene —— Scale（缩放）—— reduction（PCA，scVI）
# load(paste0(outdir, '/1.1.4_seurat_merge_normalize.RData'))
#=====================================  scVI ===============================
# scvi <- read.csv('/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/BC/R28BC/X_scVI.csv',row.names=1)
# colnames(scvi) <- paste0('scVI_',c(0:9))
# seurat_merge[['scVI']] <- CreateDimReducObject(embeddings = as.matrix(scvi), key = 'scVI_', assay = 'RNA')
# #===================================================================
seurat_merge <- FindVariableFeatures(seurat_merge, selection.method = "vst", nfeatures = 1000)     #找高变基因，让后进行scale，有了sacle才能进行pca，umap、tsne可视化之类的后续
seurat_merge <- ScaleData(seurat_merge)   ##normalization一般是对文库处理，目的消除一些技术差异；scale一般对基因表达量处理（典型的z-score：表达量减均值再除以标准差），目的是后续分析不受极值影响。
seurat_merge <- RunPCA(seurat_merge)

# #ElbowPlot 肘图  (只针对pca)
# pdf(file = "Results/1.1.4_seurat_notharmony/1.1.4_ElbowPlot.pdf", width = 18, height = 6)
# # print(ElbowPlot(seurat_merge,reduction = 'scVI'))
# print(ElbowPlot(seurat_merge,reduction = 'pca'))
# dev.off()


# #瘦身,默认会把降维矩阵（pca）删掉，添加dimreduces参数,使pca保存
# # seurat_merge <- DietSeurat(seurat_merge,dimreducs = "scVI")
# # seurat_merge <- DietSeurat(seurat_merge,dimreducs = "pca")
# #可视化PCA（umap和tsne 默认的reduction是pca），根据肘图进行维度确认
# # seurat_merge <- RunUMAP(seurat_merge, dims = 1:10,reduction = 'scVI')
# #seurat_merge <- RunTSNE(seurat_merge, dims = 1:10,check_duplicates = FALSE,reduction = 'scVI')
seurat_merge <- RunUMAP(seurat_merge, dims = 1:10,reduction = 'pca')
# #seurat_merge <- RunTSNE(seurat_merge, dims = 1:10,check_duplicates = FALSE,reduction = 'pca')
# save(seurat_merge, file= paste0(outdir, '/1.1.4_seurat_merge.RData'))
# print(dim(seurat_merge))


# ##使用seurat_merge做umap图
# # load(paste0(outdir, '/1.1.4_seurat_merge.RData'))
# #找高变基因前20
# top20genes <- head(VariableFeatures(seurat_merge),20)
# #可视化基因表达差异程度（前20）
# pdf(file = "Results/1.1.4_seurat_notharmony/1.1.4_vstgenes.pdf", width = 18, height = 6)
# plot1 <- VariableFeaturePlot(seurat_merge)
# plot2 <- LabelPoints(plot = plot1, points = top20genes, repel = TRUE)
# print(plot1 + plot2)
# dev.off()
# #=========================
# # #查看差异效果是否显著
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
print(DimPlot(seurat_merge, reduction = "umap", pt.size = 0.1, group.by = "Patients", label = TRUE))
print(DimPlot(seurat_merge, reduction = "umap", pt.size = 0.1, group.by = "group", label = TRUE))
#print(DimPlot(seurat_merge, reduction = "tsne", pt.size = 0.1, group.by = "Patients", label = TRUE))
#print(DimPlot(seurat_merge, reduction = "tsne", pt.size = 0.1, group.by = "group", label = TRUE))
dev.off()
#
pdf(file = "Results/1.1.4_seurat_notharmony/1.1.4_RNA_featureplots.pdf", width = 12, height = 5)
FeaturePlot(seurat_merge, features = c("nFeature_RNA", "nCount_RNA"))
dev.off()


##resolution(分辨率越低，聚类的数量越少)
# load(paste0(outdir, '/1.1.4_seurat_merge.RData'))
# seurat_merge <- RunUMAP(seurat_merge, dims = 1:10,reduction="pca")
# seurat_merge <- RunTSNE(seurat_merge, dims = 1:10,check_duplicates = FALSE,reduction="pca")
seurat_merge <- FindNeighbors(seurat_merge, dims = 1:10,reduction = "pca")
# seurat_merge <- RunUMAP(seurat_merge, dims = 1:10,reduction="scVI")
# seurat_merge <- RunTSNE(seurat_merge, dims = 1:10,check_duplicates = FALSE,reduction="scVI")
# seurat_merge <- FindNeighbors(seurat_merge, dims = 1:10,reduction = "scVI")
seurat_merge <- FindClusters(seurat_merge, resolution = seq(0.4,1.5,0.1),algorithm = 1)


#cluster_umap_Plot，聚类后的umap图
pdf(file = "Results/1.1.4_seurat_notharmony/1.1.4_umap_cluster.pdf",width = 7, height = 6)
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

save(seurat_merge, file = paste0(outdir, '/1.1.4_seurat_merge.RData'))

Fibroblast_cell <- list(fibroblast = c("PDGFRA", "ACTA2", "SULF1", "CTGF", "TAGLN","TPM1", "GINS1" ,"THY1", "RBP1", "COL1A2", "COL1A1", "C1R", "IGFBP7", "SFRP2", "MGP","C1S", "DCN", "CXCL14", "COL3A1","LUM"),
                       endothelial = c("VWF","PECAM1", "CD31", "CD34", "HSPG2", "LDB2", "GPR116", "PTPRB", "DOCK9", "CDH5", "SELE"),
                       muscle = c("ACTA2","MYH11","NR2F2","CRYAB","TPM2","GUCA2B","GUCA2A","LMOD1","TPPP3","TAGLN"),
                       MSC = c("CD44", "ITGA1","NT5E","THY1"),
                       mycaf = c("ACTA2","VIM","CCN2", "COL1A1","COL5A1","COL6A1","TNC","TGFB1","THY1","TAGLN","COL12A1","PDGFRB"), #TGFβ/SMAD2/3 
                       icaf = c("IL1A","IL1B","IL6","IL11","LIF","CLEC3B","COL14A1","GSN","LY6C1","CXCL12","CXCL14","PDGFRA"), ##IL-1/JAK-STAT3
                       apcaf = c("SLPI","SAA3","CD74","H2-Ab1","NKAIN4", "IRF5"),
                       psc = c("DES", "GFAP", "CHRNA1"),
                       PVL = c("MCAM", "CD146", "ACTA2", "PDGFRB", "LYVE1","RGS5","COL4A1","CALD1","SPARCL1","SPARC","NID1"),
                       Pericyte = c("KCNJ8","ACTA2","PDGFRB","MCAM","HIGD1B","ANGPT2","VIM","MFGE8","MYO1B","NOTCH3","COX4I2"),
                       Cyclegene = c("MKI67","CDC20", "CENPF","PTTG1","TOP2A","CCNB1","PCNA")
)
setwd("/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC_magic/Futher_Legend/Mesenchymal/FeaturePlot/notharmony")
n <- 1
while( as.numeric(n) <= 11 ){
  a <- unlist(Fibroblast_cell[n])
  # a <- PC_marker[n]
  b <- FeaturePlot(seurat_merge, features =a)
  c <- paste0(names(Fibroblast_cell)[n],".png")
  ggsave(b,filename = c,width = 25,height = 25,device = "png")
  n= n+1
}

#之后寻找每个cluster的差异基因，之后根据差异基因进行确定亚型，还有featurePlot进行画图
# #group_by就是按照变量进行分组，ungroup就是取消分组信息，回到默认
# seurat_merge <- FindClusters(seurat_merge, resolution = 1.1,algorithm = 1)
# seurat_merge_marker <- FindAllMarkers(seurat_merge, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# seurat_merge_marker %>%
#             group_by(cluster) %>%
#             dplyr::filter(avg_log2FC > 1) %>%
#             arrange(desc(avg_log2FC)) %>%
#             slice_head(n = 20) %>%
#             ungroup() -> top20
# write.table(seurat_merge_marker, file =paste0(outdir, "seurat_merge_marker.txt") , sep = "\t", row.names = TRUE)  #file后面写入保存路径
# save(seurat_merge_marker,file= paste0(outdir, 'seurat_merge_marker.RData'))

# seurat_merge_marker <- read.table(paste0(outdir, "seurat_merge_marker.txt"))
# #根据Marker gene 画热图，看那几个cluster是不是一类，根据表达情况
# pdf(file = "Results/1.1.4_seurat_notharmony/cluster_Marker_heatmap.pdf",width = 14, height = 18)
# print(DoHeatmap(seurat_merge, features = top20$gene))
# dev.off()


#聚类后的feature PLot看不出来的话 就单独看簇的marker gene
# seurat_merge0_marker <- FindMarkers(seurat_merge,ident.1 =0)
# write.table(seurat_merge0_marker, file = "/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC/Data/seurat_harmony/seurat_merge_marker_0.txt", sep = "\t", row.names = TRUE) 
# head(seurat_merge0_marker, n = 20)

#==================================去批次======================= 
if(!dir.exists(paste0(workdir, "/Results/1.1.4_seurat_harmony"))){dir.create(paste0(workdir, "/Results/1.1.4_seurat_harmony"), recursive = TRUE)}
if(!dir.exists(paste0(workdir, "/Data/seurat_harmony"))){dir.create(paste0(workdir, "/Data/seurat_harmony"), recursive = TRUE)}
outdir  = paste0(topdir, "/Data/seurat_harmony")
setwd(workdir)

#读取数据
# load('/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC/Data/seurat_notharmony/1.1.4_seurat_merge.RData')

#去批次
seurat_merge <- seurat_merge %>% RunHarmony("Patients", plot_convergence = FALSE,reduction = "pca") 
#然后后面所有的reduction的pca全部换成harmony，除了后面的umap或tsne不换以外。而且下面两行的run~也加上reduction+harmony。
seurat_merge <- RunUMAP(seurat_merge, dims = 1:10,reduction = "harmony")
seurat_merge <- RunTSNE(seurat_merge, dims = 1:10,check_duplicates = FALSE,reduction="harmony")
save(seurat_merge, file = paste0(outdir,'/1.1.4_seurat_merge_harmony.RData'))


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

save(seurat_merge, file = paste0(outdir,"/1.1.4_seurat_merge_harmony.RData"))
Immune_cell <- list(immuneCells= c("PTPRC"),
                      macrophage = c("CD163","ADGRE1","CD163", "SLC11A1", "APOC1", "CD86", "CSF1R", "SLCO2B1", "CD68", "F13A1", "CD14", "AIF1", "CD80","FCER1G", "FCGR3A", "TYROBP"),
                       Tcell = c("CD3E","CD2", "CD3D",  "CD3G","CD4", "CD8A", "CD8B", "GZMK", "FOXP3","TRDC","NKG7","CD79A"),
                       Treg = c("FOXP3","CTLA4","MTNFRSF4","IRF4", "BATF", "TNFRSF18","TOX2","PRDMI"),
                       CD4NV_CM_rest = c("LEF1","ATM","SELL","KLF2","ITGA6"),
                      CD4CD8_rest = c("IL7R","CD52","S100A4","TGFB3","AQP3","NLRP3","KLF2","ITGB7"),
                      T_IFNresp = c("IFIT3", "IFIT2","STAT1", "MX1","IRF7", "ISG15","IFITM3","OAS2","JAK2","SOCS1","TRIM2"),
                      T_prof = c("LIF", "IL2","CENPV","NME1","FABP5","ORC6","G0S2","GCK"),
                      T_cytotoxic = c("CCL5","GZMK","GNLY","EOMES","ZNF683","KLRG1","NKG7","ZEB2"),
                      T_cytokine = c("CCL3", "IFNG","CCL4", "XCL1","XCL2", "CSF2","IL10","HOPX","TIM3", "LAG3","PRF1","TNFRSF9","NKG7","IL26"),
                      CD4conv = c("LTB", "SPOCK2","SOCS3","PBXIP1"),
                      Treg = c("FOXP3","TNFRSF4","CARD16","IL2RA","TBC1D4","CD25"), ## IL2RA CD25
                      CD8 = c("CCL4", "CST7", "GZMA", "GZMH"),
                      Texh = c("CTLA4", "LAG3", "PTMS", "PDCD1", "TIGIT"),
                      Tmem = c("CD44", "GZMK","IL7R", "CD69", "CD27"),
                      Tresid = c("SELL","ITGAE"),
                      Teff = c("CTSW","PRF1","GZMB","FGFBP2", "FCGR3A","TRGC2","GNLY","TYROBP", "IFNG", "CXCL13"),
                      Tgd = c("TYROBP","FCER1G","CMC1","KLRD1","KLRG1","KLRF1","TRDC","CTSW","PRF1","GZMB","FGFBP2", "FCGR3A","GNLY"),
                      NK = c("KLRB1","KLRC1", "GNLY", "XCL1", "AREG", "FCER1G","FCGR3A","KLRF1", "PRF1"),
                      NKT = c("NKG7","FCER1G","FCGR3A","GZMA","KLRB1","KLRC1","GZMB", "NCAM1"), 
                      ILC = c("SOX4", "KLRB1","BCL6","LST1","RORC"),
                       Bcell = c("KIF4A","FANCI","CHAF1A", "GTSE1", "ASPM", "SPC25","NCAPG2","POLA2","NCAPD3", "CD19", "CD21", "MS4A1", "CD79A", "CD79B", "BLNK"),
                       plasmocyte = c("IGLC2","IGLC3","IGHG3","IGHG1","IGHG4","JCHAIN","IGLC1"),
                       DCs = c("HLA-DRB1","ITGAX","CD83", "HLA-DMA","HLA-DQB1", "HLA-DPA1","HLA-DPB1","AIF1","LST1","FTL","HLA-C"),
                       NKcell = c("GNLY", "FCER1G", "KLRB1", "KLRC1", "AREG", "XCL1"),
                       mastCells = c("RHOH","BTK","FER","GATA2","IL4R","KIT","LCP2","ENPP3","RAC2","LAT2","CD84","LAT","ADGRE2","UNC13D","NDEL1", 
                                     "CPA3","TPSB2","TPSAB1","MS4A2","SLC18A2","IL1RL1","PTGS1","HPGDS"),
                       neutrophils = c("CSF3R","S100A8", "S100A9", "LTF", 'CEBPE', 'CEBPA', 'ETV6', 'FOXP1', 'GFI1', 'SPI1', 'STAT3', 'ELANE', 'GSTM1', 'LCN2', 'LYZ2', 'MPO', 'PRTN3')
)

setwd("/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC_magic/Futher_Legend/Immune/FeaturePlot/harmony")
n <- 1
while( as.numeric(n) <= 27 ){
  a <- unlist(Immune_cell[n])
  # a <- PC_marker[n]
  b <- FeaturePlot(seurat_merge, features =a)
  c <- paste0(names(Immune_cell)[n],".png")
  ggsave(b,filename = c,width = 25,height = 25,device = "png")
  n= n+1
}

#之后寻找每个cluster的差异基因，之后根据差异基因进行确定亚型，还有featurePlot进行画图
#group_by就是按照变量进行分组，ungroup就是取消分组信息，回到默认
# # setwd("/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC")
# seurat_merge <- FindClusters(seurat_merge, resolution = 0.7,algorithm = 1)
# seurat_merge_marker <- FindAllMarkers(seurat_merge, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# seurat_merge_marker %>%
#             group_by(cluster) %>%
#             dplyr::filter(avg_log2FC > 1) %>%
#             arrange(desc(avg_log2FC)) %>%
#             slice_head(n = 20) %>%
#             ungroup() -> top20
# write.table(seurat_merge_marker, file = paste0(outdir,"seurat_merge_marker.txt"), sep = "\t", row.names = TRUE)  #file后面写入保存路径
# save(seurat_merge_marker,file= paste0(outdir,'seurat_merge_marker.RData'))
# write.table(top20, file =paste0(outdir, "top20.txt"), sep = "\t", row.names = TRUE) 


# #根据Marker gene 画热图，看那几个cluster是不是一类，根据表达情况
# # seurat_merge_marker <- read.table('/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC/Data/seurat_harmony/seurat_merge_marker.txt')
# pdf(file = "Results/1.1.4_seurat_harmony/cluster_Marker_heatmap.pdf",width = 14, height = 18)
# print(DoHeatmap(seurat_merge, features = top20$gene))
# dev.off()







