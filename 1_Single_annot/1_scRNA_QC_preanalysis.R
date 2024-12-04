#!/usr/bin/Rscript
#nohup sh work.sh &
#qsub -cwd -l vf=200g,p=40 -P P21H28400N0232_super -q st.q,st_supermem.q spaGCN.sh

# if (length(args) != 11) {
#   stop ("Usage: Rscript scRNA_QC_preanalysis.R <threads> <maxGenes> <id> <minUMI> <minGene> <maxMt> <pc> <res> <minDoublets> <out_dir> <tumor type>\n    
#         }
args=commandArgs(T)
if (length(args) != 11) {
stop (" Rscript /jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/code/Single_cell/codes/1_scRNA_QC_preanalysis.R 10 8000 R231011001 1000 200 20 30 0.5 0.075 /jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/BC/sample_2.1/  BC")

}
#/mount/WD12T1/fanfengpu/codes/1_scRNA_QC_preanalysis.R
### Seurat clustering
library(Seurat)
library(dplyr)
library(patchwork)
library(DoubletFinder)
library(ggplot2)
library(future)
library(Matrix)
library(Cairo)
options(bitmapType='cairo')

maxGenes <- as.numeric(args[2]) #8000
SampleID <- args[3] #R
minUMIs <- as.numeric(args[4]) #1000
minGenes <- as.numeric(args[5]) #200
maxPercent.mt <- as.numeric(args[6]) #10
dim.usage <- as.numeric(args[7]) #20
res.usage <- as.numeric(args[8]) #0.5
doublets.percentage <- as.numeric(args[9]) #0.05
workdir <- args[10] #workdir
organs <- args[11]

source("/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/code/Single_cell/codes/scRNA_primary.R")
source("/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/code/Single_cell/codes/colors.R")
source("/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/code/Single_cell/codes/signatureGenes.R") #human signature gene 
#source("/jdfsbjcas1/ST_BJ/P21H28400N0232/wenfeng/project/codes/signatureGene_mouse.R")
source("/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/code/Single_cell/codes/0_initNEWanalysis.R")
source("/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/code/Single_cell/codes/DEG_wilcox.R")

setwd(workdir)
createObjdir(workdir = workdir)

plan()#和apply下，进行并行运行
plan("multisession", workers = as.numeric(args[1]))
plan()
options(future.globals.maxSize= 429496729600)
#列出路径下的文件名
sampList <- list.dirs(paste0(workdir, "/rawData"), recursive = FALSE, full.names = FALSE)#recursive如果为 TRUE，则递归地列出所有子目录;fullnames为 TRUE:返回完整的路径,False:返回相对路径

if(length(sampList) == 1){
  mat1 <- readRawMatrix_sc(matrix_dir = paste0(workdir, "/rawData/",sampList[1],"/FilterMatrix/"))
  colnames(mat1) <- paste0(sampList[1], "_", colnames(mat1))
  mat_comb <- mat1
  metadata_comb <- sapply(colnames(mat_comb), function(x) unlist(strsplit(x,"_", fixed = TRUE))[1])
  metadata_comb <- data.frame(group = metadata_comb)
  rownames(metadata_comb) <- colnames(mat_comb)
  save(mat1, metadata_comb, file = paste0(workdir,"/saveData/RawMat.RData"))
  rm(mat1)
  
}else{mat <- list()
for(i in 1:length(sampList)){
  mat[[i]] <- readRawMatrix_sc(matrix_dir = paste0(workdir, "rawData/",sampList[i],"/FilterMatrix/"))
  colnames(mat[[i]]) <- paste0(sampList[i], "_", colnames(mat[[i]]))
}
commongenes <- rownames(mat[[1]]) 
  for(i in 2:length(sampList)){
    commongenes <- intersect(commongenes, rownames(mat[[i]]))
  }
matfilt <- list()
for(i in 1:length(sampList)){
  matfilt[[i]] <- mat[[i]][commongenes,]
}
  
  mat_comb <- Reduce(cbind, matfilt)
  
  metadata_comb <- sapply(colnames(mat_comb), function(x) unlist(strsplit(x,"_CELL", fixed = TRUE))[1])
  metadata_comb <- data.frame(group = metadata_comb)
  rownames(metadata_comb) <- colnames(mat_comb)
  
  save(mat, mat_comb, metadata_comb, file = paste0(workdir,"/saveData/RawMat.RData"))
  rm(mat,matfilt)
  
}

seurat_comb <- CreateSeuratObject(counts = as.matrix(mat_comb), meta.data = metadata_comb, names.field = FALSE,project = SampleID, min.cells = 1, min.features = 1)

seurat_comb[["percent.mt"]] <- PercentageFeatureSet(seurat_comb, pattern = "^MT-")
cat("1.raw matrix:",dim(seurat_comb),"\n")


#### plot the qc plots ####
#VlnPlot(seurat_comb, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "orig.ident")
plot1 <- FeatureScatter(seurat_comb, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_comb, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

plot3 <- FeatureScatter(seurat_comb, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "group")
plot4 <- FeatureScatter(seurat_comb, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "group")

if(!is.null(dev.list())){dev.off()} #是否有任何打开的图形设备

pdf(file = paste0(workdir, "/Results/Figures/1_qc.pdf"), width = 10, height = 6)
print(VlnPlot(seurat_comb, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "group"))
print(VlnPlot(seurat_comb, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = NULL))
print(plot1 + plot2)
print(plot3 + plot4)
dev.off()

qc_stat <- rbind(summary(seurat_comb@meta.data$nCount_RNA),summary(seurat_comb@meta.data$nFeature_RNA),summary(seurat_comb@meta.data$percent.mt)) #按行进行合并形成一个矩阵，行名和列名分别是[.,.]
rownames(qc_stat) <- c("nUMIs","nGenes","percent.mt")
write.table(qc_stat, file = paste0(workdir,"/Results/", SampleID,"_raw_qc_stat.txt"), quote = F, sep = "\t",row.names = T,col.names=NA)


#### filter the data QC  ####
seurat_comb_filter <- subset(seurat_comb, subset = nCount_RNA>minUMIs & nFeature_RNA > minGenes & nFeature_RNA < maxGenes & percent.mt < maxPercent.mt)

cat("2.Filter Low quality cell matrix:",dim(seurat_comb_filter),"\n")

#### plot the plots of post plot  ####
plot1 <- FeatureScatter(seurat_comb_filter, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "orig.ident")
plot2 <- FeatureScatter(seurat_comb_filter, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident")

plot3 <- FeatureScatter(seurat_comb_filter, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "group")
plot4 <- FeatureScatter(seurat_comb_filter, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "group")

if(!is.null(dev.list())){dev.off()}

pdf(file = paste0(workdir, "/Results/Figures/1_qc_post.pdf"), width = 10, height = 6)
print(VlnPlot(seurat_comb_filter, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "group"))
print(VlnPlot(seurat_comb_filter, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "orig.ident"))

print(plot1 + plot2)
print(plot3 + plot4)

dev.off()

#### statistic the QC results  ####
qc <- c(ncol(seurat_comb), nrow(seurat_comb), median(seurat_comb$nFeature_RNA), median(seurat_comb$nCount_RNA), median(seurat_comb$percent.mt),
        mean(seurat_comb$nFeature_RNA), mean(seurat_comb$nCount_RNA), mean(seurat_comb$percent.mt))
qcPost <- c(ncol(seurat_comb_filter), nrow(seurat_comb_filter), median(seurat_comb_filter$nFeature_RNA), median(seurat_comb_filter$nCount_RNA), 
            median(seurat_comb_filter$percent.mt),mean(seurat_comb_filter$nFeature_RNA), mean(seurat_comb_filter$nCount_RNA), mean(seurat_comb_filter$percent.mt))
stat <- data.frame(qcBefore = qc, qcPost = qcPost)
rownames(stat) <- c("cells", "genes", "medianFeatures", "medianCounts", "medianMT", "meanFeatures", "meanCounts", "meanMT")
stat$percent <- 100*qcPost/qc
write.csv(stat, file = paste0(workdir, "/Results/QCstat.csv"), quote = FALSE)

save(seurat_comb, file = paste0(workdir, "/saveData/seurat_", SampleID, ".RData"))
rm(seurat_comb)

qc_stat <- rbind(summary(seurat_comb_filter@meta.data$nCount_RNA),summary(seurat_comb_filter@meta.data$nFeature_RNA),summary(seurat_comb_filter@meta.data$percent.mt))
rownames(qc_stat) <- c("nUMIs","nGenes","percent.mt")
write.table(qc_stat, file = paste0(workdir,"/Results/", SampleID,"_post_qc_stat.txt"), quote = F, sep = "\t",row.names = T,col.names=NA)

#### normalize and reduce the data ####
seurat_comb_filter <- NormalizeData(seurat_comb_filter)
seurat_comb_filter <- FindVariableFeatures(seurat_comb_filter, selection.method = "vst", nfeatures = 2000)

top10genes <- head(VariableFeatures(seurat_comb_filter),10)
plot1 <- VariableFeaturePlot(seurat_comb_filter)
plot2 <- LabelPoints(plot = plot1, points = top10genes, repel = TRUE)

pdf(file = paste0(workdir, "/Results/Figures/2_vstgenes.pdf"), width = 18, height = 6)
print(plot1 + plot2)
dev.off()

all.genes <- rownames(seurat_comb_filter)
seurat_comb_filter <- ScaleData(seurat_comb_filter, features = all.genes)

seurat_comb_filter <- RunPCA(seurat_comb_filter, features = VariableFeatures(object = seurat_comb_filter))

pdf(file = paste0(workdir, "/Results/Figures/2_pca_reduce.pdf"), width = 8, height = 6)
print(DimPlot(seurat_comb_filter, reduction = "pca", group.by = "orig.ident"))
print(DimPlot(seurat_comb_filter, reduction = "pca", group.by = "group"))
print(DimHeatmap(seurat_comb_filter,dims = 1:10, cells = 500, balanced = TRUE))
dev.off()

pdf(file = paste0(workdir, "/Results/Figures/2_ElbowPlot.pdf"), width = 18, height = 6)
print(ElbowPlot(seurat_comb_filter))
dev.off()

## umap and tsne reduction
seurat_comb_filter <- RunUMAP(seurat_comb_filter, dims = 1:dim.usage)
seurat_comb_filter <- RunTSNE(seurat_comb_filter, dims = 1:dim.usage)

pdf(file = paste0(workdir, "/Results/Figures/3_umap_tsne.pdf"),width = 5, height = 4)
print(DimPlot(seurat_comb_filter, reduction = "umap", pt.size = 0.2, group.by = "orig.ident", label = TRUE))
print(DimPlot(seurat_comb_filter, reduction = "umap", pt.size = 0.2, group.by = "group", label = TRUE))
print(DimPlot(seurat_comb_filter, reduction = "tsne", pt.size = 0.2, group.by = "orig.ident", label = TRUE))
print(DimPlot(seurat_comb_filter, reduction = "tsne", pt.size = 0.2, group.by = "group", label = TRUE))
dev.off()

pdf(file = paste0(workdir, "/Results/Figures/3_RNA_featureplots.pdf"), width = 11, height = 5)
FeaturePlot(seurat_comb_filter, features = c("nFeature_RNA", "nCount_RNA"))
dev.off()

cat("3. reduce the data:",dim(seurat_comb_filter),"\n")

#### to remove the doublets  ####
library(tidyverse)
library(Seurat)
library(spam)
library(fields)
library(DoubletFinder)
Find_doublet <- function(data, doublets.percentage=0.05, dim.usage=30){
  #表示不使用 SCT 标准化方法,而是使用默认的对数标准化方法。
  sweep.res.list <- paramSweep_v3(data, PCs = 1:dim.usage, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE) #没有真实标签(ground truth)可用
  bcmvn <- find.pK(sweep.stats) ### output plot,pK定义了用于计算 pANN 的 PC 邻域大小
  nExp_poi <- round(doublets.percentage*ncol(data))
  p <- as.numeric(as.vector(bcmvn[bcmvn$MeanBC==max(bcmvn$MeanBC),]$pK)) ### pK Selection,MeanBC 是一个评估指标,它越大表示 pK 值越合适
  #pN 用于计算 pANN 的邻域大小,reuse.. 不重复使用之前计算的 pANN 值
  data <- doubletFinder_v3(data, PCs = 1:dim.usage, pN = 0.25, pK = p, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE) ### output plot
  colnames(data@meta.data)[ncol(data@meta.data)] = "doublet_info"
  return(data)
}

pdf(file = paste0(workdir,"/Results/Figures/4_find_doublet.pdf"))
seurat_comb_filter <- Find_doublet(seurat_comb_filter, doublets.percentage = doublets.percentage)
DimPlot(seurat_comb_filter, reduction = "umap", group.by = "doublet_info")
dev.off()

write.table(seurat_comb_filter@meta.data,paste0(workdir, "/Results/",SampleID,"_doublets_info.txt"),quote = F, sep = "\t",row.names = T,col.names=NA)

cat("4. remove the doublets ",table(seurat_comb_filter@meta.data$doublet_info),"\n")
#### to analysis with the singlets ####
seurat_comb_singlet<- subset(seurat_comb_filter,subset=doublet_info=="Singlet")  #去除双胞

save(seurat_comb_filter, file = paste0(workdir, "/saveData/seurat_", SampleID, "_filter.RData"))
rm(seurat_comb_filter)

#### normalize and reduce the data ####
seurat_comb_singlet <- NormalizeData(seurat_comb_singlet)
seurat_comb_singlet <- FindVariableFeatures(seurat_comb_singlet, selection.method = "vst", nfeatures = 2000)

top10genes <- head(VariableFeatures(seurat_comb_singlet),10)
plot1 <- VariableFeaturePlot(seurat_comb_singlet)
plot2 <- LabelPoints(plot = plot1, points = top10genes, repel = TRUE)

pdf(file = paste0(workdir, "/Results/Figures/5_vstgenes_singlet.pdf"), width = 18, height = 6)
print(plot1 + plot2)
dev.off()

all.genes <- rownames(seurat_comb_singlet)
seurat_comb_singlet <- ScaleData(seurat_comb_singlet, features = all.genes)

seurat_comb_singlet <- RunPCA(seurat_comb_singlet, features = VariableFeatures(object = seurat_comb_singlet))


pdf(file = paste0(workdir, "/Results/Figures/5_pca_reduce_singlet.pdf"), width = 8, height = 6)
print(DimPlot(seurat_comb_singlet, reduction = "pca", group.by = "orig.ident"))
print(DimPlot(seurat_comb_singlet, reduction = "pca", group.by = "group"))
print(DimHeatmap(seurat_comb_singlet,dims = 1:10, cells = 500, balanced = TRUE))
dev.off()

## umap and tsne reduction
seurat_comb_singlet <- RunUMAP(seurat_comb_singlet, dims = 1:dim.usage)
seurat_comb_singlet <- RunTSNE(seurat_comb_singlet, dims = 1:dim.usage)

pdf(file = paste0(workdir, "/Results/Figures/5_RNA_featureplots_singlet.pdf"), width = 11, height = 5)
FeaturePlot(seurat_comb_singlet,reduction = "umap", features = c("nFeature_RNA", "nCount_RNA"))
dev.off()

pdf(file = paste0(workdir, "/Results/Figures/5_umap_tsne_singlet.pdf"),width = 5, height = 4)
print(DimPlot(seurat_comb_singlet, reduction = "umap", pt.size = 0.2, group.by = "orig.ident", label = TRUE))
print(DimPlot(seurat_comb_singlet, reduction = "umap", pt.size = 0.2, group.by = "group", label = TRUE))
print(DimPlot(seurat_comb_singlet, reduction = "tsne", pt.size = 0.2, group.by = "orig.ident", label = TRUE))
print(DimPlot(seurat_comb_singlet, reduction = "tsne", pt.size = 0.2, group.by = "group", label = TRUE))
dev.off()

cat("5. reduce the singlets ",table(seurat_comb_singlet@meta.data$doublet_info),"\n")
#### cluster all the cells  ####
seurat_comb_singlet <- FindNeighbors(seurat_comb_singlet, dims = 1:dim.usage)

seurat_comb_singlet <- FindClusters(seurat_comb_singlet, resolution = 0.4)
seurat_comb_singlet <- FindClusters(seurat_comb_singlet, resolution = 0.5)
seurat_comb_singlet <- FindClusters(seurat_comb_singlet, resolution = 0.8)
seurat_comb_singlet <- FindClusters(seurat_comb_singlet, resolution = 1)
seurat_comb_singlet <- FindClusters(seurat_comb_singlet, resolution = 1.2)


pdf(file = paste0(workdir, "/Results/Figures/6_umap_cluster_singlet.pdf"),width = 7, height = 6)
print(DimPlot(seurat_comb_singlet, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.4", label = TRUE))
print(DimPlot(seurat_comb_singlet, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.5", label = TRUE))
print(DimPlot(seurat_comb_singlet, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.8", label = TRUE))
print(DimPlot(seurat_comb_singlet, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.1", label = TRUE))
print(DimPlot(seurat_comb_singlet, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.1.2", label = TRUE))
dev.off()

cat("6. cluster the singlets ",table(seurat_comb_singlet@meta.data$doublet_info),"\n")
save(seurat_comb_singlet, file = paste0(workdir, "/saveData/seurat_", SampleID, "_single.RData"))

#### cluster find DEGs, and marker genes plots  ####
seurat_comb_singlet_markers <- FindAllMarkers(seurat_comb_singlet, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

save(seurat_comb_singlet_markers, file = paste0(workdir, "/saveData/seurat_comb_singlet_markers.RData"))
seurat_comb_singlet_markers %>% group_by(cluster) %>% top_n(n=10, wt = avg_log2FC) -> top10

height <- 0.8 + nrow(top10)*0.11
pdf(file = paste0(workdir, "Results/Figures/7_ClusterMarkers_heatmap_singlet.pdf"), width = 10, height = height)
print(DoHeatmap(seurat_comb_singlet, features = top10$gene) + NoLegend())
dev.off()
###
genemarker <- SigGeneral_all
genemarker <- unique(unlist(genemarker)) %>% as.character()
genemarker <- genemarker[genemarker %in% rownames(seurat_comb_singlet)]

# hei <- ceiling(length(strommarker)/4) * 3
# pdf(file = paste0(workdir, "Results/Figures/4_strommarker_featureplots.pdf"), width = 14.5, height = hei)
# print(FeaturePlot(seurat_comb_singlet, features = strommarker, reduction = "umap", ncol = 4))
# dev.off()
# hei <- ceiling(length(immunemarker)/4) * 3
# pdf(file = paste0(workdir, "Results/Figures/4_immunemarker_featureplots.pdf"), width = 14.5, height = hei)
# print(FeaturePlot(seurat_comb_singlet, features = immunemarker, reduction = "umap", ncol = 4))
# dev.off()
hei <- ceiling(length(genemarker)/4) * 3 #向上取整
pdf(file = paste0(workdir, "/Results/Figures/7_genemarker_featureplots.pdf"), width = 14.5, height = hei)
print(FeaturePlot(seurat_comb_singlet, features = genemarker, reduction = "umap", ncol = 4))
dev.off()

wid <- length(unique(Idents(seurat_comb_singlet)))*0.3 + 1

# hei <- length(strommarker)*0.3 + 1
# pdf(file = paste0(workdir, "Results/Figures/4_strommarker_vlnplots.pdf"), width = wid, height = hei)
# print(VlnPlot(seurat_comb_singlet, features = strommarker, log = F, stack = TRUE, flip = TRUE) + NoLegend())
# dev.off()
# hei <- length(immunemarker)*0.3 + 1
# pdf(file = paste0(workdir, "Results/Figures/4_immunemarker_vlnplots.pdf"), width = wid, height = hei)
# print(VlnPlot(seurat_comb_singlet, features = immunemarker, log = F, stack = TRUE, flip = TRUE) + NoLegend())
# dev.off()
hei <- length(genemarker)*0.3 + 1
sc_names = rownames(seurat_comb_singlet)[rowSums(seurat_comb_singlet) > 0]
pdf(file = paste0(workdir, "/Results/Figures/7_genemarker_vlnplots.pdf"),, width = wid, height = hei)
print(VlnPlot(seurat_comb_singlet, features = intersect(genemarker, sc_names), log = F, stack = TRUE, flip = TRUE) + NoLegend())
dev.off()

cat("7. cluster cell type enrichment" ,table(seurat_comb_singlet@meta.data$doublet_info),"\n")
# 
# t1 <- Sys.time()
# seurat_comb_singlet_sigmarker_av <- Seurat_cluster_wilcoxRankEnrichment(seurat_obj = seurat_comb_singlet, genesets_len = 1, rank_meth = "averag", workdir = workdir)
# t2 <- Sys.time()
# t2-t1
# seurat_comb_singlet_sigmarker_fc <- Seurat_cluster_wilcoxRankEnrichment(seurat_obj = seurat_comb_singlet, genesets_len = 1, rank_meth = "FCrank", workdir = workdir)
# t3 <- Sys.time()
# t3-t2
# save(seurat_comb_singlet_sigmarker_av, file = paste0(workdir, "saveData/seurat_comb_singlet_sigmarker_av.RData"))
# save(seurat_comb_singlet_sigmarker_fc, file = paste0(workdir, "saveData/seurat_comb_singlet_sigmarker_fc.RData"))
# 
# pdf(paste0(workdir, "/Results/Figures/6_celltype_averag_enrichment.pdf"))
# ggplot(seurat_comb_singlet_sigmarker_av$Enrichment_all, aes(clust, Description)) + geom_point(aes(size=NES, fill=p.Val), alpha=0.9, pch=21, col="grey25") + theme_classic() +
#   scale_fill_gradient2(low = "white", high = "red") + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, colour = "black"),
#                                                             axis.text.y = element_text(colour = "black"))
# dev.off()
# 
# pdf(paste0(workdir, "/Results/Figures/6_celltype_FCrank_enrichment.pdf"))
# ggplot(seurat_comb_singlet_sigmarker_fc$Enrichment_all, aes(clust, Description)) + geom_point(aes(size=NES, fill=p.Val), alpha=0.9, pch=21, col="grey25") + theme_classic() +
#   scale_fill_gradient2(low = "white", high = "red") + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, colour = "black"),
#                                                             axis.text.y = element_text(colour = "black"))
# dev.off()



cat("8. finished" ,table(seurat_comb_singlet@meta.data$doublet_info),"\n")

#### to enrichment the clusters







