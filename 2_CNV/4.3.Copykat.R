library(Seurat)
library(copykat)
library(tidyverse)
library(Matrix)
library(ggplot2)
rm(list = ls())

setwd("/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC/copykat")
 ##这儿我们直接导入Seurat标准化，聚类的pbmc数据
load('/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC/Data/seurat_harmony/seurat_merge_magic_harmony_subannot.RData')

seurat_merge_infercnv <- subset(seurat_merge, subset = cellTypes_new %in% c('Epithelial','Tcell', 'Bcell','Macrophage', 'Endothelial'))
seurat_lists = SplitObject(seurat_merge_infercnv, split.by = 'Patients')
save(seurat_lists, file = "/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC/copykat/1.4.1_seurat_merge_infercnv.RData")

load('/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC/copykat/1.4.1_seurat_merge_infercnv.RData')
seurat_i = seurat_lists[["HPCP1"]]
# current.cluster.ids <- c(0:8)
# new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")
# pbmc@meta.data$celltype <- plyr::mapvalues(x = pbmc@meta.data[,"seurat_clusters"], from = current.cluster.ids, to = new.cluster.ids)
# head(pbmc@meta.data)
scRNA <- seurat_i
counts <- as.matrix(scRNA@assays$RNA@counts)

copykat <- copykat(rawmat = counts, ngene.chr=5, sam.name="PC", n.cores=8)
# ngene.chr参数是过滤细胞的一个标准，它要求被用来做CNV预测的细胞，一个染色体上至少有5个基因。
# sam.name定义样本名称 (sample name)，会给出来的文件加前缀
save(copykat, "copykat.RData")
# SCLC_copykat_prediction.txt是上一步生成的结果

mallignant <- read.delim("PC_copykat_prediction.txt", row.names = 1) #参数指定将文件的第一列作为行名
# 把细胞的良恶性信息加入metadata
scRNA <- AddMetaData(scRNA, metadata = mallignant)
# pdf('pred_mallignant.pdf',width = 12, height = 5)
p1 <- DimPlot(scRNA, group.by = "cellTypes_new", label = T)
p2 <- DimPlot(scRNA, group.by = "copykat.pred") #+ scale_color_manual(values = c("red", "#17133F"."#17133F"))
pc <- p1 + p2
# dev.off()
ggsave("pred_mallignant.pdf", pc, width = 12, height = 5)


