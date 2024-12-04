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

topdir = "/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/BC/R28BC"
workdir = topdir
setwd(workdir)


#多线程运行
plan()
plan("multisession", workers = as.numeric(2))
# Enable parallelization
plan()
options(future.globals.maxSize= 429496729600)

# raw Data
load('/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/BC/R28BC/seurat_merge_Immune.RData')
seurat_T <- subset(seurat_merge1,subset = cellTypes_new %in% c("Tcell"))
dim(seurat_T)

#标准化
seurat_T <- NormalizeData(seurat_T)
save(seurat_T, file = paste0(workdir, '/seurat_T.raw.RData'))

#load(paste0(workdir, '/seurat_merge3.raw.RData'))
if(!dir.exists(paste0(workdir, "/Results/Immune_futher"))){dir.create(paste0(workdir, "/Results/Immune_futher"), recursive = TRUE)}


#查看QC结果
pdf(file = paste0(workdir,"/Results/Immune_futher/1.1_qc_post.pdf"), width = 30, height = 6)
print(VlnPlot(seurat_T, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "Patients"))
dev.off()


#基于cluster类型进一步重新重新进行（标准化）——找高变——缩放——降维（PCA，scVI）——寻找临近距离——聚类
seurat_T <- FindVariableFeatures(seurat_T, selection.method = "vst", nfeatures = 1000)
seurat_T <- ScaleData(seurat_T)
seurat_T <- RunPCA(seurat_T)


#降维以后查看肘图 后续的可视化维度设置
pdf(file = paste0(workdir,"/Results/Immune_futher/ElbowPlot_1000_Epithe.pdf"), width = 18, height = 6)
print(ElbowPlot(seurat_T))
dev.off()

#可视化降维  umap和tsne 默认的reduction是pca，
# n <- 1
# while(n < 4){
#     a <- paste0("seurat_merge", n)
#     b <- get(as.name(a)) #b就是一个seurat对象了
# }   

seurat_T <- RunUMAP(seurat_T, dims = 1:10)
seurat_merge1 <- RunTSNE(seurat_merge1, dims = 1:10,check_duplicates = FALSE)
save(seurat_T, file='seurat_T_notharmony.RData')


# # #去批次 然后后面所有的reduction的pca全部换成harmony，除了后面的umap或tsne不换以外。而且下面两行的run~也加上reduction+harmony。
# seurat_merge3 <- seurat_merge3 %>% RunHarmony("Patients", plot_convergence = FALSE,reduction = "pca") 
# seurat_merge3 <- RunUMAP(seurat_merge3, dims = 1:10,reduction = "harmony")
# seurat_merge3 <- RunTSNE(seurat_merge1, dims = 1:10,check_duplicates = FALSE,reduction="harmony")
# save(seurat_merge3, file='seurat_merge3_harmony.RData')
# # print(dim(seurat_merge))

##使用seurat_merge做umap图
load('seurat_T_notharmony.RData')
#寻找在所有Tcell里面在每个样本的前20的高差异表达Gene
top20genes <- head(VariableFeatures(seurat_T),20)
pdf(file = "Results/Immune_futher/vstgenes.pdf", width = 18, height = 6)
plot1 <- VariableFeaturePlot(seurat_T)
plot2 <- LabelPoints(plot = plot1, points = top20genes, repel = TRUE)
print(plot1 + plot2)
dev.off()
# Plot：pca
pdf(file = "Results/Immune_futher/pca_reduce.pdf", width = 8, height = 6)
print(DimPlot(seurat_T, reduction="pca", group.by = "orig.ident"))
print(DimPlot(seurat_T, reduction="pca", group.by = "Patients"))
print(DimHeatmap(seurat_T,dims = 1:10, cells = 500, balanced = TRUE))
dev.off()


#Plot：umap_tsne_patients
pdf(file = "Results/Immune_futher/umap_tsne_patients.pdf", width = 9, height = 8)
print(DimPlot(seurat_T, reduction = "umap", pt.size = 0.1, group.by = "Patients", label = TRUE))
print(DimPlot(seurat_T, reduction = "umap", pt.size = 0.1, group.by = "group", label = TRUE))
#print(DimPlot(seurat_merge, reduction = "tsne", pt.size = 0.1, group.by = "Patients", label = TRUE))
#print(DimPlot(seurat_merge, reduction = "tsne", pt.size = 0.1, group.by = "group", label = TRUE))
dev.off()

pdf(file = "Results/Immune_futher/RNA_featureplots.pdf", width = 12, height = 5)
FeaturePlot(seurat_T, features = c("nFeature_RNA", "nCount_RNA"))
dev.off()
#pdf(file = "Results/1.1.4_seurat_notharmony/1.1.4_umap_tsne_patients_p4.pdf", width = 9, height = 8)
#print(DimPlot(seurat_merge, reduction = "umap", pt.size = 0.1, group.by = "Patients", label = TRUE, cols=c('P1T'='red', 'P22T'='green', 'P35T'='blue', 'P38T'='purple')))
#dev.off()

##resolution
seurat_T <- FindNeighbors(seurat_T, dims = 1:10,reduction = "pca")
seurat_T <- FindClusters(seurat_T, resolution = seq(0.5,1.3,0.1),algorithm = 1)

save(seurat_T, file = "seurat_T_cluster.RData")
#上皮1.1，其他两0.7
#查看一下聚类后的质量
pdf(file = "Results/Immune_futher/qc_post_Epithe.pdf", width = 30, height = 6)
plot1 <- FeatureScatter(seurat_T, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_T, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(seurat_T, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "Patients")
plot4 <- FeatureScatter(seurat_T, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "Patients")
print(VlnPlot(seurat_T, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "Patients"))
print(VlnPlot(seurat_T, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3),group.by ="orig.ident")
print(plot1 + plot2)
print(plot3 + plot4)
dev.off()

load('seurat_T_cluster.RData')
pdf(file = "Results/Immune_futher/umap_cluster.pdf",width = 7, height = 6)
#print(DimPlot(seurat_merge3, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.4", label = TRUE))
print(DimPlot(seurat_T, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.5", label = TRUE))
print(DimPlot(seurat_T, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.6", label = TRUE))
print(DimPlot(seurat_T, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.7", label = TRUE))
print(DimPlot(seurat_T, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.8", label = TRUE))
print(DimPlot(seurat_T, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.9", label = TRUE))
print(DimPlot(seurat_T, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.1",   label = TRUE))
print(DimPlot(seurat_T, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.1.1", label = TRUE))
print(DimPlot(seurat_T, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.1.2", label = TRUE))
print(DimPlot(seurat_T, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.1.3", label = TRUE))
# print(DimPlot(seurat_merge3, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.1.4", label = TRUE))
#print(DimPlot(seurat_merge3, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.1.5", label = TRUE))
dev.off()
save(seurat_T, file = "seurat_T_cluster.RData")


#Feature Plot
load('/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/BC/R28BC/Tcell/seurat_T_cluster.RData')
SigGeneral_all <- list(immuneCells= c("PTPRC"),
                    Tcell = c("CD2", "CD3D", "CD3E", "CD3G","CD4","CD8A","CD8B", "GZMK"),
                    CD4-C1 = c("FOXP3", "IL2RA", "IKZF2", "Icos", "RTKN2"),
                    CD4-C2 = c("CXCL13", "BTLA", "CTSB", "TMEM173"),
                    naive_T = c("IL7R", "CCR7", "TCF7","GPR183"),
                    CD8-C1 = c("HAVCR2", "LAG3", "RAB27A", "LYST", "IFNG"),
                    CD8-C2 = c("GZMK", "ITM2C", "CD52"),
                    CD8-C3 = c("IFIT3", "ISG15", "MX1", "IFI6", "IFIT1"),
                    NK-C1 = c("FCER1G", "TYROBP", "XCL1", "XCL2", "IL2RB"),
                    NK-C2 = c("KLRF1", "FCGR3A", "FGFBP2", "EFHD2", "SPON2")   
)

SigGeneral_all_1 <- list(immuneCells= c("PTPRC"),
                    Tcell= c("CD2", "CD3D", "CD3E", "CD3G","CD4","CD8A","CD8B", "GZMK"),
                    CD4_C1= c("FOXP3", "IL2RA", "IKZF2", "Icos", "RTKN2"),
                    CD4_C2= c("CXCL13", "BTLA", "CTSB", "TMEM173"),
                    naive_T= c("IL7R", "CCR7", "TCF7","GPR183"),
                    CD8_C1= c("HAVCR2", "LAG3", "RAB27A", "LYST", "IFNG"),
                    CD8_C2= c("GZMK", "ITM2C", "CD52"),
                    CD8_C3= c("IFIT3", "ISG15", "MX1", "IFI6", "IFIT1"),
                    NK_C1= c("FCER1G", "TYROBP", "XCL1", "XCL2", "IL2RB"),
                    NK_C2= c("KLRF1", "FCGR3A", "FGFBP2", "EFHD2", "SPON2"),
                    Treg = c("FOXP3","CTLA4","MTNFRSF4","IRF4", "BATF", "TNFRSF18","TOX2","PRDMI"),
                    CD4NV_CM_rest = c("LEF1","ATM","SELL","KLF2","ITGA6"),
                    CD4CD8_rest = c("IL7R","CD52","S100A4","TGFB3","AQP3","NLRP3","KLF2","ITGB7"),
                    T_IFNresp = c("IFIT3", "IFIT2","STAT1", "MX1","IRF7", "ISG15","IFITM3","OAS2","JAK2","SOCS1","TRIM2"),
                    T_prof = c("LIF", "IL2","CENPV","NME1","FABP5","ORC6","G0S2","GCK"),
                    T_cytotoxic = c("CCL5","GZMK","GNLY","EOMES","ZNF683","KLRG1","NKG7","ZEB2"),
                    T_cytokine = c("CCL3", "IFNG","CCL4", "XCL1","XCL2", "CSF2","IL10","HOPX","TIM3", "LAG3","PRF1","TNFRSF9","NKG7","IL26"),
                    CD8 = c("CCL4", "CST7", "GZMA", "GZMH"),
                    Texh = c("CTLA4", "LAG3", "PTMS", "PDCD1", "TIGIT"),
                    Tmem = c("CD44", "GZMK","IL7R", "CD69", "CD27"),
                    Tresid = c("SELL","ITGAE"),
                    Teff = c("CTSW","PRF1","GZMB","FGFBP2", "FCGR3A","TRGC2","GNLY","TYROBP", "IFNG", "CXCL13"),
                    Tgd = c("TYROBP","FCER1G","CMC1","KLRD1","KLRG1","KLRF1","TRDC","CTSW","PRF1","GZMB","FGFBP2", "FCGR3A","GNLY"),
                    NK = c("KLRB1","KLRC1", "GNLY", "XCL1", "AREG", "FCER1G","FCGR3A","KLRF1", "PRF1"),
                    NKT = c("NKG7","FCER1G","FCGR3A","GZMA","KLRB1","KLRC1","GZMB", "NCAM1"), ##FCGR3A cd16, NCAM1 cd56,
                    ILC = c("SOX4", "KLRB1","BCL6","LST1","RORC")
                    )
length(names(SigGeneral_all_1))
n <- 1
setwd("/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/BC/R28BC/Results/Immune_futher/Tcell/FeaturePlot")

while( as.numeric(n) < 27 ){
  a <- unlist(SigGeneral_all_1[n])
  b <- FeaturePlot(seurat_T, features =a)
  c <- paste0(names(SigGeneral_all_1)[n],".png")
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

seurat_merge1_marker <- FindAllMarkers(seurat_merge1, only.pos = TRUE)
seurat_merge1_marker %>%
            group_by(cluster) %>%
            dplyr::filter(avg_log2FC > 1) %>%
            slice_head(n = 10) %>%
            ungroup() -> top10
            
pdf(file = "Results/1.1.4_seurat_futher_cluser/cluster_1.1_Marker_heatmap.pdf",width = 15, height = 16)
print(DoHeatmap(seurat_merge1, features = top10$gene))
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
save(seurat_merge1, file = "seurat_merge_Immune.RData")


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

