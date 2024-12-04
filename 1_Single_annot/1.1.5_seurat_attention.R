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

#多线程分析
plan()
plan("multisession", workers = as.numeric(2))
# Enable parallelization
plan()
options(future.globals.maxSize= 429496729600)
#Feature Plot
length(names(SigGeneral_all_1))
n <- 1
setwd("/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/BC/R28BC/seurat_merge_FeaturePlot")  #FeaturePlot的保存路径

while( as.numeric(n) < 30 ){
  cellType <- unlist(SigGeneral_all_1[n])
  plot <- FeaturePlot(seurat_merge, features = cellType)
  name <- paste0(names(SigGeneral_all_1)[n],".png")
  ggsave(plot,filename = c,width = 25,height = 25,device = "png")
  n <- n + 1                                                         
}


#之后寻找每个cluster的差异基因，之后根据差异基因进行确定亚型，还有featurePlot进行画图
#group_by就是按照变量进行分组，ungroup就是取消分组信息，回到默认
setwd("/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC")
load('/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC_1/Data/seurat_merge_Immune_harmony.RData')
load('/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC/Data/seurat_harmony/1.1.4_seurat_merge_harmony.RData')
seurat_merge <- FindClusters(seurat_merge, resolution = 1.1,algorithm = 1)
seurat_merge_marker <- FindAllMarkers(seurat_merge, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
seurat_merge_marker %>%
            group_by(cluster) %>%
            dplyr::filter(avg_log2FC > 1) %>%
            arrange(desc(avg_log2FC)) %>%
            slice_head(n = 20) %>%
            ungroup() -> top20
#根据Marker gene 画热图，看那几个cluster是不是一类，根据表达情况
load('/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC/Data/seurat_harmony/seurat_merge_marker.RData')
seurat_merge_marker <- read.table('/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC/Data/seurat_harmony/seurat_merge_marker.txt')
pdf(file = "Results/1.1.4_seurat_harmony/cluster_Marker_heatmap_1.pdf",width = 14, height = 18)
print(DoHeatmap(seurat_merge, features = top20$gene))
dev.off()
write.table(seurat_merge_marker, file = "/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC/Data/seurat_harmony/seurat_merge_marker.txt", sep = "\t", row.names = TRUE)  #file后面写入保存路径
save(seurat_merge_marker,file= '/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC/Data/seurat_harmony/seurat_merge_marker.RData')

#聚类后的feature PLot看不出来的话 就单独看簇的marker gene
seurat_merge0_marker <- FindMarkers(seurat_merge,ident.1 =0)
write.table(seurat_merge0_marker, file = "/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC/Data/seurat_harmony/seurat_merge_marker_0.txt", sep = "\t", row.names = TRUE) 
head(seurat_merge0_marker, n = 20)
###改变gene类型
gene_n.df <- bitr(seurat_merge_marker[,7][1],
        fromType = "ENSEMBL",
        toType = "SYMBOL",
        OrgDb = org.Hs.eg.db
    )


#细胞注释
# seurat_merge <- FindVariableFeatures(seurat_merge, selection.method = "vst", nfeatures = 2000)     #找高变基因，让后进行scale，有了sacle才能进行pca，umap、tsne可视化之类的后续
# seurat_merge <- ScaleData(seurat_merge)   ##normalization一般是对文库处理，目的消除一些技术差异；scale一般对基因表达量处理（典型的z-score：表达量减均值再除以标准差），目的是后续分析不受极值影响。
# seurat_merge <- RunPCA(seurat_merge)
# seurat_merge <- FindNeighbors(seurat_merge, dims = 1:10,reduction = "pca")
seurat_merge <- FindClusters(seurat_merge, resolution =0.9)

cellTypes_new <- c("Epithelial","Epithelial","Epithelial","Epithelial","Epithelial","Epithelial","Epithelial","Epithelial","Mesenchymal","Mesenchymal","Epithelial","Mesenchymal","Epithelial","Epithelial","Mesenchymal","Immune","Immune","Mesenchymal","Epithelial","Epithelial","Mesenchymal","Epithelial","Epithelial")
length(cellTypes_new)
names(cellTypes_new) <- levels(seurat_merge)   #names是将函数名称提取出来
seurat_merge4 <- RenameIdents(seurat_merge4, cellTypes_new)  #通过 Idents 函数查看或修改。而 RenameIdents 函数则用于重命名聚类标识符
#新增一列，将细胞类型放入数据框新的一列中
seurat_merge4$cellTypes_new <- Idents(seurat_merge4)  #Ident() 函数是用于获取或设置每个细胞的聚类（cluster）标识符的函数

#细胞注释图
pdf(file = "/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/BC/R28BC/Results/1.1.4_seurat_notharmony/umap_cluster.pdf",width = 7, height = 6)
print(DimPlot(seurat_merge, reduction = "umap", pt.size = 0.2, label = TRUE))
dev.off()
#保存数据
save(seurat_merge, file = "seurat_merge_big_Type.RData")


#补充：修改内容
#根据Patients修改orig.ident里面的内容
seurat_merge@meta.data[seurat_merge@meta.data$Patients == "HBCP1" , "Patients"] = "HPCP4"  
seurat_merge@meta.data[seurat_merge@meta.data$Patients == "980755-1" , "Patients"] = "HPCP1" 
seurat_merge@meta.data[seurat_merge@meta.data$Patients == "980942-5" , "Patients"] = "HPCP2" 
seurat_merge@meta.data[seurat_merge@meta.data$Patients == "965789-3" , "Patients"] = "HPCP3" 
seurat_merge@meta.data[seurat_merge@meta.data$Patients == "10086106-4" , "Patients"] = "HPCP5"  
seurat_merge@meta.data[seurat_merge@meta.data$Patients == "10086106-6" , "Patients"] = "HPCP6" 
   #找的是列名而不是列里面的值
seurat_merge1@meta.data[seurat_merge1@meta.data$Patients == "1007877" , "orig.ident"] = "R231009005"
#更换细胞类型名字
seurat_merge@meta.data[seurat_merge@meta.data$cellTypes_new == "Neutrophils" , "cellTypes_new"] = "Undefined_Epi"
seurat_merge@meta.data[ seurat_merge$cellTypes_new=="Bcell", "cellTypes_new"] = "Immune"
seurat_merge@meta.data[ seurat_merge$cellTypes_new=="Tcell", "cellTypes_new"] = "Immune"
seurat_merge@meta.data[ seurat_merge$cellTypes_new=="Plasmocyte", "cellTypes_new"] = "Immune"
seurat_merge@meta.data[ seurat_merge$cellTypes_new=="Macrophage", "cellTypes_new"] = "Immune"
seurat_merge@meta.data[ seurat_merge$cellTypes_new=="Unsure_Immnue", "cellTypes_new"] = "Immune"

seurat_merge@meta.data[ seurat_merge$cellTypes_new=="Pericyte", "cellTypes_new"] = "Mesenchymal"
seurat_merge@meta.data[ seurat_merge$cellTypes_new=="Endothelial", "cellTypes_new"] = "Mesenchymal"
seurat_merge@meta.data[ seurat_merge$cellTypes_new=="Fibroblast", "cellTypes_new"] = "Mesenchymal"
seurat_merge@meta.data[ seurat_merge$cellTypes_new=="Muscle", "cellTypes_new"] = "Mesenchymal"

seurat_merge@meta.data[ seurat_merge$cellTypes_new=="Epithelial", "cellTypes_new"] = "Epithelial"

seurat_merge@meta.data[ seurat_merge$cellTypes_new=="Unsure", "cellTypes_new"] = "Unsure"
seurat_merge@meta.data[ seurat_merge$cellTypes_new=="Neutrophils", "cellTypes_new"] = "Undefined_Epi"





seurat_merge@meta.data[ rownames(seurat_merge4@meta.data)[seurat_merge4$cellTypes =="Tumor"], "cellTypes_new"] = "Tumor_Epithelial_cell"
seurat_merge@meta.data[ rownames(seurat_merge4@meta.data)[seurat_merge4$cellTypes =="Undefined"], "cellTypes_new"] = "Undefined"
#删数据框里面的对象
seurat_merge1 <- c[,seurat_merge$Patients != "HBCP~"]
#将Patients按照因子顺序排序
a <- c("HBCP1","HBCP2","HBCP3A","HBCP3B","HBCP3C","HBCP4A","HBCP4B","HBCP5A","HBCP5B","HBCP6","HBCP7","HBCP8","HBCP9","HBCP10","HBCP11","HBCP12","HBCP13","HBCP14","HBCP15","HBCP16","HBCP17","HBCP18","HBCP19","HBCP20","HBCP21","HBCP22")
seurat_merge1$Patients <- factor(seurat_merge1$Patients,levels = a )



#将Rdata文件转变成h5ad文件（Python）

library(MuDataSeurat)
load('/jdfsbjcas1/ST_BJ/P21H28400N0232/xieguixiang/BLCA/TSanalysis/SCdata/monocle/tumor_moduleScore.rds')
saveRDS(seurat_merge1,file='/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC_1/Data/seurat_merge_Immune.rds')

seurat <- readRDS('/jdfsbjcas1/ST_BJ/P21H28400N0232/xieguixiang/BLCA/TSanalysis/SCdata/monocle/tumor_moduleScore.rds')
WriteH5AD(seurat ,file = "/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/BC/Velocyto/SC_Tumor_Epi1.h5ad")
seurat_merge1<- DietSeurat(seurat_merge1) #将seurat进行瘦身
WriteH5AD(bm,file = "/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/BC/Velocyto/bm.h5ad")

i <- sapply(seurat_merge1@meta.data, is.factor)  #返回的结果取决于is.factor函数（判断检查对象是否为因子）
seurat_merge1@meta.data[i] <- lapply(seurat_merge1@meta.data[i], as.character) #和sapply差不多，但是lappy返回的是一个列表



#####寻找每两个cluster的Markergene
setwd('/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/BC/SingleCell/R28BC/Marker')
load('/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/BC/SingleCell/R28BC/seurat_merge_annot_subtype.RData')
c <- c('Bcell','Tcell','DCS','Endothelial','Fibroblast','Macrophage','Mastcell','Myocyte','Pericyte','Plasmocyte')
i <- 1
while(as.numeric(i) < 11){
    names <- c[i]
    # cluster.markers <- FindMarkers(seurat_merge, ident.1 = "Epithelial", ident.2 = names, min.pct = 0.25, logfc.threshold = 0.25 )
    load(paste0('/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/BC/SingleCell/R28BC/Marker','/',names,'_Marker.RData'))
    # cluster.markers %>%
    # # dplyr::filter(avg_log2FC > 1) %>%
    # dplyr::filter(pct.1 < 0.2 ) %>%
    # dplyr::filter(pct.2 > 0.6) %>%
    # arrange(desc(avg_log2FC)) %>%
    # slice_head(n = 200)  -> top200
    write.table(top200, file = paste0('/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/BC/SingleCell/R28BC/Marker','/',names,'_Marker.txt'),, sep = "\t", row.names = TRUE)
    # save(cluster.markers, file = paste0('/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/BC/SingleCell/R28BC/Marker','/',names,'_Marker.RData'))
    i <- i +1
}

