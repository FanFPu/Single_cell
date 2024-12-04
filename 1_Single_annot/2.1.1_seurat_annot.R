options(stringsAsFactors = FALSE)
library(Seurat)
library(dplyr)
library(DoubletFinder)
library(limma)
library(Cairo)
library(harmony)
library(car)
library(SingleCellExperiment)
library(ggplot2)
library(data.table)
options(bitmapType = "cairo")

#多线程分析
plan()
plan("multisession", workers = as.numeric(2))
# Enable parallelization
plan()
options(future.globals.maxSize= 4294967296000)

#
source("/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/code/Single_cell/codes/scRNA_primary.R")
source("/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/code/Single_cell/codes/colors.R")
source("/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/code/Single_cell/codes/signatureGenes.R") #human signature gene 
#source("/jdfsbjcas1/ST_BJ/P21H28400N0232/wenfeng/project/codes/signatureGene_mouse.R")
source("/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/code/Single_cell/codes/0_initNEWanalysis.R")
source("/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/code/Single_cell/codes/DEG_wilcox.R")


###注释

load('/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC/Data/seurat_harmony/1.1.4_seurat_merge_harmony.RData')

##大类 "Epithelial","Immune","Mesenchymal",

#下面两行是同等作用，第一行是简写
new.cluster.ids <- rep("Epithelial", times = 27)
new.cluster.ids <- c("Epithelial","Epithelial","Epithelial","Epithelial","Epithelial","Epithelial","Epithelial","Epithelial","Epithelial","Epithelial","Epithelial","Epithelial","Epithelial","Epithelial","Epithelial","Epithelial","Epithelial","Epithelial","Epithelial","Epithelial")
new.cluster.ids <-c("Epithelial","Mesenchymal","Epithelial","Epithelial","Epithelial","Immune","Epithelial","Immune","Mesenchymal","Mesenchymal","Mesenchymal","Immune","Epithelial","Epithelial","Epithelial","Mesenchymal","Epithelial")
length(new.cluster.ids)
new.cluster.ids <- c("Unsure","Epithelial","Epithelial","Mesenchymal","Mesenchymal","Epithelial","Mesenchymal", "Epithelial", "Epithelial","Mesenchymal","Epithelial","Immune","Immune","Immune", "Epithelial", "Immune", "Immune", "Epithelial", "Unsure", "Epithelial", "Unsure", "Mesenchymal", "Immune", "Unsure",  "Immune", "Epithelial",  "Epithelial", "Epithelial",  "Immune")

new.cluster.ids <- c("Undefined","Mesenchymal","Epithelial","Epithelial","Mesenchymal","Immune","Epithelial","Epithelial","Epithelial","Mesenchymal","Epithelial","Immune","Mesenchymal","Epithelial","Immune","Mesenchymal","Epithelial","Immune","Undefined","Immune","Epithelial","Mesenchymal")
new.cluster.ids <- c("Fibroblast","Fibroblast","Fibroblast","Fibroblast","Fibroblast","Pericyte","Fibroblast", "Muscle", "Endothelial","Endothelial","Fibroblast","Endothelial","Endothelial","Endothelial","Endothelial","Endothelial","Endothelial","Fibroblast","Fibroblast","Pericyte")

new.cluster.ids <- c("Nervous","Epithelial","Epithelial","Epithelial","Mesenchymal","Immune","Epithelial","Immune","Immune","Mesenchymal","Mesenchymal","Mesenchymal","Epithelial","Immune","Epithelial","Mesenchymal","Immune","Mesenchymal","Immune","Mesenchymal","Mesenchymal","Mesenchymal","Mesenchymal","Epithelial","Melanin")



#seurat_merge <- FindNeighbors(seurat_merge, dims = 1:30)
seurat_merge <- FindClusters(seurat_merge, resolution =0.4,algorithm = 1)
names(new.cluster.ids) <- levels(seurat_merge)
seurat_merge <- RenameIdents(seurat_merge, new.cluster.ids)
seurat_merge$cellTypes_new <- Idents(seurat_merge) 

# Idents(seurat_merge4) = seurat_merge$RNA_snn_res.0.4
# anno_mannual <- c("Tcell","Plasmocyte","Bcell","Plasmocyte","Macrophage","Neutrophils","Neutrophils","Unsure_Immnue","Tcell","Macrophage","Plasmocyte","Macrophage","Neutrophils","Unsure_Immnue","Bcell","Neutrophils","Tcell","Neutrophils","Unsure_Immnue","Neutrophils","Bcell","Tcell","Neutrophils","Unsure_Immnue","Unsure_Immnue")
# names(anno_mannual) <- seq(0,24)
# seurat_merge <- RenameIdents(seurat_merge, anno_mannual)#进行重新命名
pdf(file = "/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC/Test/sctrans_CCA/sctrans/umap_bigType1.pdf",width = 7, height = 6)
print(DimPlot(seurat_merge, reduction = "umap",  pt.size = 0.2, raster = FALSE,label = TRUE))
print(DimPlot(seurat_merge, reduction = "umap", pt.size = 0.2, raster = FALSE,group.by = "cellTypes_new", label = T,label.size = 4,cols=c('#a2cffe', '#87a922', '#ffa62b', '#f8481c', '#cffdbc'))) #cols=c(Bcell='red')'#a6814c', '#a484ac', '#fc86aa', '#952e8f', '#02ccfe','#2000b1', '#009337', '#ad0afd'
dev.off()



seurat_merge@meta.data[ rownames(seurat_merge1@meta.data)[seurat_merge1$cellTypes_new =="Macrophage"] , "cellTypes_new"] = "Macrophage"
seurat_merge@meta.data[ rownames(seurat_merge1@meta.data)[seurat_merge1$cellTypes_new =="Bcell"], "cellTypes_new"] = "Bcell"
seurat_merge@meta.data[ rownames(seurat_merge1@meta.data)[seurat_merge1$cellTypes_new =="Tcell"], "cellTypes_new"] = "Tcell"
seurat_merge@meta.data[ rownames(seurat_merge1@meta.data)[seurat_merge1$cellTypes_new =="Plasmocyte"], "cellTypes_new"] = "Plasmocyte"
# seurat_merge@meta.data[ rownames(seurat_merge1@meta.data)[seurat_merge1$cellTypes_new =="Neutrophils"], "cellTypes_new"] = "Neutrophils"
seurat_merge@meta.data[ rownames(seurat_merge1@meta.data)[seurat_merge1$cellTypes_new =="DCS"], "cellTypes_new"] = "DCS"
seurat_merge@meta.data[ rownames(seurat_merge1@meta.data)[seurat_merge1$cellTypes_new =="Mastcell"], "cellTypes_new"] = "Mastcell"

seurat_merge@meta.data[ rownames(seurat_merge2@meta.data)[seurat_merge2$cellTypes_new =="Undefined"], "cellTypes_new"] = "Undefined"

seurat_merge@meta.data[ rownames(seurat_merge4@meta.data)[seurat_merge4$cellTypes =="Normal"], "cellTypes_new"] = "Normal_Epithelial_cell"
seurat_merge@meta.data[ rownames(seurat_merge4@meta.data)[seurat_merge4$cellTypes =="Tumor"], "cellTypes_new"] = "Tumor_Epithelial_cell"
seurat_merge@meta.data[ rownames(seurat_merge4@meta.data)[seurat_merge4$cellTypes =="Undefined"], "cellTypes_new"] = "Undefined"

seurat_merge@meta.data[  rownames(seurat_merge4@meta.data), "cellTypes_new"] = "Unsure"
 

pdf(file = "/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC/Results/1.1.4_seurat_harmony/umap_subType1.pdf",width = 7, height = 6)
print(DimPlot(seurat_merge, reduction = "umap", group.by = "cellTypes_new", pt.size = 0.2, label = TRUE,raster=FALSE))

print(DimPlot(seurat_merge, reduction = "umap", pt.size = 0.1, group.by = "cellTypes_new", label = T,label.size = 2,cols=c("#02ccfe", "#ceda3f", "#7e39c9", "#72d852", "#d849cc","#5e8f37", "#5956c8", "#dd3e34", "#d293c8", "#DABE8B","#8d378c", "#68d9a3", "#2A2D5E36"))) #cols=c(Bcell='red')dd3e34
dev.off()
#"#b377cf", "#899140", "#564d8b", "#ddb67f", "#292344", "#d0cdb8", "#421b28", "#5eae99", "#a03259", "#2000b1", "#e598d7", "#406024", "#bbb5d9"

#DimPlot(ifnb, reduction = "umap", split.by = "stim")


# ===================================3D___UMAP===================
library(plotly)
library(Seurat)
# Re-run UMAPs that you have accurate calculations for all UMAP(s)
seurat_merge2 <- RunUMAP(seurat_merge,
                            dims = 1:10,
                            n.components = 4L)

umap_1 <- seurat_merge1[["umap"]]@cell.embeddings[,1]
umap_2 <- seurat_merge1[["umap"]]@cell.embeddings[,2]
umap_3 <- seurat_merge1[["umap"]]@cell.embeddings[,3]
umap_4 <- seurat_merge1[["umap"]]@cell.embeddings[,4]
pdf(file = "FeaturePlot2_KRT17.pdf", width = 9, height = 8)
print(FeaturePlot(seurat_merge, features = c("KRT5", "KRT14", "KRT17"),raster = TRUE,dims= c(2,3)))
dev.off()
# Prepare a dataframe for cell plotting
plot.data <- FetchData(object = seurat_merge1, vars = c("UMAP_1", "UMAP_2", "UMAP_3", "Patients"))

# Make a column of row name identities (these will be your cell/barcode names)
plot.data$label <- paste(rownames(plot.data))

plot.data <- FetchData(object = seurat_merge1, vars = c("UMAP_1", "UMAP_2", "UMAP_3", "KRT17"), slot = 'data')

# Say you want change the scale, so that every cell having an expression >1 will be one color
# Basically, you are re-adjusting the scale here, so that any cell having a certain expression will light up on your 3D plot
# First make another column in your dataframe, where all values above 1 are re-assigned a value of 1
# This information is stored in the 'changed' column of your dataframe
plot.data$changed <- ifelse(test = plot.data$BID <1, yes = plot.data$KRT17, no = 1)
# Add the label column, so that now the column has 'cellname-its expression value'
plot.data$label <- paste(rownames(plot.data)," - ", plot.data$KRT17, sep="")
# Plot your data, in this example my Seurat object had 21 clusters (0-20), and cells express a gene called ACTB
plot_ly(data = plot.data, 
        x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
        color = ~changed, # you can just run this against the column for the gene as well using ~ACTB, the algorith will automatically scale in that case based on maximal and minimal values
        opacity = .5,
        colors = c('darkgreen', 'red'), 
        type = "scatter3d", 
        mode = "markers",
        marker = list(size = 5, width=2), 
        text=~label,
        hoverinfo="text"
)

# create a dataframe
goi <- "KRT17"
plotting.data <- FetchData(object = seurat_merge1, vars = c("UMAP_1", "UMAP_2", "UMAP_3", "Expression"=goi), slot = 'data')
# Say you want change the scale, so that every cell having an expression >1 will be one color
# Basically, you are re-adjusting the scale here, so that any cell having a certain expression will light up on your 3D plot

# First make another column in your dataframe, where all values above 1 are re-assigned a value of 1
# This information is stored in the 'Expression' column of your dataframe
# Cutoff <- 2
Cutoff <- quantile(plotting.data[,goi], probs = .95)
plotting.data$"ExprCutoff" <- ifelse(test = plotting.data[,goi] <Cutoff, yes = plotting.data[,goi], no = Cutoff)

# Add the label column, so that now the column has 'cellname-its expression value'
plotting.data$label <- paste(rownames(plotting.data)," - ", plotting.data[,goi], sep="")

# Plot your data
pdf(file = "3d_KRT17.pdf", width = 9, height = 8)
print(plot_ly(data = plotting.data,
        # name = goi,
        x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
        color = ~ExprCutoff, # you can just run this against the column for the gene as well using ~ACTB, the algorith will automatically scale in that case based on maximal and minimal values
        opacity = .5,
        colors = c('darkgrey', 'red'), 
        type = "scatter3d", 
        mode = "markers",
        marker = list(size = 10), 
        text=~label,
        hoverinfo="text"
) %>% layout(title = goi))
dev.off()

save(seurat_merge1,file= '/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/BC/SingleCell/Epithelial/Epithelial_PCA_3D.RData' )
save(plotting.data,file= '/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/BC/SingleCell/Epithelial/plotting.data_3D.RData' )