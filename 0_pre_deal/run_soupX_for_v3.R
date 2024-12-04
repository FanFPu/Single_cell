library(SoupX)
library(Seurat)
args=commandArgs(T)

toc <- Read10X(args[1],gene.column=1) #filter expression
tod <- Read10X(args[2],gene.column=1) #backgroundd "02.cDNAAnno/RawMatrix/"
setwd(args[3]) #outdir
sample <- args[4] #id
#merge barcode
out_barcode <- read.table(args[5],header = FALSE) #03.M280UMI_stat/*_barcodeTranslate_16.txt
out_barcode <- out_barcode$V1

old <- CreateSeuratObject(toc)
tod <- CreateSeuratObject(tod)
#remove merge barcode
tod
tod <- subset(tod,cells=out_barcode,invert = TRUE)
tod
#merge cell and background
tod <- merge(tod,y=old)
tod
tod <- tod@assays$RNA@counts
tod <- tod[rownames(old),]
toc <- old@assays$RNA@counts
##get clusters
##old <- CreateSeuratObject(toc)
old <- NormalizeData(old,verbose = FALSE)
old <- FindVariableFeatures(old, selection.method = "vst", nfeatures = 3000)
old <- ScaleData(old)
old <- RunPCA(old,verbose = F)
old <- FindNeighbors(old, dims = 1:20)
old <- FindClusters(old, resolution = 0.5)
old <- RunUMAP(old, dims = 1:20)
#my_features <- c("Hba-a1","Hba-a2", "Hbb-bs", "Hbb-bt","Mki67","Top2a","Cd4","Cd8a","Klrd1","Cd19","Pecam1",'S100a8')
#pdf(paste0(args[4],"_marker_raw.pdf"),4.5*length(which(my_features %in% rownames(old))),4)
#p <- FeaturePlot(old, features = my_features, ncol=length(which(my_features %in% rownames(old))), order = T, reduction = "umap")
#print(p)
#dev.off()

pdf(paste0(args[4],"_Cluster.pdf"),8,7)
p <- DimPlot(object = old, reduction = "umap",pt.size = 0.1,label=T)
print(p)
dev.off()
#run soupX
matx <- old@meta.data
sc = SoupChannel(tod, toc)
sc = setClusters(sc, setNames(matx$seurat_clusters, rownames(matx)))
sc = autoEstCont(sc,forceAccept =T,doPlot =F)
rho <- unique(sc$metaData$rho)

#auto
out = adjustCounts(sc)
auto <- CreateSeuratObject(out, min.cells = 1, min.features=1)
saveRDS(auto,paste0(args[4],"_soupX_auto.rds"))
auto <- NormalizeData(auto,verbose = FALSE)
auto <- FindVariableFeatures(auto, selection.method = "vst", nfeatures = 3000)
auto <- ScaleData(auto)
auto@reductions <- old@reductions

#pdf(paste0(args[4],"_marker_rm_auto.pdf"),4.5*length(which(my_features %in% rownames(auto))),4)
#p <- FeaturePlot(auto, features = my_features, ncol=length(which(my_features %in% rownames(auto))), order = T, reduction = "umap")
#print(p)
#dev.off()

#add 10
rho <- rho + 0.1
sc = setContaminationFraction(sc, rho,forceAccept=T)
print(unique(sc$metaData$rho))
out = adjustCounts(sc)
add10 <- CreateSeuratObject(out, min.cells = 1, min.features=1)
saveRDS(add10,paste0(args[4],"_soupX_add10.rds"))
add10 <- NormalizeData(add10,verbose = FALSE)
add10 <- FindVariableFeatures(add10, selection.method = "vst", nfeatures = 3000)
add10 <- ScaleData(add10)
add10@reductions <- old@reductions


#pdf(paste0(args[4],"_marker_add10.pdf"),4.5*length(which(my_features %in% rownames(add10))),4)
#p <- FeaturePlot(add10, features = my_features, ncol=length(which(my_features %in% rownames(add10))), order = T, reduction = "umap")
#print(p)
#dev.off()

	
