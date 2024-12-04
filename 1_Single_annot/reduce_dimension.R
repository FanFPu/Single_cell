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
scelist = list()
i <- "RNA"
setwd('/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC/Sample')
for(i in 1:length(dirs)){
  load(data.dir =dirs[[i]])
  scelist[[i]] <- CreateSeuratObject(counts = seurat_comb_singlet@assays$RNA@counts, project = paste0("R",i))
  scelist[[i]][["percent.mt"]] <- PercentageFeatureSet(scelist[[i]], pattern = "^MT-")
  scelist[[i]] <- subset(scelist[[i]], subset = percent.mt < 10)
}
load('/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC/Sample/HPCP6/seurat_HPCP6_single.RData')
i <- "HPCP6"

names(scelist)  = paste0("R",1:4)
sum(sapply(scelist, function(x)ncol(x@assays$RNA@counts)))
# scelist[[i]] <- CreateSeuratObject(counts = seurat_merge@assays$RNA@counts)
# scelist[[i]][["percent.mt"]] <- PercentageFeatureSet(scelist[[i]], pattern = "^MT-")
# scelist[[i]] <- subset(scelist[[i]], subset = percent.mt < 10)

# seurat_merge <- FindVariableFeatures(seurat_merge, selection.method = "vst", nfeatures = 1000)
# normalize and identify variable features for each dataset independently
scelist <- lapply(X = scelist, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
  }) 
features <- SelectIntegrationFeatures(object.list = scelist)
scelist1 <- FindIntegrationAnchors(object.list = scelist, anchor.features = features)
seurat_merge <- IntegrateData(anchorset = scelist1)
save(seurat_merge,file = '/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC/Test/sctrans_CCA/CCA/seurat_merge_CCA.RData')
DefaultAssay(seurat_merge) <- "integrated"
#####============================sub_hvg_pc==========================================
#针对于RNA
geneset_test <- c()
for (i in c('HPCP1','HPCP2','HPCP3','HPCP4','HPCP5','HPCP6')){
    sc_test <- subset(seurat_merge, Patients %in% c(i))
    dim(sc_test)
    sc_test <- NormalizeData(sc_test)
    sc_test <- FindVariableFeatures(sc_test, selection.method = "vst",nfeatures = 2000)
    geneset_test <- c(geneset_test, VariableFeatures(sc_test))
}
print(length(geneset_test))
# geneset_test <- na.omit(geneset_test)
geneset <- c()
for (i in geneset_test){
    if(table(geneset_test)[i] >= 4){
    geneset <- c(geneset, i)
    }
}
table(geneset_test)['CLIC2']
geneset <- unique(geneset)
print(length(geneset))
#针对于SCT
geneset_test <- c()
for (i in c('HPCP1','HPCP2','HPCP3','HPCP4','HPCP5','HPCP6')){
    sc_test <- subset(seurat_merge, Patients %in% c(i))
    sc_test@assays$SCT@scale.data <- sc_test@assays$RNA@scale.data
    dim(sc_test)
    sc_test <- NormalizeData(sc_test)
    sc_test <- PercentageFeatureSet(sc_test, pattern = "^MT-", col.name = "percent.mt")
    sc_test <- SCTransform(sc_test, vars.to.regress = "percent.mt", verbose = FALSE)
    geneset_test <- c(geneset_test, sc_test@assays$SCT@var.features)
}
print(length(geneset_test))
# geneset_test <- na.omit(geneset_test)
geneset <- c()
for (i in geneset_test){
    if(table(geneset_test)[i] >= 4){
    geneset <- c(geneset, i)
    }
}
table(geneset_test)['CLIC2']
geneset <- unique(geneset)
print(length(geneset))