# usage: args1 h_matrix args2 w_matrix args3 seuratobj arg4 outdir arg5 prefix
library(dplyr)
library(Seurat)
# library(ppsr)
library(org.Hs.eg.db)
library(harmony)
args <- commandArgs(T)
h_matrix <- read.table(args[1], header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
w_matrix <- read.table(args[2], header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
# setwd("/jdfsbjcas1/ST_BJ/P21H28400N0232/gongchanghao/project/yixianai/SC_reanalysis/result/2_NMF_filter/Bcell")
# w_matrix <- read.table("/jdfsbjcas1/ST_BJ/P21H28400N0232/gongchanghao/project/yixianai/analysis/merge_split2/seurat_Bcell/Results/Figures/NMF/nmf_w_matrix.xls", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
# h_matrix <- read.table("/jdfsbjcas1/ST_BJ/P21H28400N0232/gongchanghao/project/yixianai/analysis/merge_split2/seurat_Bcell/Results/Figures/NMF/nmf_h_matrix.xls", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
h_matrix <- t(h_matrix)
# seurat_obj <- readRDS("/jdfsbjcas1/ST_BJ/P21H28400N0232/gongchanghao/project/yixianai/SC_reanalysis/result/1_cNMF/Bcell/seurat_obj.rds")
suffix <- strsplit(args[3], "\\.")[[1]][length(unlist(strsplit(args[3], "\\.")))]
if(suffix == "RData"){
    tmp <- load(args[3])
    seurat_obj <- get(tmp[1])
}else{
    seurat_obj <- readRDS(args[3])
}
# seurat_obj <- readRDS(args[3])
outdir <- args[4]
prefix <- args[5]

dir.create(paste0(outdir,"/2_NMF_filter/",prefix), recursive = TRUE)
dir.create(paste0(outdir,"/3_cluster/",prefix), recursive = TRUE)

setwd(paste0(outdir,"/2_NMF_filter/",prefix))
print(dim(h_matrix))
print(dim(w_matrix))
# 只保留Top 25% cell累计>0.5
# h_matrix <- apply(h_matrix, 2, function(x) {x/sum(x)})
filter_index_cell <- c()
for(i in 1:ncol(h_matrix)){
    tmp <- h_matrix[,i]
    sum_percent <- sum(tmp[order(tmp,decreasing = TRUE)][1:round(length(tmp)*0.25)])
    print(sum_percent)
    if(sum_percent < 0.5){
        filter_index_cell <- c(filter_index_cell,i)
    }
}
# 只保留Top 100 gene累计>0.5
# w_matrix <- apply(w_matrix, 1, function(x) {x/sum(x)})
# filter_index_gene <- c()
# for(i in 1:ncol(w_matrix)){
#     tmp <- w_matrix[,i]
#     sum_gene <- sum(tmp[order(tmp,decreasing = TRUE)][1:1000])
#     print(sum_gene)
#     if(sum_gene < 0.5){
#         filter_index_gene <- c(filter_index_gene,i)
#     }
# }

# 去除与批次相关性强的NMF
# batch_info <- seurat_obj$Patients
# filter_index_batch <- c()
# for(i in 1:ncol(h_matrix)){
#     tmp <- h_matrix[,i]
#     data <- data.frame(weight = tmp ,group = batch_info)
#     pps <- ppsr::score(data, x = 'group', y = 'weight', algorithm = 'tree')[['pps']]
#     print(pps)
#     if(pps > 0.5){
#         filter_index_batch <- c(filter_index_batch,i)
#     }
# }

# 去掉上述nmf
filter_nmf <- c(filter_index_cell)
# filter_nmf <- c(filter_index_cell,filter_index_gene,filter_index_batch)
filter_nmf <- filter_nmf[!duplicated(filter_nmf)]
print(filter_nmf)
if(!is.null(filter_nmf)){
 w_matrix <- w_matrix[,-filter_nmf]
}


# NMF W matrix Correlation>0.85的，合并取中位数



# 根据过滤的w_matrix，求解h_matrix
library(RcppML)
library(Matrix)
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 100)
data <- GetAssayData(seurat_obj, "data")
w <- as.matrix(w_matrix)
h_matrix_new <- RcppML::predict.nmf(w=w,data=data)
h_matrix_new <- t(apply(h_matrix_new, 1, function(x) {x/sum(x)}))
write.table(h_matrix_new,"h_matrix_new.txt",sep = '\t', row.names = TRUE, col.names = NA)
write.table(w_matrix,"w_matrix_new.txt",sep = '\t', row.names = TRUE, col.names = NA)
# 对top100gene进行enrichGO和enrich HALLMARK
library(clusterProfiler)
library(msigdbr)
m_df <- msigdbr(species = "Homo sapiens", category = "H")%>% 
      dplyr::select(gs_name, gene_symbol)
m_df$gs_name <- gsub("HALLMARK_","",m_df$gs_name)

gc_gsea <- list()
gc_hyper <- list()
NMF_list <- c()
for (i in c(1:nrow(h_matrix_new))) {
    #eval(parse(text = paste0("seurat_obj$NMF_", i, " <- as.numeric(model$h[", i, ",])")))
    # eval(parse(text = paste0("seurat_obj$NMF_",i," <- rep(0,ncol(seurat_obj))")))
    # eval(parse(text = paste0("seurat_obj$NMF_",i,"[seurat_obj$nFeature_RNA >100]<- as.numeric(h_matrix_new[",i,",])")))
    # pdf(file = paste0("../Results/Figures/NMF/nmf_activaty_umap_featureplots_", i, ".pdf"), width = 7, height = 7)
    # print(FeaturePlot(seurat_obj, features = paste0("NMF_", i), reduction = "umap"))
    # dev.off()
    NMF_list <- c(NMF_list, paste0("NMF_", i))
    tmp <- w_matrix[, i]
    names(tmp) <- rownames(w_matrix)
    tmp <- tmp[order(tmp, decreasing = TRUE)]
    tmp <- tmp[tmp > 0]
    tmp_n <- tmp[1:100]
    tmp_n <- list(tmp_n)
    names(tmp_n)[1] <- i
    tmp_n_id <- names(tmp_n[[1]])
    tmp_n_id <- list(tmp_n_id)
    names(tmp_n_id)[1] <- i
    gc_gsea <- append(gc_gsea, tmp_n)
    gc_hyper <- append(gc_hyper, tmp_n_id)

    #myenrich(tmp_n, tmp, "figures/NMF", paste0("NMF_", i))
}
print(gc_hyper)
ehm <- compareCluster(gc_hyper, fun = "enricher", TERM2GENE = m_df, pvalueCutoff = 0.05)
write.table(as.data.frame(ehm), "compare_hyper_HALLMARK.txt", sep = "\t", row.names = FALSE)
pdf("compare_hyper_HALLMARK_dotplot.pdf", width = 10, height = 10)
print(dotplot(ehm))
dev.off()
# ego <- compareCluster(gc_hyper, fun = "enrichGO", keyType = "SYMBOL", OrgDb = org.Hs.eg.db, pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, ont = "BP")
# write.table(as.data.frame(ego), "compare_hyper_GO_BP.txt", sep = "\t", row.names = FALSE)
# pdf("compare_hyper_GO_BP_dotplot.pdf", width = 10, height = 15)
# print(dotplot(ego))
# dev.off()

w_matrix_top100 <- data.frame()
for (i in c(1:ncol(w_matrix))) {
    tmp <- w_matrix[, i]
    names(tmp) <- rownames(w_matrix)
    tmp <- tmp[order(tmp, decreasing = TRUE)]
    tmp <- tmp[tmp > 0]
    tmp_n <- tmp[1:100]
    w_matrix_top100 <- rbind(w_matrix_top100, names(tmp_n))
}
rownames(w_matrix_top100) <- colnames(w_matrix)
write.table(t(w_matrix_top100), "nmf_top100_genes.xls", sep = "\t", row.names = FALSE)

# 根据H matrix，使用kNN进行clustering，得到cell subtype，并根据marker和DEG enrichment进行注释，对每个subtype命名

setwd(paste0(outdir,"/3_cluster/",prefix))
rownames(h_matrix_new) <- paste0("NMF_",1:nrow(h_matrix_new))
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all.genes)
# seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
seurat_obj[["nmf"]] <- CreateDimReducObject(embeddings = as.matrix(t(h_matrix_new)), loadings = as.matrix(w_matrix), key = "NMF_", assay = DefaultAssay(seurat_obj))
pdf("seurat_NMF_reduction.pdf")
DimPlot(seurat_obj, reduction = "nmf",group.by = "Patients", pt.size = 0.5)
dev.off()
seurat_obj <- seurat_obj %>% 
    RunHarmony("Patients", plot_convergence = FALSE,reduction = "nmf")
seurat_obj <- RunUMAP(seurat_obj, reduction = 'harmony', dims = 1:nrow(h_matrix_new))
# seurat_obj <- RunUMAP(seurat_obj, reduction = 'pca', dims = 1:30)
pdf(file = "umap_pca_batch.pdf", width = 10, height = 10)
print(DimPlot(seurat_obj, reduction = "umap", pt.size = 0.2, group.by = "Patients", label = TRUE))
print(DimPlot(seurat_obj, reduction = "umap", pt.size = 0.2, group.by = "group", label = TRUE))
dev.off()
seurat_obj <- FindNeighbors(seurat_obj, reduction = 'harmony', dims = 1:nrow(h_matrix_new))
# seurat_obj <- FindNeighbors(seurat_obj, reduction = 'pca', dims = 1:30)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.4)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.8)
seurat_obj <- FindClusters(seurat_obj, resolution = 1.2)
seurat_obj <- FindClusters(seurat_obj, resolution = 1)
pdf(file = "umap_cluster_nmf.pdf", width = 7, height = 6)
print(DimPlot(seurat_obj, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.4", label = TRUE))
print(DimPlot(seurat_obj, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.5", label = TRUE))
print(DimPlot(seurat_obj, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.8", label = TRUE))
print(DimPlot(seurat_obj, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.1", label = TRUE))
print(DimPlot(seurat_obj, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.1.2", label = TRUE))
dev.off()
height <- 2.5*(length(rownames(h_matrix_new))/4)
pdf(file = "umap_nmfs.pdf", width = 12, height = height)
print(FeaturePlot(seurat_obj, features =rownames(h_matrix_new),reduction = "umap", ncol = 4, raster = TRUE))
dev.off()
saveRDS(seurat_obj,"seurat_obj_nmf.rds")





