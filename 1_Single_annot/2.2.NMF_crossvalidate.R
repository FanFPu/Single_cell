##加载包
library(RcppML)
library(ggplot2)
library(SingleCellExperiment)
library(Matrix)
library(Seurat)
options(RcppML.threads = 3) #指定RcppML库的线程数为3，此库是提高R代码的执行效率
library(MuDataSeurat)

# args <- commandArgs(trailingOnly = TRUE) # get the arguments   #获取 R 脚本的命令行参数,返回的是一个字符串向量
# load(args[1])
# od <- dirname(args[1])
# setwd(od)
# dir.create("figures/NMF", recursive = TRUE)
workdir <- '/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC/NMF'
if(!dir.exists(workdir)){dir.create(workdir, recursive = TRUE)}
setwd(workdir)


merged <- ReadH5AD("/jdfsbjcas1/ST_BJ/P21H28400N0232/xieguixiang/BLCA/TSanalysis/Total_tumor_bigCell_filterNor.h5ad")
# merged <- readRDS('/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/BC/NMF/Data/Total_tumor_bigCell.rds')
# merged <- NormalizeData(merged)
load('/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC/Data/seurat_notharmony/1.1.4_seurat_merge_normalize.RData')

seurat_spatialObj <- seurat_merge
# seurat_spatialObj <- subset(seurat_spatialObj, subset = nFeature_Spatial > 100)
seurat_spatialObj <- subset(seurat_spatialObj, subset = nFeature_RNA > 100)
table(seurat_spatialObj@meta.data$merged_cluster)
# sce <- as.SingleCellExperiment(seurat_spatialObj)
# seurat_spatialObj <- NormalizeData(seurat_spatialObj, assay = "Spatial", verbose = FALSE)
# data <- assay(sce, "logcounts")

####################数据处理
# seurat_spatialObj@assays$RNA@data <- seurat_spatialObj@assays$RNA@data[, colSums(!is.finite(seurat_spatialObj@assays$RNA@data)) == 0] #删掉含有nan的列
# str(seurat_spatialObj)
# seurat_spatialObj@assays$RNA@data@x <- seurat_spatialObj@assays$RNA@data@x[,!is.nan(seurat_spatialObj@assays$RNA@data@x)]  #删掉含有nan的
# seurat_spatialObj@assays$RNA@data@x <- na.omit(seurat_spatialObj@assays$RNA@data@x) #将其中含有NaN的行删掉

data <- seurat_spatialObj@assays$RNA@data
data <- GetAssayData(seurat_spatialObj@assays$RNA@data, slot = 'scale.data')  #获取特定的数据，此处是data的测定数据
#用于执行交叉验证，会将数据集分成若干份（通常是 k 折,由参数 k 控制），在每一份上训练模型，并在剩余的数据上进行评估
cv_mse <- crossValidate(data, method = "predict", k = 1:40, seed = 123)
cv_robust <- crossValidate(data, method = "robust", k = 1:40, seed = 123) #鲁棒性强的方法进行模型的训练和预测
pdf("NMF_CV.pdf", width = 10, height = 7)
P1 <- plot(cv_mse)
P2 <- plot(cv_robust)
print(P1 + P2)
dev.off()
