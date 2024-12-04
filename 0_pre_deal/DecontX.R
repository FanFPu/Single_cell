#==============================decontX=============
#R版本需要4.3以上
library(Seurat)
library(decontX)
library(ggplot2)

# 读取seurat对象
seurat <- readRDS('test.rds')
# 提取矩阵(genes x cells)
# 需要注意的是，你应该查看你的矩阵是否已经注释了行名和列明，以及如果你的矩阵是scanpy对象中提取的矩阵应该行列转置
counts <- seurat@assays$RNA@counts
# 计算
decontX_results <- decontX(counts) 
str(decontX_results)
# RNA污染的计算结果储存在：decontX_results$contamination
# 他是一个数值型向量，长度与细胞数量一致，顺序也与矩阵的colnames一致
# 我们可以直接把他写入metadata
seurat$Contamination =decontX_results$contamination
head(seurat@meta.data) 

FeaturePlot(seurat, 
            features = 'Contamination', 
            raster=FALSE       # 细胞过多时候需要加这个参数
           ) + 
  scale_color_viridis_c()+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  xlab('scVI_UMAP_1')+
  ylab('scVI_UMAP_2')
save(seurat,file= 'decontX.RData')
# 在这里看到，对比前文最开始的状态，我们经过了严格的质控之后，散在的细胞已经少了很多，但是仍然有亚群之间不清晰的界限
# 我们可视化Contamination，惊奇地发现这些边缘的毛躁就是Contamination(污染)较高地细胞
