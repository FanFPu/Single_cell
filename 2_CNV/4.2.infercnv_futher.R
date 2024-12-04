library('Seurat')
library('infercnv')
library('dendextend')
library('phylogram')
library('miscTools')
library('ggthemes')
library('RColorBrewer')
library('umap')
library('ggplot2')
library('car')
library('limma')
library('tibble')
library('dplyr')
library('stringr')
options(bitmapType = "cairo")

rm(list = ls())
gc()


# args = commandArgs(T)
sample = "HPCP3"  # args[1]

#工作路径
topdir = '/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC/SC_infercnv/SC_infercnv1.5'
indir = paste0(topdir, '/Results/1.4.1_run_inferCNV/ref_T_B_macrophage_endothelial/', sample)
workdir=paste0(topdir, '/Results/1.4.2_analyze_inferCNV/ref_T_B_macrophage_endothelial/', sample)
if(!dir.exists(indir)){dir.create(indir, recursive = TRUE)}
if(!dir.exists(workdir)){dir.create(workdir, recursive = TRUE)}
setwd(workdir)

#
k1 = as.numeric(20)
tree = read.dendrogram(paste0(indir, '/infercnv.observations_dendrogram.txt'))
labels_hclust2 = cutree(tree, k=k1)   
labels_hclust2 = data.frame(row.names = names(labels_hclust2), 'subclone'=paste0('Subclo_', labels_hclust2))  #行名细胞，列名subclone

#
load('/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/BC/Spatial_inferCNV/Results/1.4.1_run_inferCNV/ref_T_B_macrophage_endothelial/bigcell_merged.rds')
seurat_merge1 <- subset(seurat_merge,subset = merged_cluster %in% c("TS1","TS2","TS3","TS4","TS5","TS6"))
seurat_i <- subset(seurat_i,subset = cellTypes_new %in% c('Epithelial'))
labels_hclust2 = data.frame(row.names = rownames(seurat_merge1@meta.data), 'subclone' = seurat_merge1@meta.data$merged_cluster)

## 将cluster再分一下
# subdends = get_subdendrograms(as.dendrogram(tree), k=k1, order_clusters_as_data = T)[[16]]
# labels_sub = cutree(subdends, k=4, order_clusters_as_data = T)
# labels_sub2 = paste0("Subclo_16.", labels_sub)
# names(labels_sub2) = names(labels_sub)
# labels_hclust2[names(labels_sub2),] <- labels_sub2

load(paste0(topdir,'/Results/1.4.1_run_inferCNV/ref_T_B_macrophage_endothelial/', sample, "/infercnv_obj2.RData"))  #infercnv_obj
infercnv_obj2 <-  infercnv_obj
#将对照细胞提取出来
ref_ids = setdiff(colnames(infercnv_obj2@expr.data), rownames(labels_hclust2))   # setdiff(A, B) 返回在集合 A 中出现但不在集合 B 中出现的元素。换句话说，它返回属于 A 但不属于 B 的元素
#将其细胞名进行合并，不漏掉一个细胞
infercnv_obj2@expr.data = infercnv_obj2@expr.data[,union(rownames(labels_hclust2), ref_ids)]   #union()计算两个或多个向量的并集,且不包含重复的元素

name = list() #建立空列表
index = list()
index_length = list()
hc_list = list()


for (i in unique(labels_hclust2[,1])) {
    name[[i]] <- rownames(labels_hclust2)[labels_hclust2$subclone == i]  #将为subclone【K】的细胞名取出来
    index[[i]] <- which(colnames(infercnv_obj2@expr.data) %in% name[[i]])  #which建立满足括号里的索引
    index_length[[i]] <- seq(from=1,to=length(name[[i]]))   #建立等差数列，总共是subclone【K】的数量
    hc_list[[i]] = list(order=seq(from=1,to=length(name[[i]])),labels=name[[i]])  # 创建两个列表，一个是index_length的列表，一个是将每个亚型的细胞贴上标签，内容是细胞名，
}

infercnv_obj2@observation_grouped_cell_indices <- index_length
infercnv_obj2@tumor_subclusters$hc <- hc_list
infercnv_obj2@reference_grouped_cell_indices = list('ref' = which(colnames(infercnv_obj2@expr.data) %in% ref_ids))  #对照的细胞名，建立索引，然后单独拿出来,内容是索引号
length(table(infercnv_obj2@reference_grouped_cell_indices)) #对照细胞类型的细胞数量
saveRDS(infercnv_obj2, file=paste0(workdir, '/infercnv_obj2.rds'))

source('/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/My_code/Single_cell/SC_1.4.10_infercnv_heatmap_plot.R')  #里面有plot_cnv函数
plot_cnv(infercnv_obj2,
        out_dir=workdir,
        title="inferCNV",
        obs_title="Observations (Cells)",
        ref_title="References (Cells)",
        cluster_by_groups=T,
        cluster_references=FALSE,
        plot_chr_scale=FALSE,
        chr_lengths=NULL,
        k_obs_groups = 1,
        contig_cex=1,
        x.center=mean(infercnv_obj2@expr.data),
        x.range="auto", #NA,
        hclust_method='ward.D',
        custom_color_pal=NULL,
        color_safe_pal=FALSE,
        output_filename="infercnv",
        output_format="pdf" , #pdf, png, NA
        png_res=300,
        dynamic_resize=0,
        ref_contig = NULL,
        write_expr_matrix=F,
        useRaster=TRUE)
                   

#*根据cnv树状图区分亚型
#load(paste0(workdir, '/infercnv_obj2.RData'))
#sample = "HBCP10"
tree = read.dendrogram(paste0(indir, '/infercnv.observations_dendrogram.txt'))  #从文件或其他数据源中读取聚类树，其实并不直接读文件，而是对特定R对象创建dendrogram（聚类数）对象
labels = cutree(tree, k=k1, order_clusters_as_data = F)  #k具有所需组数的整数标量或向量


## 将cluster再分一下
# subdends = get_subdendrograms(as.dendrogram(tree), k=k1, order_clusters_as_data = T)[[16]]
# labels_sub = cutree(subdends, k=4, order_clusters_as_data = T)
# labels_sub2 = paste0("16.", labels_sub)
# names(labels_sub2) = names(labels_sub)
# labels[names(labels_sub2)] <- labels_sub2
# print(table(labels))

# #将labels重新编码
cell_names <- names(labels) #提取细胞名
labels = Recode(labels, "c(1,2,3,4,5,6,7,10,11,12,13,14,15,16,17,18,19,20)= 1;c(8) =2;c(9)=3")  #基于dplyr包，用于修改或重新编码变量的取值
names(labels) = cell_names
print(table(labels))

bars = as.data.frame(c(brewer.pal(n = 10, name = "Set1"), brewer.pal(n = 11, name = "Set2"))[-11][labels])  #生成颜色调色板，特别是在数据可视化中使用。 #其生成的数据框列名就是后面的一大串
colnames(bars) = 'inferCNV' #将其列名从上面的c（....）改成infercnv
bars$inferCNV = as.character(bars$inferCNV)  #将其变成字符串的形式

#生成重新编码clusters后的调色板
tiff(file = paste0(workdir, '/', sample, 'inferCNV.dend.tiff'), height = 5, width = 14, res = 300, compression="lzw", units="in")
print(tree %>% dendextend::set('labels', rep('', nobs(tree))) %>%  # nobs() 函数通常用于获取对象的观测数
  plot(main='inferCNV') %>%  #main 是一个用于设置主标题的参数，可以是一个字符向量，可以理解为图像的标题
  colored_bars(colors = as.data.frame(bars), dend=tree, add = T, y_scale = 100, sort_by_labels_order = F, y_shift = 0, horiz = F))
dev.off()
#horiz:条形的方向T水平，F垂直；dend是树状图对象；add是添加元素的函数，T就是可以添加，sort是对其标签进行排序，y_scale:条形图应在 y 轴上拉伸多少？如果没有提供 dend - 默认值为 1


#*画umap图--根据cnv matrix进行umap降维，贴上亚型标签

# 提取单个样本的数据
seurat_i_path = paste0(workdir, '/1.4.2_seurat_', sample, '.RData')
if (file.exists(seurat_i_path)) {
    load(seurat_i_path)
}else {
    seurat_path = paste0(topdir, "/saveData/1.4.1_seurat_merge_infercnv.RData")  #全部样本的注释信息，按照Patients进行分开，细胞类型只有免疫、内皮和上皮
    load(seurat_path)
    seurat_i = seurat_lists[[sample]]     #seurat_lists是每个病人放在一个list中，其中包含了所有的注释信息
    save(seurat_i, file=seurat_i_path)
}
# load(seurat_i_path)
# 单个样本的cnv打分数据提取出来，打分数据存在infercnv.observations.txt 文件中
seurat_Epithe = subset(seurat_i, subset = cellTypes_new %in% 'Epithelial')  #可以使用%in%, 将一个样本的上皮提取出来
if (file.exists(paste0(workdir, '/1.4.2_cnv_umap_', sample, '.RData'))) {
    load(paste0(workdir, '/1.4.2_cnv_umap_', sample, '.RData'))
} else {
    cnv_mtx = read.delim(paste0(indir, '/infercnv.observations.txt'), sep=' ')  #读完之后列名是细胞，行名是gene   该文件是对gene进行打分#取以制表符分隔的文本文件，一般读取数据框，以“\t”制表符分隔
	cnv_mtx2 = as(t(cnv_mtx), "dgCMatrix")   #as() 是一个通用的转换函数，dense general," 即一般的密集矩阵；dgc稀疏矩阵，压缩存储零元素，可以显著减少内存占用
	cnv_umap = umap(t(cnv_mtx)) #umap非线性降维技术；t() 函数是用于转置矩阵的函数。它可以用于交换矩阵的行和列，将矩阵的行列互换。
	save(cnv_umap, file = paste0(workdir, '/1.4.2_cnv_umap_', sample, '.RData'))
}

# load(paste0(workdir, '/1.4.2_cnv_umap_', sample, '.RData'))
#画gene打分的可视化图，umap_cnv
#cnv_umap$layout = t(cnv_umap$layout)
all(rownames(cnv_umap$layout) == names(labels))   #all() 函数是用于检查给定条件是否对向量中的所有元素都为真的函数
cnv_dim = data.frame(as.data.frame(cnv_umap$layout), cluster = factor(labels, levels=1:length(unique(labels)), ordered=T)) #新增一列，将其聚类类别K个设为因子改成一列
colnames(cnv_dim) = c("UMAP_1", "UMAP_2", "cellSubtype") #重新定义列名,UMAP_1和UMAP_2是对Gene的打分
#画两种图。表达内容一样 pdf和tif
pdf(file = paste0(workdir, '/1.4.2_', sample, '_umap_cnv.pdf'),width = 9, height = 8)
print(SingleDimPlot(cnv_dim, dims = c("UMAP_1","UMAP_2"), pt.size = 1, col.by = "cellSubtype", cols = c(brewer.pal(n = 9, name = "Set1"), brewer.pal(n = 8, name = "Set2"))[-11][unique(labels)]) + theme(plot.title = element_text(hjust = 0.5)))
dev.off()
tiff(file = paste0(workdir, '/1.4.2_', sample, '_umap_cnv.tiff'), res=300, width=9, height=8, compression="lzw", units="in")
print(SingleDimPlot(cnv_dim, dims = c("UMAP_1","UMAP_2"), pt.size = 1, col.by = "cellSubtype", cols = c(brewer.pal(n = 9, name = "Set1"), brewer.pal(n = 8, name = "Set2"))[-11][unique(labels)]) + theme(plot.title = element_text(hjust = 0.5)))
dev.off()

#*保存亚型标签信息
seurat_Epithe$cellSubtype = labels[colnames(seurat_Epithe)]  #seurat_Epithe的列名是细胞，行名是Gene; labels是cnv的分类
save(seurat_Epithe, file = paste0(workdir, '/1.4.2_seurat_Epithe_', sample, '.RData'))  #此处保存的是样本的仅上皮信息

##读取infercnv的cluster信息
#load(paste0(workdir, '/1.4.2_seurat_Epithe_', sample, '.RData'))
# seurat_Epithe$cellSubtype_infercnv = NULL   
# seurat_Epithe$cellSubtype_infercnv_ifNormal=NULL

#将seurat对象的列名修改，并用因子进行排序
names(seurat_Epithe@meta.data) = replace(names(seurat_Epithe@meta.data), names(seurat_Epithe@meta.data)=='cellSubtype', 'cellSubtype_infercnv')
seurat_Epithe$cellSubtype_infercnv = factor(paste0('CNV-', seurat_Epithe$cellSubtype_infercnv), levels = paste0('CNV-', sort(unique(seurat_Epithe$cellSubtype_infercnv))), ordered = T)
print(table(seurat_Epithe$cellSubtype_infercnv))

save(seurat_Epithe,file= paste0(workdir, "/1.4.2_seurat_Epithe_scaled_", sample, ".RData"))
##根据expr对epithe区分cluster
seurat_Epithe <- NormalizeData(seurat_Epithe)
seurat_Epithe <- FindVariableFeatures(seurat_Epithe, selection.method = "vst", nfeatures = 2000)
seurat_Epithe <- ScaleData(seurat_Epithe)
seurat_Epithe <- RunPCA(seurat_Epithe)
seurat_Epithe <- RunUMAP(seurat_Epithe, dims = 1:30)
# seurat_Epithe <- RunTSNE(seurat_Epithe, dims = 1:30)
# save(seurat_Epithe, file = paste0(workdir, "/1.4.2_seurat_Epithe_scaled_", sample, ".RData"))

seurat_Epithe <- FindNeighbors(seurat_Epithe, dims = 1:30)
seurat_Epithe <- FindClusters(seurat_Epithe, resolution = seq(0.1,1.2,0.1))  
pdf(file = paste0(workdir, "/1.4.2_umap_cluster.pdf"),width = 7, height = 6)
print(DimPlot(seurat_Epithe, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.1", label = TRUE))
print(DimPlot(seurat_Epithe, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.2", label = TRUE))
print(DimPlot(seurat_Epithe, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.3", label = TRUE))
print(DimPlot(seurat_Epithe, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.4", label = TRUE))
print(DimPlot(seurat_Epithe, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.5", label = TRUE))
print(DimPlot(seurat_Epithe, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.6", label = TRUE))
print(DimPlot(seurat_Epithe, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.7", label = TRUE))
print(DimPlot(seurat_Epithe, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.8", label = TRUE))
print(DimPlot(seurat_Epithe, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.0.9", label = TRUE))
print(DimPlot(seurat_Epithe, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.1", label = TRUE))
print(DimPlot(seurat_Epithe, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.1.1", label = TRUE))
print(DimPlot(seurat_Epithe, reduction = "umap", pt.size = 0.2, group.by = "RNA_snn_res.1.2", label = TRUE))
dev.off()
save(seurat_Epithe, file = paste0(workdir, "/1.4.2_seurat_Epithe_scaled_", sample, ".RData"))

#因为res.0.7将所有离散的群都标注上了label，所以res=0.7。根据res=0.7时的umap图，将cluster 0-6,8,13分到一类。
# load(paste0(workdir, "/1.4.2_seurat_Epithe_scaled_", sample, ".RData"))
res = 0.4  # args[2]
seurat_Epithe$'cellSubtype_expr' = seurat_Epithe@meta.data[,paste0('RNA_snn_res.',res)]
# seurat_Epithe$cellSubtype_expr = Recode(seurat_Epithe$cellSubtype_expr, "c(0,1,2,3,4,10)=1;8=2;9=3;11=4")
seurat_Epithe$cellSubtype_expr = factor(paste0('Expr-', seurat_Epithe$cellSubtype_expr), levels = paste0('Expr-', sort(unique(as.numeric(seurat_Epithe$cellSubtype_expr)))), ordered = T)
# save(seurat_Epithe, file = paste0(workdir, "/1.4.2_seurat_Epithe_scaled_", sample, ".RData"))

#使用infercnv的label标记expr的umap图
P_expr = DimPlot(seurat_Epithe, reduction = "umap", pt.size = 0.2, group.by = "cellSubtype_expr", label = TRUE, cols = c(brewer.pal(n = 9, name = "Set1"), brewer.pal(n = 8, name = "Set2"))[as.numeric(strsplit2(levels(seurat_Epithe$cellSubtype_expr), split='-')[,2]) +1])
pdf(file = paste0(workdir, "/1.4.2_umap_exprCluster.pdf"),width = 7, height = 6)
print(P_expr)
dev.off()
tiff(file = paste0(workdir, "/1.4.2_umap_exprCluster.tiff"), res=300, width = 7, height = 6, compression="lzw", units="in")
print(P_expr)
dev.off()

#根据infercnv的标签画umap图
P_infercnv = DimPlot(seurat_Epithe, reduction = "umap", pt.size = 0.2, group.by = "cellSubtype_infercnv", label = TRUE, cols = c(brewer.pal(n = 9, name = "Set1"), brewer.pal(n = 8, name = "Set2"))[as.numeric(strsplit2(levels(seurat_Epithe$cellSubtype_infercnv), split='-')[,2])])

pdf(file = paste0(workdir, "/1.4.2_umap_cnvCluster.pdf"),width = 7, height = 6)
print(DimPlot(seurat_Epithe, reduction = "umap", pt.size = 0.2, group.by = "cellSubtype_infercnv", label = TRUE, cols = c(brewer.pal(n = 9, name = "Set1"), brewer.pal(n = 8, name = "Set2"))[as.numeric(strsplit2(levels(seurat_Epithe$cellSubtype_infercnv), split='-')[,2]) %>% replace(11,12)]))
print(P_infercnv)
dev.off()
#res: 分辨率，单位是每英寸的像素数，units是单位，in是英寸；compression: 图像的压缩类型，lzw无损数据压缩算法，广泛应用于图像和文本等数据的压缩
tiff(file = paste0(workdir, "/1.4.2_umap_cnvCluster.tiff"), res=300, width = 7, height = 6, compression="lzw", units="in")
print(P_infercnv)
dev.off()
save(seurat_Epithe, file = paste0(workdir, "/1.4.2_seurat_Epithe_scaled_", sample, ".RData"))


##统计信息 cnv和expr以及reference还有上皮的数量
# load(paste0(workdir, "/1.4.2_seurat_Epithe_scaled_", sample, ".RData"))
sum(table(seurat_Epithe@meta.data$cellTypes_new)) #单独样本的上皮的数量
a <- sum(table(seurat_i@meta.data$cellTypes_new))-sum(table(seurat_Epithe@meta.data$cellTypes_new))
print(a) #对照的数量
table(seurat_Epithe$cellSubtype_infercnv) #每个分出来的cnv的细胞个数
table(seurat_Epithe$cellSubtype_expr) #每个根据findcluster函数找出来的聚类的细胞数量

#将其所有上皮命名为Unselected，然后将单个样本根据cnv种类的不同将单个样本中分成正常和肿瘤等等，其他样本还是Unselected
load('/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC/Results/Futher_Legend/Epithelial/Data/seurat_harmony/1.1.4_seurat_merge_Epithelial_harmony.RData')
seurat_merge4 <- seurat_merge
seurat_merge4@meta.data <- data.frame(seurat_merge4@meta.data, cellTypes = c("Unselected"))

seurat_merge4@meta.data[rownames(seurat_Epithe@meta.data[seurat_Epithe@meta.data$cellSubtype_infercnv %in% c("CNV-1"),]) , "cellTypes"] = "HPCP2-tumor"
seurat_merge4@meta.data[rownames(seurat_Epithe@meta.data[seurat_Epithe@meta.data$cellSubtype_infercnv == "CNV-2",]) , "cellTypes"] = "HPCP2-Normal_1"
seurat_merge4@meta.data[rownames(seurat_Epithe@meta.data[seurat_Epithe@meta.data$cellSubtype_infercnv == "CNV-3",]) , "cellTypes"] = "HPCP2-Normal_2"

# seurat_merge4@meta.data[rownames(seurat_Epithe@meta.data[seurat_Epithe@meta.data$cellSubtype_infercnv %in% c("CNV-20"),]) , "cellTypes"] = paste0(sample , "-Normal")
# seurat_merge4@meta.data[rownames(seurat_Epithe@meta.data[seurat_Epithe@meta.data$cellSubtype_infercnv %in% c("CNV-1","CNV-2","CNV-3","CNV-4","CNV-5","CNV-6","CNV-7","CNV-8","CNV-9","CNV-10","CNV-11","CNV-12","CNV-13","CNV-14","CNV-15","CNV-16","CNV-17","CNV-18","CNV-19","CNV-20") , ]) , "cellTypes"] =  paste0(sample , "-tumor")  # "HBCP5A-tumor"


#画单个样本在整个上皮的表达图
pdf(file = paste0(workdir,"/umap_cluster.pdf"),width = 7, height = 6) 
print(DimPlot(seurat_merge4, reduction = "umap", pt.size = 0.2, group.by = "cellTypes", label = TRUE,cells.highlight = list("HBCP10-tumor_1"=rownames(seurat_Epithe@meta.data[seurat_Epithe@meta.data$cellSubtype_infercnv %in% c("CNV-1") , ]) )), cols.highlight = list("HBCP10-tumor_1" = "red")) #cols= c("blue","grey","red")    cols= c('HBCP4B-tumor'= "blue",'Unselected'= "grey")

print(DimPlot(seurat_merge4, reduction = "umap", pt.size = 0.2, group.by = "cellTypes", label = TRUE,cells.highlight = list("HBCP10-Normal"=rownames(seurat_Epithe@meta.data[seurat_Epithe@meta.data$cellSubtype_infercnv %in% c("CNV-2") , ]) )), cols.highlight = list("HBCP10-Normal" = "red")) #cols= c("blue","grey","red")    cols= c('HBCP4B-tumor'= "blue",'Unselected'= "grey")
# print(DimPlot(seurat_merge4, reduction = "umap", pt.size = 0.2, group.by = "cellTypes", label = TRUE,cols= c("blue","red","grey"))) #cols= c("blue","grey","red")    cols= c('HBCP4B-tumor'= "blue",'Unselected'= "grey")
dev.off()


#防止颜色覆盖,导致颜色出现不全
workdir=paste0(topdir, '/Results/1.4.2_analyze_inferCNV/ref_T_B_macrophage_endothelial/', sample)
setwd(workdir)
#读入全部上皮数据
load('/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/BC/R28BC/seurat_merge_Epithelial.RData')
seurat_merge4 <- seurat_merge3
seurat_merge4@meta.data <- data.frame(seurat_merge4@meta.data, cellTypes = c("Unselected"))
#读入样本scaled后的数据
sample ="HBCP15"
load(paste0(workdir, "/1.4.2_seurat_Epithe_scaled_", sample, ".RData"))
table(seurat_Epithe$cellSubtype_infercnv) 

#查看cnv分布
# seurat_merge4@meta.data[rownames(seurat_Epithe@meta.data[seurat_Epithe@meta.data$cellSubtype_infercnv == "CNV-1",]) , "cellTypes"] = "HBCP15-tumor"
# seurat_merge4@meta.data[rownames(seurat_Epithe@meta.data[seurat_Epithe@meta.data$cellSubtype_infercnv == "CNV-2",]) , "cellTypes"] = "HBCP15-Normal"
seurat_merge4@meta.data[rownames(seurat_Epithe@meta.data[seurat_Epithe@meta.data$cellSubtype_infercnv %in% c("CNV-1","CNV-2","CNV-3","CNV-4","CNV-5","CNV-6","CNV-7","CNV-8","CNV-9","CNV-10") , ]) , "cellTypes"] = paste0(sample , "-tumor") #  ,"CNV-11","CNV-12","CNV-13","CNV-14","CNV-15"

#画单个样本在整个上皮的表达图
pdf(file = paste0(workdir,"/umap_cluster.pdf"),width = 7, height = 6) 
print(DimPlot(seurat_merge4, reduction = "umap", pt.size = 0.1, group.by = "cellTypes", label = TRUE,cells.highlight = list("HPCP3-tumor"=rownames(seurat_Epithe@meta.data[seurat_Epithe@meta.data$cellSubtype_infercnv %in% c("CNV-1"), ]) , "HPCP3-Normal_1" = rownames(seurat_Epithe@meta.data[seurat_Epithe@meta.data$cellSubtype_infercnv %in% "CNV-2" ,] ),"HPCP3-Normal_2" = rownames(seurat_Epithe@meta.data[seurat_Epithe@meta.data$cellSubtype_infercnv %in% "CNV-3" ,] )), cols.highlight = list("HPCP3-tumor"= 'red' , "HPCP3-Normal_1"= "blue","HPCP2-Normal_3"= "#9af764")))
dev.off()



#######
seurat_merge1 <- subset(seurat_Epithe,subset = cellSubtype_infercnv %in% c("CNV-3"))
seurat_merge2 <- subset(seurat_Epithe,subset = cellSubtype_infercnv %in% c("CNV-2","CNV-4"))
pdf(file = "/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/BC/R28BC_cnv_1/Results/1.4.2_analyze_inferCNV/ref_T_B_macrophage_endothelial/HBCP19.1/feature_1.pdf",,width = 15, height = 14)
print(VlnPlot(seurat_Epithe, features = c("PPARG","FGFR3","ERBB2","EGFR","MKI67","GATA3","TP63","CD44","KRT7","HER2","KRT20","TP53"),group.by = "cellSubtype_infercnv"))
dev.off()








#对整个上皮细胞进行注释
load('/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/BC/R28BC_cnv_1/Results/1.4.2_analyze_inferCNV/ref_T_B_macrophage_endothelial/HBCP19/1.4.2_seurat_Epithe_scaled_HBCP19.RData')
load('/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/BC/R28BC_cnv_1/Results/1.4.2_analyze_inferCNV/ref_T_B_macrophage_endothelial/HBCP22/1.4.2_seurat_Epithe_scaled_HBCP22.RData')
seurat_merge4 <- seurat_merge3
seurat_merge4@meta.data <- data.frame(seurat_merge4@meta.data, cellTypes = c("Tumor"))
length(seurat_merge4@meta.data[rownames(seurat_merge3@meta.data[,seurat_merge3@meta.data$RNA_snn_res.1.3 %in% c("27")])])
a <- rownames(seurat_merge3@meta.data[seurat_merge3@meta.data$RNA_snn_res.1.3 == "27",])
b <- rownames(seurat_merge4@meta.data[seurat_merge4@meta.data$cellTypes == "Normal",])
c <- intersect(a,b)
seurat_merge4@meta.data[rownames(seurat_Epithe@meta.data[seurat_Epithe@meta.data$cellSubtype_infercnv %in% c("CNV-2","CNV-4"),]) , "cellTypes"] = "Undefined"
seurat_merge4@meta.data[c , "cellTypes"] = "Normal"
seurat_merge4$cellTypes <- Idents(seurat_merge4)
table(seurat_merge4$cellTypes)
seurat_merge4$cellTypes <- as.factor(seurat_merge4$cellTypes)

pdf(file = "/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/BC/R28BC/Results/Futher_Legend/umap_cluster3.1.pdf",width = 7, height = 6)
print(DimPlot(seurat_merge4, reduction = "umap",pt.size = 0.2, group.by = "cellTypes",cols = c("red","#75bbfd","grey")))#,cols = c("grey","red")
print(DimPlot(seurat_merge4, reduction = "umap", pt.size = 0.2, group.by = "cellTypes",cells.highlight = list("Normal"=rownames(seurat_merge4@meta.data[seurat_merge4$cellTypes %in% c("Normal") , ]) )), cols.highlight = list("Normal" = "red"))
dev.off()
save(seurat_merge4,file = '/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/BC/R28BC/seurat_merge_Epithelial_annot.RData')
load('/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/BC/R28BC/seurat_merge_Epithelial_annot.RData')

