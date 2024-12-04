#### inferCNV copy number for each tumor samples, including all cells or only tumor cells####
#包不能在2-4，2-5，login-2上面安装，因为连不上网,只能在softwware上面安装

###整条染色体或大片段染色体的增加或丢失(gain or deletions)。
###工作原理是:以一组"正常"细胞作为参考，分析肿瘤基因组上各个位置的基因表达量强度变化. 通过热图的形式展示每条染色体上的基因相对表达量，相对于正常细胞，肿瘤基因组总会过表达或者低表达。
library(rjags)
library(infercnv)
library(AnnoProbe)
library(tidyverse)
library(Seurat)
library(parallel)
library(Cairo)
library(limma)
options(bitmapType = "cairo")


args = commandArgs(T)
sample_i = args[1] #P13

library(rjags)
library(infercnv)
library(AnnoProbe)
library(tidyverse)
library(Seurat)
library(parallel)
library(Cairo)
library(limma)
options(bitmapType = "cairo")


#工作路径
topdir = "/jdfsbjcas1/ST_BJ/P21H28400N0232/xieguixiang/BLCA/SC/Mouse/yangzongzheng/infercnv/6w_test"
workdir = "/jdfsbjcas1/ST_BJ/P21H28400N0232/xieguixiang/BLCA/SC/Mouse/yangzongzheng/infercnv/6w_test"
if(!dir.exists(workdir)){dir.create(workdir, recursive = TRUE)}
setwd(workdir)
#写入数据
# load("/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/BC/R23BC/seurat_merge_fibroblast.RData")
# load("/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/BC/R28BC/seurat_merge_Immune.RData")
# load('/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/BC/R28BC/seurat_merge_Epithelial.RData')
# seurat_merge =merge(seurat_merge1,y=seurat_merge2)
# seurat_merge = merge(seurat_merge,y=seurat_merge3)
# table(seurat_merge$cellTypes_new)
##normalize,标准化
load('/jdfsbjcas1/ST_BJ/P21H28400N0232/xieguixiang/BLCA/SC/Mouse/yangzongzheng/infercnv/6w/seurat_Epithe.RData')
seurat_merge <- seurat_i
seurat_merge <- NormalizeData(seurat_merge)
save(seurat_merge, file = paste0(topdir, "/saveData/seurat_merge_annot.RData"))
#读入数据，全部样本的注释信息
# load(paste0(topdir, "/saveData/seurat_merge_annot.RData"))
#保存文件路径
seurat_path = paste0(topdir, "/saveData/1.4.1_seurat_merge_infercnv.RData")  #按照样本进行分开，里面具有每个样本的细胞类型

#将数据筛选，其次按照样本进行分开，按照样本进行分析
if (!file.exists(seurat_path)) {
    #load(paste0(topdir, "/saveData/1.2.1_seurat_merge_annot.RData"))
    seurat_merge_infercnv <- subset(seurat_merge, subset = cellTypes_new %in% c('Epithelial','Tcell', 'Bcell','Macrophage', 'Endothelial'))
    seurat_lists = SplitObject(seurat_merge_infercnv, split.by = 'Patients') #按照split.by这一列进行拆分
    save(seurat_lists, file = seurat_path)
}
# table(seurat_merge_infercnv$cellTypes_new) #查看merge后的信息是否subset正确

load(seurat_path)

for (sample in names(seurat_lists)) {
    if(!dir.exists(paste0(workdir, '/', sample))){dir.create(paste0(workdir, '/', sample), recursive = TRUE)}
}

##run inferCNV
run_infercnv = function(sample) {
  cat(sample)
  options(bitmapType = "cairo")  #用于指定绘图设备（device）使用 Cairo 图形库进行绘制。Cairo 是一个2D图形库，用于绘制矢量图形，支持多种输出格式，包括PDF、PNG等，前提安装了cairo图形库
  seurat_i = seurat_lists[[sample]]
  save(seurat_i, file=paste0(workdir, '/', sample, '/seurat_Epithe.RData'))
  print(table(seurat_i$cellTypes_new))
  # rm(seurat_lists)
  gc() # gc函数用来执行垃圾收集，自动处理内存管理
  mat_comb <- as.matrix(GetAssayData(seurat_i, assay="RNA",slot = "counts"))
  
  geneinfo <- annoGene(rownames(mat_comb), "SYMBOL", "human")
  geneinfo <- geneinfo[with(geneinfo,order(chr,start)),c(1,4:6)] #with(对象，命令)
  geneinfo <- geneinfo[!duplicated(geneinfo[,1]),]
  rownames(geneinfo) <- geneinfo$SYMBOL#意思应该是将geneinfo中SYMBOL的信息作为行名，SYMBOL里是基因名
  
  ## for inferCNV the geneorder file just have 3 colums: chr, star, end, there is no name colum
  geneinfo$SYMBOL <- NULL
  geneinfo = geneinfo %>% mutate(chr_n=strsplit2(geneinfo$chr, split='chr')[,2]) %>%  #mutate函数是新增一列
    mutate(chr_n = replace(chr_n, chr_n=='X', 23)) %>%#对已有列的值进行变换、计算，或者添加新的列
    mutate(chr_n = replace(chr_n, chr_n=='Y', 24)) %>%
    mutate(chr_n = replace(chr_n, chr_n=='M', 25)) %>% #将chr_n=M换成了chr_n=25
    mutate(chr_n = as.numeric(chr_n)) %>% #将其换成数值型
    arrange(chr_n, start, end) %>%  #进行排序，默认是升序，降序需加上desc()
    select(-chr_n) #此函数可以删除，选择，重命名
  #head(geneinfo)
  
  seurat_i$cluster <- Idents(seurat_i) #将其细胞类型筛选出来新增到cluster这一列中
  metada <- seurat_i@meta.data
  # metada <- metada[metada$cellTypes_new != "Endothelial", ]   # 删掉数量较少的细胞类型，防止被全部过滤掉
  metadata <- metada[,c("cellTypes_new", "cellTypes_new")]
  colnames(metadata) <- c("cellType", "group")
  #head(metadata)
  
  #quote 的默认值是 TRUE，表示字符型数据将被引号包裹。这有助于确保在读取文件时能够正确地区分字符型数据，特别是当数据中包含分隔符或特殊字符时
  write.table(metadata, file=paste0(workdir,'/',sample,"/cellAnnotations.txt"), sep="\t", col.names = FALSE,quote = FALSE)
  write.table(geneinfo, file=paste0(workdir,'/',sample,"/gene_ordering_file.txt"), sep = "\t", col.names = FALSE, quote = FALSE)
  
  count_mat <- mat_comb[rownames(geneinfo),]
  rm(mat_comb)
  rm(seurat_i)
  rm(metada)
  gc()
  
  infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = count_mat, 
                                       gene_order_file = paste0(workdir,'/',sample,"/gene_ordering_file.txt"), 
                                       delim = "\t", min_max_counts_per_cell = c(10, +Inf),
                                       #  ref_group_names = c("Tcell", "Macrophage", "Endothelial"),
                                        ref_group_names = c("Tcell", "Bcell", "Macrophage", "Endothelial"),
                                      #  ref_group_names = c("Tcell","Macrophage"),
                                       #  ref_group_names = c("Endothelial","Macrophage"),
                                       annotations_file = paste0(workdir,'/',sample,"/cellAnnotations.txt"))
  
  rm(count_mat)
  gc()
  save(infercnv_obj, file = paste0(workdir,'/',sample,"/infercnv_obj1.RData"))
  # perform infercnv operations to reveal cnv signal
#sd_amplifier 的作用是调整差异性基因的标准差阈值，以影响差异性基因的检测结果。增加 sd_amplifier 的值会导致更多的基因被认为是差异性基因，而减小它的值会减少差异性基因的数量
  infercnv_obj = infercnv::run(infercnv_obj,
                               cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics,过滤低表达基因
                               out_dir=sample,  # dir is auto-created for storing outputs
                               cluster_by_groups=F,   # cluster
                               scale_data=FALSE,
                               denoise=T, #denoise：默认FALSE，对CNV矩阵进行降噪
                               sd_amplifier=1.5,  
                               HMM=T, #隐马尔可夫模型，运行它会导致时间变长
                               #analysis_mode='subclusters',
                               output_format = "png",  #如果细胞数量过大，会导致png图片出现空白，那么就需要保存为pdf格式
                               num_threads=50 #线程数
  )
  save(infercnv_obj, file = paste0(workdir,'/',sample,"/infercnv_obj2.RData"))
#### the medianfitered process for inferCNV
# infercnv_obj_medianfiltered = infercnv::apply_median_filtering(infercnv_obj)
# rm(infercnv_obj)
# save(infercnv_obj_medianfiltered, file = "infercnv_obj_medianfiltered.RData")
}

#clus <- makeCluster(4)
#clusterExport(clus,c('workdir', 'seurat_lists'),envir = environment())
#clusterEvalQ(clus, library('Seurat'))
#clusterEvalQ(clus, library('infercnv'))
#clusterEvalQ(clus, library('AnnoProbe'))
#clusterEvalQ(clus, library('tidyverse'))
#clusterEvalQ(clus, library('Cairo'))
#parLapply(clus, names(seurat_lists), fun = run_infercnv)
#stopCluster(clus)

run_infercnv("HBCP21")










