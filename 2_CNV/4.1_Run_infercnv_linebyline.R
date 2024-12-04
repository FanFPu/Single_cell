#### inferCNV copy number for each tumor samples, including all cells or only tumor cells####
#包不能在2-4，2-5，login-2上面安装，因为连不上网
library(rjags)
library(infercnv)
library(AnnoProbe)
library(tidyverse)
library(Seurat)
#library(parallel)
library(Cairo)
library(limma)
options(bitmapType = "cairo")

# args = commandArgs(T)
# sample_i = args[1] #P13

#多线程运行
plan()
plan("multisession", workers = as.numeric(2))
# Enable parallelization
plan()
options(future.globals.maxSize= 429496729600)

#工作路径
topdir = "/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC/SC_infercnv1"
workdir = "/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC/SC_infercnv1/test"
if(!dir.exists(workdir)){dir.create(workdir, recursive = TRUE)}
setwd(workdir)
#写入数据  #返回seurat_merge
# load("/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/BC/R28BC_cnv/saveData/seurat_merge_annot.RData")
# #保存文件路径
seurat_path = paste0(topdir, "/saveData/1.4.1_seurat_merge_infercnv.RData")  

# # seurat_merge <- subset(seurat_merge, subset = nCount_RNA > 10 & nCount_RNA < +Inf )  
# if (!file.exists(seurat_path)) {
#     #load(paste0(topdir, "/saveData/1.2.1_seurat_merge_annot.RData"))
#     seurat_merge_infercnv <- subset(seurat_merge, subset = cellTypes_new %in% c('Epithelial','Tcell', 'Bcell','Macrophage', 'Endothelial'))
#     seurat_lists = SplitObject(seurat_merge_infercnv, split.by = 'Patients') #按照split.by这一列进行拆分
#     save(seurat_lists, file = seurat_path)
# }

load(seurat_path)
for (sample in names(seurat_lists)) {
    if(!dir.exists(paste0(workdir, '/', sample))){dir.create(paste0(workdir, '/', sample), recursive = TRUE)}
}

##run inferCNV
sample ="HPCP1"

cat(sample)
options(bitmapType = "cairo")
seurat_i = seurat_lists[[sample]]
save(seurat_i, file=paste0(workdir, '/', sample, '/seurat_Epithe.RData'))



print(table(seurat_i$cellTypes_new))
gc()
mat_comb <- as.matrix(GetAssayData(seurat_i, assay="RNA",slot = "counts"))

geneinfo <- annoGene(rownames(mat_comb), "SYMBOL", "human")
# geneinfo <- annoGene(rownames(mat_comb), "SYMBOL", "mouse")
geneinfo <- geneinfo[with(geneinfo,order(chr,start)),c(1,4:6)]
geneinfo <- geneinfo[!duplicated(geneinfo[,1]),]

rownames(geneinfo) <- geneinfo$SYMBOL

## for inferCNV the geneorder file just have 3 colums: chr, star, end, there is no name colum
geneinfo$SYMBOL <- NULL
geneinfo = geneinfo %>% mutate(chr_n=strsplit2(geneinfo$chr, split='chr')[,2]) %>%
 mutate(chr_n = replace(chr_n, chr_n=='X', 23)) %>%
 mutate(chr_n = replace(chr_n, chr_n=='Y', 24)) %>%
 mutate(chr_n = replace(chr_n, chr_n=='M', 25)) %>%
 mutate(chr_n = as.numeric(chr_n)) %>%
 arrange(chr_n, start, end) %>%
 select(-chr_n)
head(geneinfo)


seurat_i$cluster <- Idents(seurat_i)
metada <- seurat_i@meta.data
# 删掉数量较少的细胞类型，防止被全部过滤掉
# metada <- metada[metada$cellTypes_new != "Bcell", ]  

metadata <- metada[,c("cellTypes_new", "cellTypes_new")]
colnames(metadata) <- c("cellType", "group")

head(metadata)

write.table(metadata, file=paste0(workdir,'/',sample,"/cellAnnotations.txt"), sep="\t", col.names = FALSE,quote = FALSE)
write.table(geneinfo, file=paste0(workdir,'/',sample,"/gene_ordering_file.txt"), sep = "\t", col.names = FALSE, quote = FALSE)

count_mat <- mat_comb[rownames(geneinfo),]
dim(count_mat)
gc()
infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = count_mat, 
                                     gene_order_file = paste0(workdir,'/',sample,"/gene_ordering_file.txt"), 
                                     delim = "\t",
                                     min_max_counts_per_cell = c(10, +Inf),
									 #  ref_group_names = c("Tcell", "Macrophage"),
                                     ref_group_names = c("Tcell","Bcell","Macrophage","Endothelial"),
                                    #  ref_group_names = c("Tcell","Macrophage","Endothelial"),
                                    #  ref_group_names = c("Macrophage","Endothelial"),
                                     annotations_file = paste0(workdir,'/',sample,"/cellAnnotations.txt"))

gc()
save(infercnv_obj, file = paste0(workdir,'/',sample,"/infercnv_obj1.RData"))
# perform infercnv operations to reveal cnv signal
##cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics，表示基因在所有参照细胞中，表达count的平均值的最小阈值，10X数据更稀疏，所以这个值小 
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir=sample,  # dir is auto-created for storing outputs  ,out_dir表示输出文件夹的名字，没有会自动创建
                             cluster_by_groups = F,   # cluster   cluster_by_groups：先区分细胞来源，再做层次聚类 
                             scale_data = FALSE,
                             denoise=T,
                             sd_amplifier = 1.0,
                             window_length = 150,
                             HMM=T,
                             #analysis_mode='subclusters',  #默认是"samples"
			                 num_threads=50,
                             output_format = "pdf" # "png"
                             )
save(infercnv_obj, file = paste0(workdir,'/',sample,"/infercnv_obj2.RData"))









#数据加载
mouse_geneOrderingFile<-read.table('/share/nas1/Data/DataBase/Genome/Mus_musculus/Transcriptome/Mm10/genes/genes.gtf', header=F, sep = '\t')
View(mouse_geneOrderingFile)

#列选择
mouse_geneOrderingFile1<-mouse_geneOrderingFile[, c(9, 1, 4, 5)]

#基因名提取
mouse_geneOrderingFile1$V9<-gsub('.*gene_name ', '', mouse_geneOrderingFile1$V9)
mouse_geneOrderingFile1$V9<-gsub(';.*', '', mouse_geneOrderingFile1$V9)

head(mouse_geneOrderingFile1)

##基因排序
View(mouse_geneOrderingFile1)

mouse_geneOrderingFile2<-subset(mouse_geneOrderingFile1, subset=V1%in%c(1:19, "MT", "X", "Y" ))

dim(mouse_geneOrderingFile1)
dim(mouse_geneOrderingFile2)


order_df<-data.frame(chr=c(1:19,'MT', 'X', 'Y'), order=c(1:22), stringsAsFactors = F)

mouse_geneOrderingFile2$order<-mapvalues(mouse_geneOrderingFile2$V1, as.character(order_df$chr), as.character(order_df$order))

mouse_geneOrderingFile2$order<-as.numeric(mouse_geneOrderingFile2$order)
View(mouse_geneOrderingFile2)


mouse_geneOrderingFile2<-mouse_geneOrderingFile2[sort(mouse_geneOrderingFile2$order, index.return=T)$ix,]
View(mouse_geneOrderingFile2)
mouse_geneOrderingFile2$order<-NULL

##去重
mouse_geneOrderingFile3<-mouse_geneOrderingFile2[!duplicated(mouse_geneOrderingFile2$V9),] 

##保存
write.table(mouse_geneOrderingFile3, file = '/share/nas1/Data/Users/luohb/tools_20201026/mouse_inferCNV/mm_geneOrderingFile.txt',quote = F, sep = '\t', row.names = F, col.names = F)







