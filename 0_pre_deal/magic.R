#qsub -cwd -l vf=200g,p=40 -P P21H28400N0232 -q st.q,st_supermem.q work.sh
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
library(sctransform)

#Work directory
indir = '/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC/Sample'
topdir = '/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC'
workdir = topdir
setwd(workdir)


#多线程分析
plan()
plan("multisession", workers = as.numeric(2))
# Enable parallelization
plan()
options(future.globals.maxSize= 429496729600)

##Merge
pids = read.delim('SID-PID.txt', col.names=1)
pid_0 = rownames(pids)[1]

#load(paste0(indir, '/', '1.1.4_seurat_merge_raw.RData'))
#首先获得一个数据作为第一个seurat_merge
load(paste0(indir, '/', pid_0, '/seurat_', pid_0, '_single.RData'))  #里面的对象是seurat_comb_singlet
seurat_merge = seurat_comb_singlet
seurat_merge$Patients = pids[pid_0,1]


for (sample_id in rownames(pids)[-1]) {
	print(sample_id)
	load(paste0(indir, '/', sample_id, '/seurat_', sample_id, '_single.RData'))
  seurat_comb_singlet$Patients = pids[sample_id,1]
	print('begin merge...')
	seurat_merge = merge(seurat_merge, seurat_comb_singlet)
	print(c(table(seurat_merge$orig.ident), table(seurat_merge$Patients), dim(seurat_merge)))
}
save(seurat_merge,file = "/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC/Data/rawData/1.1.4_seurat_merge_raw_NO_QC.RData")

#QC after Merger，合并后的质控
maxGenes <- as.numeric(8000) #8000
minUMIs <- as.numeric(1000) #1000
minGenes <- as.numeric(500) #500
maxPercent.mt <- as.numeric(10) #10
seurat_merge <- subset(seurat_merge, subset = nCount_RNA>minUMIs & nFeature_RNA > minGenes & nFeature_RNA < maxGenes & percent.mt < maxPercent.mt)
save(seurat_merge, file = paste0(workdir, '/1.1.4_seurat_merge_raw.RData'))

load('1.1.4_seurat_merge_raw.RData')
if(!dir.exists(paste0(workdir, "/Results/1.1.4_seurat_notharmony"))){dir.create(paste0(workdir, "/Results/1.1.4_seurat_notharmony"), recursive = TRUE)}

#质控图
pdf(file = "Results/1.1.4_seurat_notharmony/1.1.4_qc_post.pdf", width = 30, height = 6)
plot1 <- FeatureScatter(seurat_merge, feature1 = "nCount_RNA", feature2 = "percent.mt")#一般呈现负相关比较好
plot2 <- FeatureScatter(seurat_merge, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")#一般呈现正相关比较好，feature代表捕获的gene的数量，越多代表转录活性越高
plot3 <- FeatureScatter(seurat_merge, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "Patients")
plot4 <- FeatureScatter(seurat_merge, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "Patients")
print(VlnPlot(seurat_merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "Patients"))
print(VlnPlot(seurat_merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3),group.by ="orig.ident")
print(plot1 + plot2)
print(plot3 + plot4)
dev.off()

##normalize,标准化
# load('1.1.4_seurat_merge_raw.RData')
seurat_merge <- NormalizeData(seurat_merge)

##seurat_merge，去基因
print(dim(seurat_merge))
removeBiasGenes <- function(genes){
  RPgenes <- genes[intersect(grep("^RP", genes), grep("-", genes))]
  RPgenes2 <- genes[grep("^RP[SL]", genes)] #后面跟S或L的gene
  MTgenes <- genes[grep("^MT-", genes)]
  # MTgenes2 <- genes[grep("^MTRNR", genes)]
  # HSPgenes <- genes[grep("^HSP[ABDE]",genes)]
  # DNAgenes <- genes[grep("^DNAJ",genes)]
  CTCgenes <- genes[intersect(grep("^CTC", genes), grep("-", genes))]
  MIRgenes <- genes[grep("^MIR", genes)]
  ACgenes <- genes[intersect(grep("^AC[0-9]", genes), grep(".", genes))]
  CTgenes <- genes[intersect(grep("^CT", genes), grep("-", genes))]
  LINCgenes <- genes[grep("^LINC[0-9]", genes)]
  ALgenes <- genes[intersect(grep("^AL", genes), grep(".", genes))]
  HEMOgene <- c("HBB", "HBA1", "HBA2")
  # Othergene <- c('ATF3','BTG2','CEBPB','CEBPB-AS1','CEBPD','CXCL1','EGR1','FOS','FOSB','FOSL1','FOSL1P1','FOSL2','ID3',	'IER2','JUN','JUNB','JUND','MT1A','MT1B','MT1E','MT1F','MT1G','MT1H','MT1L','MT1M','MT1X','MT2A','NFKBIA','NR4A1',	'PPP1R15A','SOCS3','UBC','ZFP36')
  rmgenes <- c(RPgenes, RPgenes2, MTgenes, CTCgenes, MIRgenes, ACgenes, CTgenes, LINCgenes, ALgenes, HEMOgene,MTgenes2,HSPgenes,DNAgenes,Othergene)
  return(rmgenes)
}
rmgenes = removeBiasGenes(rownames(seurat_merge))
seurat_merge = seurat_merge[setdiff(rownames(seurat_merge), rmgenes),]
save(seurat_merge, file = '1.1.4_seurat_merge_normalize.RData')



#magic处理
library(Rmagic)
library(reticulate)
# use_python("/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/software/miniconda/envs/miniconda/bin/python")
use_python("/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/miniconda/envs/gch_base/bin/python")
seurat_merge <- magic(seurat_merge)
seurat_merge@assays$RNA@data <- seurat_merge@assays$MAGIC_RNA@data


###Find Varaiable Gene —— Scale（缩放）—— reduction（PCA，scVI）
load('1.1.4_seurat_merge_normalize.RData')
