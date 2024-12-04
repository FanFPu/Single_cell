from pathlib import Path
import os
import pandas as pd
import scanpy as sc
import infercnvpy as cnv
import matplotlib.pyplot as plt
sc.settings.set_figure_params(figsize=(5,5))


OUTPUT_DIR='/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/BC/Mouse/03.inferCNV/24w'
Path(OUTPUT_DIR).mkdir(parents=True,exist_ok=True)
os.chdir(OUTPUT_DIR)
os.getcwd()    

#读取原始数据
adata_raw = sc.read_h5ad('/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/BC/Mouse/03.inferCNV/24w/24w.h5ad')

#以下两者就是统一名字
# CellTypeKey="CellTypeS2"
# CNVCellTypeKey="CellTypeS3"

#为adata.var添加基因信息
cnv.io.genomic_position_from_gtf("/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/software/gencode.vM23.annotation.gtf.gz",adata_raw)
# adata = adata_raw
adata_raw.obs['celltype'].value_counts()

adata_raw.var.loc[:, [ 'gene_id','chromosome','start','end']].head()
# "Adipocyte","B cell","DC","Endothelial cell","Fibroblast","Macrophage","Neurocyte","Plasma cell","Smooth muscle cell","T cell","Normal epithelial cell"

#画数据的umap图
# plt.figure(figsize=(8, 6))
# plt.figure(figsize=(9,6), dpi=300)
sc.pp.neighbors(adata_raw, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata_raw)
sc.pl.umap(adata_raw, color="celltype")
plt.savefig('sc_umap.pdf', dpi=300, bbox_inches = 'tight')

# # rc_context用于指定figure大小
# with rc_context({'figure.figsize': (4, 4)}):
#     sc.pl.umap(adata, color='celltype')
# plt.savefig("sc_umap.pdf", format="pdf")

#进行cnv
# np.unique(adata_raw.obs["celltype"])
cnv.tl.infercnv(
    adata_raw,
    reference_key= "celltype",
    reference_cat=['B cell', 'Endothelial cell', 'Plasma cell', 'T cell',"Macrophage","DC"],
    window_size=250
    ,n_jobs=12
    )
#各类细胞在染色体上的基因表达量
cnv.pl.chromosome_heatmap(adata_raw, groupby="celltype")
plt.savefig('chromosome_heatmap.pdf', dpi=300, bbox_inches = 'tight')


#细胞聚类和细胞注释时的聚类类似，只是使用的数据时adata.obsm'X_cnv'中的包括3步算PCA、算距离neighbors、leiden聚类
#各cluster在染色体上的基因表达量图，前几个cluster没有显著差异的表达的基因组区域。显著差异的表达的基因组区域可能是由于CNV，因此这些cluster可能是肿瘤细胞。（前几#个cluster可能是肿瘤细胞）

cnv.tl.pca(adata_raw)
cnv.pp.neighbors(adata_raw)
cnv.tl.leiden(adata_raw)
plt.figure(figsize=(10, 8))
cnv.pl.chromosome_heatmap(adata_raw, groupby="cnv_leiden", dendrogram=True)
plt.savefig('chromosome_heatmap_leiden.pdf', dpi=300, bbox_inches = 'tight')

#CNV分数
#Clusters有高的分数说明其更可能受到CNV影响，更可能是肿瘤细胞。
#9、23分数较高而且没有聚类在一起，可能是肿瘤细胞。
cnv.tl.umap(adata_raw)
cnv.tl.cnv_score(adata_raw)

#画cnv的cnv_umap
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(11, 11))
ax4.axis("off")
cnv.pl.umap(
    adata_raw,
    color="cnv_leiden",
    legend_loc="on data",
    legend_fontoutline=2,
    ax=ax1,
    show=False,
)
cnv.pl.umap(adata_raw, color="cnv_score", ax=ax2, show=False)
cnv.pl.umap(adata_raw, color="celltype", ax=ax3)
plt.savefig('cnv_umap_score.pdf', dpi=300, bbox_inches = 'tight')


#画cnv的sc_umap的打分
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 11), gridspec_kw=dict(wspace=0.5))
ax4.axis("off")
sc.pl.umap(adata_raw, color="cnv_leiden", ax=ax1, show=False)
sc.pl.umap(adata_raw, color="cnv_score", ax=ax2, show=False)
sc.pl.umap(adata_raw, color="celltype", ax=ax3)
plt.savefig('sc_umap_score.pdf', dpi=300, bbox_inches = 'tight')


#基于这些观察，将细胞分配给“肿瘤”或“正常”。在 adata.obs 中添加一个新列 cnv_status。
adata_raw.obs["cnv_status"] = "normal"
adata_raw.obs.loc[adata_raw.obs["cnv_leiden"].isin(["5","0"]), "cnv_status"] = "tumor"

#将定义的正常和肿瘤进行UMAP展示
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), gridspec_kw=dict(wspace=0.5))
cnv.pl.umap(adata_raw, color="cnv_status", ax=ax1, show=False)
sc.pl.umap(adata_raw, color="cnv_status", ax=ax2)
plt.savefig("cnv_status.pdf", format="pdf", dpi=300, bbox_inches = 'tight')

#根据聚类后画出正常和肿瘤的染色体变化热图
cnv.pl.chromosome_heatmap(adata_raw[adata_raw.obs["cnv_status"] == "tumor", :],)
plt.savefig("chromosome_heatmap_tumor.pdf", format="pdf", dpi=300, bbox_inches = 'tight')
cnv.pl.chromosome_heatmap(adata_raw[adata_raw.obs["cnv_status"] == "normal", :])
plt.savefig("chromosome_heatmap_normal.pdf", format="pdf", dpi=300, bbox_inches = 'tight')

adata_raw.write_h5ad('adata_raw.h5ad')

















































