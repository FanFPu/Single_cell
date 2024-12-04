# annotation for BCa_immune_cell
import os
import gc
import glob
import pandas as pd
import numpy as np
import scanpy as sc
import scanpy.external as sce
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from anndata import read_csv
from loompy import make_row_attrs_from_gene_annotations
from matplotlib.pyplot import rc_context
import warnings
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import scanpy as sc
import numpy as np
import os

warnings.filterwarnings ('ignore')
np.random.seed(42)
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Arial']
plt.rcParams['font.size'] = 6  
sc.set_figure_params(dpi=300, color_map="gist_heat_r")
os.chdir('/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/Paper/Single_cell/4.Epi_sc_ana/')
keys = sorted([f"leiden_res_{res:4.2f}" for res in np.arange(0.1, 2.1, 0.1)])


adata = sc.read_h5ad('/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/Paper/Single_cell/4.Epi_sc_ana/Normal_Epi_pyharmony/Nor_Epi_pyharmony_new.h5ad') # load data

adata1 = sc.read('/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/Paper/Single_cell/4.Epi_sc_ana/Tumor_Epi_pyharmony/Tumor_Epi_pyharmony_new.h5ad')

adata = adata.concatenate(adata1, 
)

adata.obs_names_make_unique()
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
sc.tl.pca(adata)
# sce.pp.harmony_integrate(adata, 'orig.ident')
# adata.obsm['X_pca'] = adata.obsm['X_pca_harmony']
sc.pp.neighbors(adata)
for res in np.arange(0.1, 2.1, 0.1):
    sc.tl.leiden(
        adata, key_added=f"leiden_res_{res:4.2f}", resolution=res
    )
sc.tl.umap(adata)
# 输出文件夹路径
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import os

# 输出文件夹路径
output_dir = "/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/Paper/Single_cell/4.Epi_sc_ana/"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# 设置全局字体
plt.rcParams["font.family"] = "sans-serif"
# 创建一个 PDF 文件保存所有图像
pdf_filename = os.path.join(output_dir, "All_UMAP_Leiden_Resolutions.pdf")
with PdfPages(pdf_filename) as pdf:
    for res in np.arange(0.1, 2.1, 0.1):
        leiden_key = f"leiden_res_{res:4.2f}"
        sc.pl.umap(
            adata,
            color=leiden_key,
            title=f"Leiden Clustering (Resolution={res:4.2f})",
            legend_loc="on data",
            frameon=True,
            show=False
        )
        plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
        fig = plt.gcf()
        fig.set_size_inches(10, 8)  # 设置更大的画布尺寸
        
        # 添加轴标签
        ax = plt.gca()
        ax.set_xlabel("UMAP 1")
        ax.set_ylabel("UMAP 2")
        
        # 保存当前图到PDF
        pdf.savefig(fig)
        plt.close()  # 关闭当前图，避免重复显示


# 定义输出目录
output_dir = "你的输出路径"  # 请替换为实际的输出路径
pdf_filename = os.path.join(output_dir, "Patients.pdf")
# 打开PDF文件，将所有 Patients 分组的图放入同一个PDF文件
with PdfPages(pdf_filename) as pdf:
    sc.pl.umap(
        adata,
        color="Patients",  # 使用 Patients 列进行着色
        title="UMAP Clustering by Patients",
        legend_loc="on data",
        frameon=True,
        show=False
    )
    plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
    fig = plt.gcf()
    fig.set_size_inches(10, 8)  # 设置更大的画布尺寸
    
    # 添加轴标签
    ax = plt.gca()
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    pdf.savefig(fig)
    plt.close()






# adata = sc.read_h5ad('Nor_Epi_pyharmony.h5ad')
res = 0.30
adata.obs["CellType"] = adata.obs[f"leiden_res_{res:4.2f}"].map(
    {"0" : "Nor_type1",
     "1" : "Nor_type3",
     "2":  "Nor_type2",
     "3":  "Nor_type4",
     "4": "Nor_type5"}
)
adata.obs['CellType'] = adata.obs['CellType'].cat.add_categories('unknown')
adata.obs['CellType'] = adata.obs['CellType'].fillna('unknown')
adata.write("Nor_Epi_pyharmony.h5ad") # cluster 5 have mix signal
adata.obs.to_csv('Nor_Epi_pyharmony_annotation.csv')



pdf_filename = os.path.join(output_dir, "CellType1.pdf")
# 打开PDF文件，将所有 Patients 分组的图放入同一个PDF文件
with PdfPages(pdf_filename) as pdf:
    sc.pl.umap(
        adata,
        color="CellType",  # 使用 Patients 列进行着色
        title="UMAP Clustering by Patients",
        legend_loc="on data",
        frameon=True,
        show=False
    )
    plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
    fig = plt.gcf()
    fig.set_size_inches(10, 8)  # 设置更大的画布尺寸
    
    # 添加轴标签
    ax = plt.gca()
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    legend = ax.legend(fontsize=12, title_fontsize=14)
    pdf.savefig(fig)
    plt.close()


with PdfPages(pdf_filename) as pdf:
    sc.pl.umap(
        adata,
        color="CellType",  # 使用 CellType 列进行着色
        title="UMAP Clustering by Patients",
        legend_loc="on data",
        frameon=True,
        show=False
    )
    
    # 设置画布尺寸
    plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
    fig = plt.gcf()
    fig.set_size_inches(10, 8)  # 设置画布尺寸为 10x8 英寸

    # 获取当前轴对象
    ax = plt.gca()
    ax.set_xlabel("UMAP 1", fontsize=14, fontweight='bold')  # 设置 X 轴标签字体加粗
    ax.set_ylabel("UMAP 2", fontsize=14, fontweight='bold')  # 设置 Y 轴标签字体加粗

    legend = ax.legend(fontsize=12, title_fontsize=14)  # 控制图例的字体大小和标题大小
    # 保存图像到 PDF
    pdf.savefig(fig)
    plt.close()