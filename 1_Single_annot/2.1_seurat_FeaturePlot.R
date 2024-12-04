options(stringsAsFactors = FALSE)
library(Seurat)
library(dplyr)
library(DoubletFinder)
library(limma)
library(Cairo)
library(harmony)
library(car)
library(SingleCellExperiment)
library(ggplot2)
library(data.table)
options(bitmapType = "cairo")

#
source("/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/code/Single_cell/codes/scRNA_primary.R")
source("/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/code/Single_cell/codes/colors.R")
source("/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/code/Single_cell/codes/signatureGenes.R") #human signature gene 
#source("/jdfsbjcas1/ST_BJ/P21H28400N0232/wenfeng/project/codes/signatureGene_mouse.R")
source("/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/code/Single_cell/codes/0_initNEWanalysis.R")
source("/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/code/Single_cell/codes/DEG_wilcox.R")


#Marker Gene
SigGeneral_all_1 <- list(immuneCells= c("PTPRC"),
                       macrophage = c("ADGRE1","CD163", "SLC11A1", "APOC1", "CD86", "CSF1R", "SLCO2B1", "CD68", "F13A1", "CD14", "AIF1", "CD80","FCER1G", "FCGR3A", "TYROBP"),
                       fibroblast = c("ACTA2", "SULF1", "CTGF", "TAGLN","TPM1", "GINS1" ,"THY1", "RBP1", "COL1A2", "COL1A1", "C1R", "IGFBP7", "SFRP2", "MGP",
                                      "C1S", "DCN", "CXCL14", "COL3A1","LUM"),
                      #  Tcell = c("CD2", "CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD8B", "GZMK"),
                      Tcell = c("CD2","CD14","CD163","CD68", "CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD8B", "GZMK"),
                       Bcell = c("KIF4A","FANCI","CHAF1A", "GTSE1", "ASPM", "SPC25","NCAPG2","POLA2","NCAPD3", "CD19", "CD21", "MS4A1", "CD79A", "CD79B", "BLNK"),
                       DCs = c("HLA-DRB1","ITGAX","CD83", "HLA-DMA","HLA-DQB1", "HLA-DPA1","HLA-DPB1","AIF1","LST1","FTL","HLA-C"),
                       Epithelial = c("EPCAM", "KRT3", "KRT14", "MUC1", "TP63","CDH1"),
                       Epithe_Lum = c("HIF1A","CEBPD","BTG1","KRT8","CD9","AQP3","MUC1","KRT18","SLPI","AGR2","LCN2","CLDN4","ANXA1","CD74"),
                       Epithe_Mam = c("PRLR","CLDN4","CSN3","CSN1S1","KRT19","KRT7","KRT8","KRT18"),
                       Adipocytes = c("ADIPOQ","PNPLA2","CFD","CIDEC","APOC1","LRP1","	CIDEA","PLIN1","PCK1","CDO1", "STAT6", "TFE3"),
                       PVL = c("MCAM", "CD146", "ACTA2", "PDGFRB", "LYVE1","RGS5","COL4A1","CALD1","SPARCL1","SPARC","NID1"),
                       Pericyte = c("ACTA2","PDGFRB","MCAM","HIGD1B","ANGPT2","VIM","MFGE8","MYO1B","NOTCH3","COX4I2"),
                       endothelial = c("PECAM1", "CD31", "CD34", "HSPG2", "LDB2", "GPR116", "PTPRB", "VWF", "DOCK9", "CDH5", "SELE"),
                       endocrine = c("CHGB", "CHGA", "TTR", "SCG5", "SLC30A8", "GCG", "CLU","CPE","SCG3","CRYBA2","TM4SF4","SCGN"),
                       tuftCells = c("AZGP1", "PLCG2", "HPGDS", "AVIL", "PAEP", "SH2D6", "BMX", "LRMP"),
                       neutrophils = c("A1BG","ALOX5","ASAH1","CD33","CD44","CD63","CTSG","DOCK2","HSPA1B","HSP90AA1"),
                       mastCells = c("RHOH","BTK","FER","GATA2","IL4R","KIT","LCP2","ENPP3","RAC2","LAT2","CD84","LAT","ADGRE2","UNC13D","NDEL1", 
                                     "CPA3","TPSB2","TPSAB1","MS4A2","SLC18A2","IL1RL1","PTGS1","HPGDS"),
                       NKcell = c("GNLY", "FCER1G", "KLRB1", "KLRC1", "AREG", "XCL1"),
                       NKT = c("NKG7", "GNLY", "GZMA", "GZMB", "FCGR3A", "KLRB1"),
                       PAADcancer = c("TSPAN8","AGR2","KRT19","CEACAM6","MMP7","TFF2","CLDN18","TFF1","S100P","SPINK1", "KRT7", "KRT18", "MUC1", "FXYD3"),
                       ducat1Cell = c("AMBP", "FXYD2"),
                       muscle = c("MYH11","NR2F2","CRYAB","TPM2","GUCA2B","GUCA2A","LMOD1","TPPP3","TAGLN","ACTA2"),
                       plasmocyte = c("IGLC2","IGLC3","IGHG3","IGHG1","IGHG4","JCHAIN","IGLC1"),
                       DuctalCell = c("CAPS","TPPP3","MIA","RSPH1","PIFO","LCN2","GDF15","AGR3","CETN2"),
                       Cyclegene = c("MKI67","CDC20", "CENPF","PTTG1","TOP2A","CCNB1","PCNA"),
                       astrocyte = c("GFAP", "BMPR1B", "CD44", "SLC1A2", "AQP4", "S100B", "GJB7", "ALDH1L1", "ALDOC", "MLC1"),
                       oligodendrocyte = c("MBP", "SOX10", "MOG", "CA2", "CNP", "RTN4", "PLP1", "PLP2", "OPALIN", "OMG","OLIG1", "TNR", "ALCAM", "PLLP"),
                       neuron = c("DCX", "NCAM1", "NCAM2", "NEUROD1", "NEUROD2", "ENO2", "RBFOX3", "RBFOX2", "RBFOX1","MAP2", "TUBB3", "NEFM", "NEFH", 
                                  "GAP43", "TH", "GAD1", "GAD2", "SLC17A7", "SLC17A6","FEV", "PET1", "SLC6A4", "DLG4", "SYP", "BSN"),
                       neuronStem = c("POU3F2","NES", "FABP7","FGFR1","FOXA2", "NEFM","NEFL","NR4A2","PAX3", "PAX6", "TUBB3","SOX21", "EOMES","SOX2",
                                      "CD133","PROML1","SLC1A3", "FABP7")
)
Immune_cell <- list(immuneCells= c("PTPRC"),
                      macrophage = c("CD163","ADGRE1","CD163", "SLC11A1", "APOC1", "CD86", "CSF1R", "SLCO2B1", "CD68", "F13A1", "CD14", "AIF1", "CD80","FCER1G", "FCGR3A", "TYROBP"),
                       Tcell = c("CD3E","CD2", "CD3D",  "CD3G","CD4", "CD8A", "CD8B", "GZMK", "FOXP3","TRDC","NKG7","CD79A"),
                      NK = c("KLRB1","KLRC1", "GNLY", "XCL1", "AREG", "FCER1G","FCGR3A","KLRF1", "PRF1"),
                      NKT = c("NKG7","FCER1G","FCGR3A","GZMA","KLRB1","KLRC1","GZMB", "NCAM1"), 
                      ILC = c("SOX4", "KLRB1","BCL6","LST1","RORC"),
                       Bcell = c("KIF4A","FANCI","CHAF1A", "GTSE1", "ASPM", "SPC25","NCAPG2","POLA2","NCAPD3", "CD19", "MS4A1", "CD79A", "CD79B", "BLNK"), #"CD21", 
                       plasmocyte = c("IGLC2","IGLC3","IGHG3","IGHG1","IGHG4","JCHAIN","IGLC1"),
                       DCs = c("HLA-DRB1","ITGAX","CD83", "HLA-DMA","HLA-DQB1", "HLA-DPA1","HLA-DPB1","AIF1","LST1","FTL","HLA-C"),
                       NKcell = c("GNLY", "FCER1G", "KLRB1", "KLRC1", "AREG", "XCL1"),
                       mastCells = c("RHOH","BTK","FER","GATA2","IL4R","KIT","LCP2","ENPP3","RAC2","LAT2","CD84","LAT","ADGRE2","UNC13D","NDEL1", 
                                     "CPA3","TPSB2","TPSAB1","MS4A2","SLC18A2","IL1RL1","PTGS1","HPGDS"),
                       neutrophils = c("CSF3R","S100A8", "S100A9", "LTF", 'CEBPE', 'CEBPA', 'ETV6', 'FOXP1', 'GFI1', 'SPI1', 'STAT3', 'ELANE', 'GSTM1', 'LCN2',  'MPO', 'PRTN3') #'LYZ2',
)
T_marker <- list(CD4NV_CM_rest = c("LEF1","ATM","SELL","KLF2","ITGA6"),
                  CD4CD8_rest = c("IL7R","CD52","S100A4","TGFB3","AQP3","NLRP3","KLF2","ITGB7"),
                  # Treg = c("FOXP3","CTLA4","MTNFRSF4","IRF4", "BATF", "TNFRSF18","TOX2","PRDMI"),
                  T_IFNresp = c("IFIT3", "IFIT2","STAT1", "MX1","IRF7", "ISG15","IFITM3","OAS2","JAK2","SOCS1","TRIM2"),
                  T_prof = c("LIF", "IL2","CENPV","NME1","FABP5","ORC6","G0S2","GCK"),
                  T_cytotoxic = c("CCL5","GZMK","GNLY","EOMES","ZNF683","KLRG1","NKG7","ZEB2"),
                #   T_cytokine = c("CCL3", "IFNG","CCL4", "XCL1","XCL2", "CSF2","IL10","HOPX","TIM3", "LAG3","PRF1","TNFRSF9","NKG7","IL26"),
                  CD4conv = c("LTB", "SPOCK2","SOCS3","PBXIP1"),
                #   Treg = c("FOXP3","TNFRSF4","CARD16","IL2RA","TBC1D4","CD25"), ## IL2RA CD25
                  CD8 = c("CCL4", "CST7", "GZMA", "GZMH"),
                  Texh = c("CTLA4", "LAG3", "PTMS", "PDCD1", "TIGIT"),
                  Tmem = c("CD44", "GZMK","IL7R", "CD69", "CD27"),
                  Tresid = c("SELL","ITGAE"),
                  Teff = c("CTSW","PRF1","GZMB","FGFBP2", "FCGR3A","TRGC2","GNLY","TYROBP", "IFNG", "CXCL13"),
                  Tgd = c("TYROBP","FCER1G","CMC1","KLRD1","KLRG1","KLRF1","TRDC","CTSW","PRF1","GZMB","FGFBP2", "FCGR3A","GNLY"))

Fibroblast_cell <- list(fibroblast = c("PDGFRA", "ACTA2", "SULF1", "TAGLN","TPM1", "GINS1" ,"THY1", "RBP1", "COL1A2", "COL1A1", "C1R", "IGFBP7", "SFRP2", "MGP","C1S", "DCN", "CXCL14", "COL3A1","LUM"), #"CTGF",
                       endothelial = c("VWF","PECAM1", "CD34", "HSPG2", "LDB2",  "PTPRB", "DOCK9", "CDH5", "SELE"), #"CD31","GPR116",
                       muscle = c("ACTA2","MYH11","NR2F2","CRYAB","TPM2","GUCA2B","GUCA2A","LMOD1","TPPP3","TAGLN"),
                       MSC = c("CD44", "ITGA1","NT5E","THY1"),
                       mycaf = c("ACTA2","VIM","CCN2", "COL1A1","COL5A1","COL6A1","TNC","TGFB1","THY1","TAGLN","COL12A1","PDGFRB"), #TGFβ/SMAD2/3 
                       icaf = c("IL1A","IL1B","IL6","IL11","LIF","CLEC3B","COL14A1","GSN","CXCL12","CXCL14","PDGFRA"), ##IL-1/JAK-STAT3,"LY6C1"
                       #apcaf = c("SLPI","CD74","NKAIN4", "IRF5"),#"SAA3","H2-Ab1",
                       psc = c("DES", "GFAP", "CHRNA1"),
                       PVL = c("MCAM",  "ACTA2", "PDGFRB", "LYVE1","RGS5","COL4A1","CALD1","SPARCL1","SPARC","NID1"), #"CD146",
                       Pericyte = c("KCNJ8","ACTA2","PDGFRB","MCAM","HIGD1B","ANGPT2","VIM","MFGE8","MYO1B","NOTCH3","COX4I2"),
                       Cyclegene = c("MKI67","CDC20", "CENPF","PTTG1","TOP2A","CCNB1","PCNA")
)
Epithe_cell <- list(Epithelial = c("EPCAM", "KRT3", "KRT14", "MUC1", "TP63","CDH1"),
                       Epithe_Lum = c("HIF1A","CEBPD","BTG1","KRT8","CD9","AQP3","MUC1","KRT18","SLPI","AGR2","LCN2","CLDN4","ANXA1","CD74"),
                       Epithe_Mam = c("PRLR","CLDN4","CSN3","CSN1S1","KRT19","KRT7","KRT8","KRT18"),
                       Adipocytes = c("ADIPOQ","PNPLA2","CFD","CIDEC","APOC1","LRP1","PLIN1","PCK1","CDO1", "STAT6", "TFE3"),#"CIDEA"
                       tuftCells = c("AZGP1", "PLCG2", "HPGDS", "AVIL", "PAEP", "SH2D6", "BMX", "LRMP"),
                       PAADcancer = c("TSPAN8","AGR2","KRT19","CEACAM6","MMP7","TFF2","CLDN18","TFF1","S100P","SPINK1", "KRT7", "KRT18", "MUC1", "FXYD3"),
                       Cyclegene = c("MKI67","CDC20", "CENPF","PTTG1","TOP2A","CCNB1","PCNA")
)

PC_marker <- list(immuneCells= c("PTPRC"),
                  macrophage = c("CD163","ADGRE1","CD163", "SLC11A1", "APOC1", "CD86", "CSF1R", "SLCO2B1", "CD68", "F13A1", "CD14", "AIF1", "CD80","FCER1G", "FCGR3A", "TYROBP"),
                  fibroblast = c("PDGFRA","ACTA2", "SULF1", "CTGF", "TAGLN","TPM1", "GINS1" ,"THY1", "RBP1", "COL1A2", "COL1A1", "C1R", "IGFBP7", "SFRP2", "MGP","C1S", "DCN", "CXCL14", "COL3A1","LUM"),
                  Tcell = c("CD3E","CD2", "CD3D","CD3G", "CD4", "CD8A", "CD8B", "GZMK"),
                  Bcell = c("KIF4A","FANCI","CHAF1A", "GTSE1", "ASPM", "SPC25","NCAPG2","POLA2","NCAPD3", "CD19", "CD21", "MS4A1", "CD79A", "CD79B", "BLNK"),
                  DCs = c("HLA-DRB1","ITGAX","CD83", "HLA-DMA","HLA-DQB1", "HLA-DPA1","HLA-DPB1","AIF1","LST1","FTL","HLA-C"),
                  Epithelial = c("EPCAM", "KRT3", "KRT14", "MUC1", "TP63","CDH1"),
                  Epithe_Lum = c("HIF1A","CEBPD","BTG1","KRT8","CD9","AQP3","MUC1","KRT18","SLPI","AGR2","LCN2","CLDN4","ANXA1","CD74"),
                  Epithe_Mam = c("PRLR","CLDN4","CSN3","CSN1S1","KRT19","KRT7","KRT8","KRT18"),
                  PVL = c("MCAM", "CD146", "ACTA2", "PDGFRB", "LYVE1","RGS5","COL4A1","CALD1","SPARCL1","SPARC","NID1"),
                  Pericyte = c("KCNJ8","ACTA2","PDGFRB","MCAM","HIGD1B","ANGPT2","VIM","MFGE8","MYO1B","NOTCH3","COX4I2"),
                  endothelial = c("PECAM1", "CD31", "CD34", "HSPG2", "LDB2", "GPR116", "PTPRB", "VWF", "DOCK9", "CDH5", "SELE"),
                  tuftCells = c("AZGP1", "PLCG2", "HPGDS", "AVIL", "PAEP", "SH2D6", "BMX", "LRMP"),
                  neutrophils = c("A1BG","ALOX5","ASAH1","CD33","CD44","CD63","CTSG","DOCK2","HSPA1B","HSP90AA1"),
                  mastCells = c("RHOH","BTK","FER","GATA2","IL4R","KIT","LCP2","ENPP3","RAC2","LAT2","CD84","LAT","ADGRE2","UNC13D","NDEL1", 
                                "CPA3","TPSB2","TPSAB1","MS4A2","SLC18A2","IL1RL1","PTGS1","HPGDS"),
                  NKT = c("NKG7", "GNLY", "GZMA", "GZMB", "FCGR3A", "KLRB1"),
                  muscle = c("ACTA2","MYH11","NR2F2","CRYAB","LMOD1","TPPP3","TAGLN"),
                  Cyclegene = c("MKI67","CDC20", "CENPF","PTTG1","TOP2A","CCNB1","PCNA"),
                  plasmocyte = c("IGLC2","IGLC3","IGHG3","IGHG1","IGHG4","JCHAIN","IGLC1"),
                  TSK = c("MMP10","INHBA","PTHLH","MAGEA4","FEZ1","NT5E","IL24","LAMC2","KCNMA","SLITRK6")
)

PC <- list(immuneCells= c("PTPRC"),
                      macrophage = c("CD163","ADGRE1","CD163", "SLC11A1", "APOC1", "CD86", "CSF1R", "SLCO2B1", "CD68", "F13A1", "CD14", "AIF1", "CD80","FCER1G", "FCGR3A", "TYROBP"),
                       Tcell = c("CD3E","CD2", "CD3D",  "CD3G","CD4", "CD8A", "CD8B", "GZMK", "FOXP3","TRDC","NKG7","CD79A"),
                       Bcell = c("KIF4A","FANCI","CHAF1A", "GTSE1", "ASPM", "SPC25","NCAPG2","POLA2","NCAPD3", "CD19", "MS4A1", "CD79A", "CD79B", "BLNK"), #"CD21", 
                       plasmocyte = c("IGLC2","IGLC3","IGHG3","IGHG1","IGHG4","JCHAIN","IGLC1"),
                       DCs = c("HLA-DRB1","ITGAX","CD83", "HLA-DMA","HLA-DQB1", "HLA-DPA1","HLA-DPB1","AIF1","LST1","FTL","HLA-C"),
                       NKcell = c("GNLY", "FCER1G", "KLRB1", "KLRC1", "AREG", "XCL1"),
                       mastCells = c("RHOH","BTK","FER","GATA2","IL4R","KIT","LCP2","ENPP3","RAC2","LAT2","CD84","LAT","ADGRE2","UNC13D","NDEL1", 
                                     "CPA3","TPSB2","TPSAB1","MS4A2","SLC18A2","IL1RL1","PTGS1","HPGDS"),
                       neutrophils = c("CSF3R","S100A8", "S100A9", "LTF", 'CEBPE', 'CEBPA', 'ETV6', 'FOXP1', 'GFI1', 'SPI1', 'STAT3', 'ELANE', 'GSTM1', 'LCN2',  'MPO', 'PRTN3'),
                       fibroblast = c("PDGFRA", "ACTA2", "SULF1", "TAGLN","TPM1", "GINS1" ,"THY1", "RBP1", "COL1A2", "COL1A1", "C1R", "IGFBP7", "SFRP2", "MGP","C1S", "DCN", "CXCL14", "COL3A1","LUM"), #"CTGF",
                       endothelial = c("VWF","PECAM1", "CD34", "HSPG2", "LDB2",  "PTPRB", "DOCK9", "CDH5", "SELE"), #"CD31","GPR116",
                       muscle = c("ACTA2","MYH11","NR2F2","CRYAB","TPM2","GUCA2B","GUCA2A","LMOD1","TPPP3","TAGLN"),
                       Pericyte = c("KCNJ8","ACTA2","PDGFRB","MCAM","HIGD1B","ANGPT2","VIM","MFGE8","MYO1B","NOTCH3","COX4I2"),
                       Cyclegene = c("MKI67","CDC20", "CENPF","PTTG1","TOP2A","CCNB1","PCNA"),
                       Epithelial = c("EPCAM", "KRT3", "KRT14", "MUC1", "TP63","CDH1"),
                       Epithe_Lum = c("HIF1A","CEBPD","BTG1","KRT8","CD9","AQP3","MUC1","KRT18","SLPI","AGR2","LCN2","CLDN4","ANXA1","CD74"),
                       Epithe_Mam = c("PRLR","CLDN4","CSN3","CSN1S1","KRT19","KRT7","KRT8","KRT18"),
                       Adipocytes = c("ADIPOQ","PNPLA2","CFD","CIDEC","APOC1","LRP1","PLIN1","PCK1","CDO1", "STAT6", "TFE3"),#"CIDEA"
                       tuftCells = c("AZGP1", "PLCG2", "HPGDS", "AVIL", "PAEP", "SH2D6", "BMX", "LRMP") #'LYZ2'
                       )
length(names(PC_marker))
#Immune_cell 31   PC_markder 20 Fibroblast_cell 13

load('/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC_1/Data/seurat_merge_Immune.RData')
#Work directory
workdir = '/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC1/FeaturePlot'
if(!dir.exists(paste0(workdir,'/Mesenchymal'))){dir.create(paste0(workdir,'/Mesenchymal'), recursive = TRUE)}
setwd("/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC/Test/sctrans_CCA/CCA/FeaturePlot")
n <- 1


marker <- Immune_cell
marker <- Epithe_cell
marker <- PC
seurat <- seurat_merge
while( as.numeric(n) <= length(names(marker))){
    name <- paste0(names(marker)[n],".png")
    genes <- unlist(marker[n])
    genes_count <- length(genes)
    png(name, width = genes_count * 100, height = genes_count * 20, units = "px", res = 300)
    combined_plot <- ggplot() + theme_void()
    for (gene in genes) {
        p <- FeaturePlot(object = seurat, features = gene,raster = TRUE)
        # p <- p + geom_text(aes(x = 1, y = 1,label = ""), vjust = -0.5, size = 3, nudge_y = 0.1)
        combined_plot <- combined_plot + p 
    }
    combined_plot <- combined_plot + facet_wrap(~ genes, scales = "free", nrow  = 5) #combined_plot +
    ggsave(name, plot = combined_plot, width = 20, height = 14)
    dev.off()
    n <- n + 1
}
# while( as.numeric(n) <= 11 ){
#   a <- unlist(Fibroblast_cell[n])
#   # a <- PC_marker[n]
#   b <- FeaturePlot(seurat_merge2, features =a)
#   c <- paste0(names(Fibroblast_cell)[n],".png")
#   ggsave(b,filename = c,width = 25,height = 25,device = "png")
#   n= n+1
# }
###################下面没用
setwd('/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC/FeaturePlot')
load('/jdfsbjcas1/ST_BJ/P21H28400N0232/fanfengpu/PC/Data/seurat_harmony/1.1.4_seurat_merge_harmony.RData')
Immune_cell <- list(immuneCells= c("KRT14","KRT16","KRT17","KRT5","KRT6A","KRT6B","KRT6C","KRT7","KRTDAP"))
n <- 2

