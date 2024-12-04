#### to creat the dirs for a new analysis object, such as a new sample, or specific cell type 

createObjdir <- function(workdir = workdir){
  if(!dir.exists(paste0(workdir, "saveData/"))){dir.create(paste0(workdir, "/saveData/"), recursive = TRUE)}
  if(!dir.exists(paste0(workdir, "codes/"))){dir.create(paste0(workdir, "codes/"), recursive = TRUE)}
  if(!dir.exists(paste0(workdir, "CNV/infercnv/"))){dir.create(paste0(workdir, "CNV/infercnv/"), recursive = TRUE)}
  if(!dir.exists(paste0(workdir, "Results/Figures/"))){dir.create(paste0(workdir, "Results/Figures/"), recursive = TRUE)}
  if(!dir.exists(paste0(workdir, "pathwayEnrich/metab/"))){dir.create(paste0(workdir, "pathwayEnrich/metab/"), recursive = TRUE)}
  if(!dir.exists(paste0(workdir, "pathwayEnrich/pathway/"))){dir.create(paste0(workdir, "pathwayEnrich/pathway/"), recursive = TRUE)}
  if(!dir.exists(paste0(workdir, "NMF/cNMF/"))){dir.create(paste0(workdir, "NMF/cNMF/"), recursive = TRUE)}
}
