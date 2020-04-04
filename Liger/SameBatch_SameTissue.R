library(devtools)
library(liger)
library(cowplot)
library(future)
library(Matrix)
library(Seurat)
library(dplyr)
library(patchwork)
options(future.globals.maxSize = 16000 * 1024^2)
memory.limit(size = 12000000000)

# Create the Seurat object
pbmc1.data <- Read10X(data.dir = "C:\\Users\\scsac\\Desktop\\GATech\\GhosnLab\\SachinKumar_GhosnLab Shared folder\\SameTissues_SameBatch\\PBMC_sameBatch\\PBMC 5\\filtered_feature_bc_matrix\\")
pbmc1_obj <- CreateSeuratObject(counts = pbmc1.data, project = "PBMC1")

pbmc2.data <- Read10X(data.dir = "C:\\Users\\scsac\\Desktop\\GATech\\GhosnLab\\SachinKumar_GhosnLab Shared folder\\SameTissues_SameBatch\\PBMC_sameBatch\\PBMC 6\\filtered_feature_bc_matrix\\")
pbmc2_obj <- CreateSeuratObject(counts = pbmc2.data, project = "PBMC2")

pbmc3.data <- Read10X(data.dir = "C:\\Users\\scsac\\Desktop\\GATech\\GhosnLab\\SachinKumar_GhosnLab Shared folder\\SameTissues_SameBatch\\PBMC_sameBatch\\PBMC 7\\filtered_feature_bc_matrix\\")
pbmc3_obj <- CreateSeuratObject(counts = pbmc3.data, project = "PBMC3")

# Add percentage mito as a column
pbmc1_obj[["percent.mt"]] <- PercentageFeatureSet(pbmc1_obj, pattern = "^MT-")
pbmc2_obj[["percent.mt"]] <- PercentageFeatureSet(pbmc2_obj, pattern = "^MT-")
pbmc3_obj[["percent.mt"]] <- PercentageFeatureSet(pbmc3_obj, pattern = "^MT-")

pbmc_list <- c(pbmc1_obj, pbmc2_obj, pbmc3_obj)

# QC for number of features, counts and percent mito
ThreshFilter <- function(obj){
  obj <- subset(obj, subset = nFeature_RNA > 500)
  nUMI <- obj$nCount_RNA
  nGene <- obj$nFeature_RNA
  pMito <- obj$percent.mt
  # <<- operator to make variable available from within subset function 
  # (otherwise scoped only to within ThreshFilter)
  max.UMI <<- 5*sd(nUMI) + median(nUMI)
  max.Gene <<- 5*sd(nGene) + median(nGene)
  max.mito <<- 5*sd(pMito) + median(pMito)
  print(paste0(unique(obj$orig.ident), 
               "; max.umi = ", max.UMI, 
               "; max.Gene = ", max.Gene,
               "; max.mito = ", max.mito))
  # Apply thresholds to subset object
  obj <- obj %>% subset(subset = nCount_RNA < max.UMI) %>%
    subset(nFeature_RNA < max.Gene) %>%
    subset(percent.mt < max.mito)
  return(obj)
}
pbmc_list <- lapply(pbmc_list, function(obj) ThreshFilter(obj))

pbmc_list[[1]] <- RenameCells(object = pbmc_list[[1]], add.cell.id = "PBMC1")
pbmc_list[[2]] <- RenameCells(object = pbmc_list[[2]], add.cell.id = "PBMC2")
pbmc_list[[3]] <- RenameCells(object = pbmc_list[[3]], add.cell.id = "PBMC3")

lObj.pbmc <- seuratToLiger(pbmc_list, names = c('PBMC1', 'PBMC2', 'PBMC3'), num.hvg.info = 2000)
lObj.pbmc <- normalize(lObj.pbmc)
lObj.pbmc <- selectGenes(lObj.pbmc, var.thresh = c(0.3, 0.875), do.plot = F)
lObj.pbmc <- scaleNotCenter(lObj.pbmc)

lObj.pbmc <- optimizeALS(lObj.pbmc, k=20, thresh = 5e-5, nrep = 3)
lObj.pbmc <- quantile_norm(lObj.pbmc)
lObj.pbmc <- runTSNE(lObj.pbmc, dims = 1:15)

p_a <- plotByDatasetAndCluster(lObj.pbmc, return.plots = T) 
# Modify plot output slightly
print(p_a[[1]])
print(p_a[[2]])
# save the rds
saveRDS(lObj.pbmc, "C:\\Users\\scsac\\Desktop\\GATech\\GhosnLab\\SachinKumar_GhosnLab Shared folder\\SameTissues_SameBatch\\PBMC_Liger.rds")