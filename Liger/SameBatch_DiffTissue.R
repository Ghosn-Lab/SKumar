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
bm.data <- Read10X(data.dir = "C:\\Users\\scsac\\Desktop\\GATech\\GhosnLab\\SachinKumar_GhosnLab Shared folder\\DifferentTissues_sameBatch\\BM B cells\\batch\\")
bm_obj <- CreateSeuratObject(counts = bm.data$`Gene Expression`, project = "BM")
bm_obj[['ADT']] <- CreateAssayObject(counts = bm.data$`Antibody Capture`)

pbmc.data <- Read10X(data.dir = "C:\\Users\\scsac\\Desktop\\GATech\\GhosnLab\\SachinKumar_GhosnLab Shared folder\\DifferentTissues_sameBatch\\PBMC B cells\\batch\\")
pbmc_obj <- CreateSeuratObject(counts = pbmc.data$`Gene Expression`, project = "PBMC")
pbmc_obj[['ADT']] <- CreateAssayObject(counts = pbmc.data$`Antibody Capture`)

# Add percentage mito as a column
bm_obj[["percent.mt"]] <- PercentageFeatureSet(bm_obj, pattern = "^MT-")
pbmc_obj[["percent.mt"]] <- PercentageFeatureSet(pbmc_obj, pattern = "^MT-")

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

bm_obj <- ThreshFilter(bm_obj)
pbmc_obj <- ThreshFilter(pbmc_obj)

bm_obj <- RenameCells(object = bm_obj, add.cell.id = "BM")
pbmc_obj <- RenameCells(object = pbmc_obj, add.cell.id = "PBMC")

lObj.diffTissue <- seuratToLiger(list(bm_obj, pbmc_obj), names = c('BM', 'PBMC'), num.hvg.info = 2000)
lObj.diffTissue <- normalize(lObj.diffTissue)
lObj.diffTissue <- selectGenes(lObj.diffTissue, var.thresh = c(0.3, 0.875), do.plot = F)
lObj.diffTissue <- scaleNotCenter(lObj.diffTissue)

lObj.diffTissue <- optimizeALS(lObj.diffTissue, k=20, thresh = 5e-5, nrep = 3)
lObj.diffTissue <- quantile_norm(lObj.diffTissue)
lObj.diffTissue <- runTSNE(lObj.diffTissue, dims = 1:15)

p_a <- plotByDatasetAndCluster(lObj.diffTissue, return.plots = T) 
# Modify plot output slightly
print(p_a[[1]])
print(p_a[[2]])
# save the rds
saveRDS(lObj.diffTissue, "C:\\Users\\scsac\\Desktop\\GATech\\GhosnLab\\SachinKumar_GhosnLab Shared folder\\DifferentTissues_sameBatch\\BM_PBMC_Liger.rds")