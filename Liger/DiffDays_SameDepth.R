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
h1b1.data <- Read10X(data.dir = "C:\\Users\\scsac\\Desktop\\GATech\\GhosnLab\\SachinKumar_GhosnLab Shared folder\\Sample_5\\5GEX_PBMC3_46yoM_H1_B1\\")
h1b1_obj <- CreateSeuratObject(counts = h1b1.data$`Gene Expression`, project = "H1B1")
h1b1_obj[['ADT']] <- CreateAssayObject(counts = h1b1.data$`Antibody Capture`)

h1b2.data <- Read10X(data.dir = "C:\\Users\\scsac\\Desktop\\GATech\\GhosnLab\\SachinKumar_GhosnLab Shared folder\\Sample_5\\5GEX_PBMC3_46yoM_H1_B2\\")
h1b2_obj <- CreateSeuratObject(counts = h1b2.data$`Gene Expression`, project = "H1B2")
h1b2_obj[['ADT']] <- CreateAssayObject(counts = h1b2.data$`Antibody Capture`)

# Add percentage mito as a column
h1b1_obj[["percent.mt"]] <- PercentageFeatureSet(h1b1_obj, pattern = "^MT-")
h1b2_obj[["percent.mt"]] <- PercentageFeatureSet(h1b2_obj, pattern = "^MT-")

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
h1b1_obj <- ThreshFilter(h1b1_obj)
h1b2_obj <- ThreshFilter(h1b2_obj)

h1b1_obj <- RenameCells(object = h1b1_obj, add.cell.id = "H1B1")
h1b2_obj <- RenameCells(object = h1b2_obj, add.cell.id = "H1B2")

lObj.pbmc <- seuratToLiger(list(h1b1_obj, h1b2_obj), names = c('H1B1', 'H1B2'), num.hvg.info = 2000)
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
saveRDS(lObj.pbmc, "C:\\Users\\scsac\\Desktop\\GATech\\GhosnLab\\SachinKumar_GhosnLab Shared folder\\Sample_5\\H1_Liger.rds")