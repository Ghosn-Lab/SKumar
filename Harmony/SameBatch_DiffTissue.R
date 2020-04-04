library(harmony)
library(Seurat)
library(sctransform)
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

bm_pbmc <- merge(bm_obj, pbmc_obj, add.cell.ids = c("BM", "PBMC"), project = "Mixed_Tissue") # make a set for running future functions

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

dT <- ThreshFilter(bm_pbmc)

dT <- Seurat::NormalizeData(dT, verbose = FALSE)
dT <- FindVariableFeatures(dT, selection.method = "vst", nfeatures = 2000) 
dT <- ScaleData(dT, verbose = FALSE)
dT <- RunPCA(dT, npcs = 20, verbose = FALSE)

options(repr.plot.height = 2.5, repr.plot.width = 6)
dT <- dT %>% 
  RunHarmony("orig.ident", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(dT, 'harmony')

dT <- dT %>% 
  RunTSNE(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()
options(repr.plot.height = 4, repr.plot.width = 10)
DimPlot(dT, reduction = "tsne", pt.size = .1, split.by = 'orig.ident')
# save the rds
saveRDS(dT, "C:\\Users\\scsac\\Desktop\\GATech\\GhosnLab\\SachinKumar_GhosnLab Shared folder\\DifferentTissues_sameBatch\\BM_PBMC_Harmony_LogN.rds")