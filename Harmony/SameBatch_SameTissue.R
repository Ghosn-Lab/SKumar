library(harmony)
library(Seurat)
library(sctransform)
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
pbmc <- merge(pbmc_list[[1]], pbmc_list[2:length(pbmc_list)], add.cell.ids = c("PBMC1", "PBMC2", "PBMC3"), project = "PBMC") # make a set for running future functions

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

sTPbmc <- ThreshFilter(pbmc)

sTPbmc <- Seurat::NormalizeData(sTPbmc, verbose = FALSE)
sTPbmc <- FindVariableFeatures(sTPbmc, selection.method = "vst", nfeatures = 2000) 
sTPbmc <- ScaleData(sTPbmc, verbose = FALSE)
sTPbmc <- RunPCA(sTPbmc, npcs = 20, verbose = FALSE)

options(repr.plot.height = 2.5, repr.plot.width = 6)
sTPbmc <- sTPbmc %>% 
  RunHarmony("orig.ident", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(sTPbmc, 'harmony')

sTPbmc <- sTPbmc %>% 
  RunTSNE(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()
options(repr.plot.height = 4, repr.plot.width = 10)
DimPlot(sTPbmc, reduction = "tsne", pt.size = .1, split.by = 'orig.ident')
# save the rds
saveRDS(sTPbmc, "C:\\Users\\scsac\\Desktop\\GATech\\GhosnLab\\SachinKumar_GhosnLab Shared folder\\SameTissues_SameBatch\\PBMC_Harmony_LogN.rds")