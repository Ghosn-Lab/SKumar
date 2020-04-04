library(harmony)
library(Seurat)
library(sctransform)
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

h1 <- merge(h1b1_obj, h1b2_obj, add.cell.ids = c("H1B1", "H1B2"), project = "H1") # make a set for running future functions

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
dSD <- ThreshFilter(h1)

dSD <- Seurat::NormalizeData(dSD, verbose = FALSE)
dSD <- FindVariableFeatures(dSD, selection.method = "vst", nfeatures = 2000) 
dSD <- ScaleData(dSD, verbose = FALSE)
dSD <- RunPCA(dSD, npcs = 20, verbose = FALSE)

options(repr.plot.height = 2.5, repr.plot.width = 6)
dSD <- dSD %>% 
  RunHarmony("orig.ident", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(dSD, 'harmony')

dSD <- dSD %>% 
  RunTSNE(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()
options(repr.plot.height = 4, repr.plot.width = 10)
DimPlot(dSD, reduction = "tsne", pt.size = .1, split.by = 'orig.ident')
# save the rds
saveRDS(dSD, "C:\\Users\\scsac\\Desktop\\GATech\\GhosnLab\\SachinKumar_GhosnLab Shared folder\\Sample_5\\H1_Harmony_LogN.rds")