library(Seurat)
library(dplyr)
library(patchwork)

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
#h1b1_obj <- subset(h1b1_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
#h1b2_obj <- subset(h1b2_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
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

# normalize - log normalize
bm_obj <- NormalizeData(bm_obj)
pbmc_obj <- NormalizeData(pbmc_obj)

# Find variable features
bm_obj <- FindVariableFeatures(bm_obj, selection.method = "vst", nfeatures = 2000)
pbmc_obj <- FindVariableFeatures(pbmc_obj, selection.method = "vst", nfeatures = 2000)

# Scale data - default scales on only variable features and not all genes
bm_obj <- ScaleData(bm_obj)
pbmc_obj <- ScaleData(pbmc_obj)

# Do PCA for the set
bm_obj <- RunPCA(bm_obj)
pbmc_obj <- RunPCA(pbmc_obj)

# Do the jackstraw to check for dimensions to include
# bm_obj <- JackStraw(bm_obj, num.replicate = 100)
# bm_obj <- ScoreJackStraw(bm_obj, dims = 1:20)
# 
# pbmc_obj <- JackStraw(pbmc_obj, num.replicate = 100)
# pbmc_obj <- ScoreJackStraw(pbmc_obj, dims = 1:20)

# Plot the elbow plot to check the dims
# ep1 <- ElbowPlot(bm_obj)
# ep2 <- ElbowPlot(pbmc_obj)

# Find neighbors and clusters
bm_obj <- FindNeighbors(bm_obj, dims = 1:12)
bm_obj <- FindClusters(bm_obj, resolution = 0.5)

pbmc_obj <- FindNeighbors(pbmc_obj, dims = 1:12)
pbmc_obj <- FindClusters(pbmc_obj, resolution = 0.5)

# Run UMAP/tSNE
bm_obj <- RunTSNE(bm_obj, dims = 1:12)
pbmc_obj <- RunTSNE(pbmc_obj, dims = 1:12)

# get the final cluster plots
cp1 <- DimPlot(bm_obj, reduction = "tsne")
cp2 <- DimPlot(pbmc_obj, reduction = "tsne")

plot(cp1 + cp2)

# save the rds
saveRDS(bm_obj, "C:\\Users\\scsac\\Desktop\\GATech\\GhosnLab\\SachinKumar_GhosnLab Shared folder\\DifferentTissues_sameBatch\\BM.rds")
saveRDS(pbmc_obj, "C:\\Users\\scsac\\Desktop\\GATech\\GhosnLab\\SachinKumar_GhosnLab Shared folder\\DifferentTissues_sameBatch\\PBMC.rds")