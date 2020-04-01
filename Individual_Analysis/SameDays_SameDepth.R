library(Seurat)
library(dplyr)
library(patchwork)

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
pbmc1_obj <- ThreshFilter(pbmc1_obj)
pbmc2_obj <- ThreshFilter(pbmc2_obj)
pbmc3_obj <- ThreshFilter(pbmc3_obj)

# normalize - log normalize
pbmc1_obj <- NormalizeData(pbmc1_obj)
pbmc2_obj <- NormalizeData(pbmc2_obj)
pbmc3_obj <- NormalizeData(pbmc3_obj)

# Find variable features
pbmc1_obj <- FindVariableFeatures(pbmc1_obj, selection.method = "vst", nfeatures = 2000)
pbmc2_obj <- FindVariableFeatures(pbmc2_obj, selection.method = "vst", nfeatures = 2000)
pbmc3_obj <- FindVariableFeatures(pbmc3_obj, selection.method = "vst", nfeatures = 2000)

# Scale data - default scales on only variable features and not all genes
pbmc1_obj <- ScaleData(pbmc1_obj)
pbmc2_obj <- ScaleData(pbmc2_obj)
pbmc3_obj <- ScaleData(pbmc3_obj)

# Do PCA for the set
pbmc1_obj <- RunPCA(pbmc1_obj)
pbmc2_obj <- RunPCA(pbmc2_obj)
pbmc3_obj <- RunPCA(pbmc3_obj)

# Do the jackstraw to check for dimensions to include
# pbmc1_obj <- JackStraw(pbmc1_obj, num.replicate = 100)
# pbmc1_obj <- ScoreJackStraw(pbmc1_obj, dims = 1:20)
# 
# pbmc2_obj <- JackStraw(pbmc2_obj, num.replicate = 100)
# pbmc2_obj <- ScoreJackStraw(pbmc2_obj, dims = 1:20)
# 
# pbmc3_obj <- JackStraw(pbmc3_obj, num.replicate = 100)
# pbmc3_obj <- ScoreJackStraw(pbmc3_obj, dims = 1:20)

# Plot the elbow plot to check the dims
# ep1 <- ElbowPlot(pbmc1_obj)
# ep2 <- ElbowPlot(pbmc2_obj)
# ep3 <- ElbowPlot(pbmc3_obj)

# Find neighbors and clusters
pbmc1_obj <- FindNeighbors(pbmc1_obj, dims = 1:12)
pbmc1_obj <- FindClusters(pbmc1_obj, resolution = 0.5)

pbmc2_obj <- FindNeighbors(pbmc2_obj, dims = 1:12)
pbmc2_obj <- FindClusters(pbmc2_obj, resolution = 0.5)

pbmc3_obj <- FindNeighbors(pbmc3_obj, dims = 1:12)
pbmc3_obj <- FindClusters(pbmc3_obj, resolution = 0.5)

# Run UMAP/tSNE
pbmc1_obj <- RunTSNE(pbmc1_obj, dims = 1:12)
pbmc2_obj <- RunTSNE(pbmc2_obj, dims = 1:12)
pbmc3_obj <- RunTSNE(pbmc3_obj, dims = 1:12)

# get the final cluster plots
cp1 <- DimPlot(pbmc1_obj, reduction = "tsne")
cp2 <- DimPlot(pbmc2_obj, reduction = "tsne")
cp3 <- DimPlot(pbmc3_obj, reduction = "tsne")

plot(cp1 + cp2 + cp3)

# save the rds
saveRDS(pbmc1_obj, "C:\\Users\\scsac\\Desktop\\GATech\\GhosnLab\\SachinKumar_GhosnLab Shared folder\\SameTissues_SameBatch\\PBMC1.rds")
saveRDS(pbmc2_obj, "C:\\Users\\scsac\\Desktop\\GATech\\GhosnLab\\SachinKumar_GhosnLab Shared folder\\SameTissues_SameBatch\\PBMC2.rds")
saveRDS(pbmc3_obj, "C:\\Users\\scsac\\Desktop\\GATech\\GhosnLab\\SachinKumar_GhosnLab Shared folder\\SameTissues_SameBatch\\PBMC3.rds")