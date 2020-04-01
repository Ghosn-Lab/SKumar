library(Seurat)
library(dplyr)
library(patchwork)

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
h1b1_obj <- ThreshFilter(h1b1_obj)
h1b2_obj <- ThreshFilter(h1b2_obj)

# normalize - log normalize
h1b1_obj <- NormalizeData(h1b1_obj)
h1b2_obj <- NormalizeData(h1b2_obj)

# Find variable features
h1b1_obj <- FindVariableFeatures(h1b1_obj, selection.method = "vst", nfeatures = 2000)
h1b2_obj <- FindVariableFeatures(h1b2_obj, selection.method = "vst", nfeatures = 2000)

# Scale data - default scales on only variable features and not all genes
h1b1_obj <- ScaleData(h1b1_obj)
h1b2_obj <- ScaleData(h1b2_obj)

# Do PCA for the set
h1b1_obj <- RunPCA(h1b1_obj)
h1b2_obj <- RunPCA(h1b2_obj)

# Do the jackstraw to check for dimensions to include
# h1b1_obj <- JackStraw(h1b1_obj, num.replicate = 100)
# h1b1_obj <- ScoreJackStraw(h1b1_obj, dims = 1:20)
# 
# h1b2_obj <- JackStraw(h1b2_obj, num.replicate = 100)
# h1b2_obj <- ScoreJackStraw(h1b2_obj, dims = 1:20)

# Plot the elbow plot to check the dims
# ep1 <- ElbowPlot(h1b1_obj)
# ep2 <- ElbowPlot(h1b2_obj)

# Find neighbors and clusters
h1b1_obj <- FindNeighbors(h1b1_obj, dims = 1:12)
h1b1_obj <- FindClusters(h1b1_obj, resolution = 0.5)

h1b2_obj <- FindNeighbors(h1b2_obj, dims = 1:12)
h1b2_obj <- FindClusters(h1b2_obj, resolution = 0.5)

# Run UMAP/tSNE
h1b1_obj <- RunTSNE(h1b1_obj, dims = 1:12)
h1b2_obj <- RunTSNE(h1b2_obj, dims = 1:12)

# get the final cluster plots
cp1 <- DimPlot(h1b1_obj, reduction = "tsne")
cp2 <- DimPlot(h1b2_obj, reduction = "tsne")

plot(cp1 + cp2)

# save the rds
saveRDS(h1b1_obj, "C:\\Users\\scsac\\Desktop\\GATech\\GhosnLab\\SachinKumar_GhosnLab Shared folder\\Sample_5\\H1B1.rds")
saveRDS(h1b2_obj, "C:\\Users\\scsac\\Desktop\\GATech\\GhosnLab\\SachinKumar_GhosnLab Shared folder\\Sample_5\\H1B2.rds")