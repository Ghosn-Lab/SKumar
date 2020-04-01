library(Seurat)
library(dplyr)
library(patchwork)

# Create the Seurat object
f1b1.data <- Read10X(data.dir = "C:\\Users\\scsac\\Desktop\\GATech\\GhosnLab\\SachinKumar_GhosnLab Shared folder\\Sample_5\\5GEX_PBMC1_36yoF_F1_B1\\")
f1b1_obj <- CreateSeuratObject(counts = f1b1.data$`Gene Expression`, project = "F1B1")
f1b1_obj[['ADT']] <- CreateAssayObject(counts = f1b1.data$`Antibody Capture`)

f1b2.data <- Read10X(data.dir = "C:\\Users\\scsac\\Desktop\\GATech\\GhosnLab\\SachinKumar_GhosnLab Shared folder\\Sample_5\\5GEX_PBMC1_36yoF_F1_B2\\")
f1b2_obj <- CreateSeuratObject(counts = f1b2.data$`Gene Expression`, project = "F1B2")
f1b2_obj[['ADT']] <- CreateAssayObject(counts = f1b2.data$`Antibody Capture`)

# Add percentage mito as a column
f1b1_obj[["percent.mt"]] <- PercentageFeatureSet(f1b1_obj, pattern = "^MT-")
f1b2_obj[["percent.mt"]] <- PercentageFeatureSet(f1b2_obj, pattern = "^MT-")

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
f1b1_obj <- ThreshFilter(f1b1_obj)
f1b2_obj <- ThreshFilter(f1b2_obj)

# normalize - log normalize
f1b1_obj <- NormalizeData(f1b1_obj)
f1b2_obj <- NormalizeData(f1b2_obj)

# Find variable features
f1b1_obj <- FindVariableFeatures(f1b1_obj, selection.method = "vst", nfeatures = 2000)
f1b2_obj <- FindVariableFeatures(f1b2_obj, selection.method = "vst", nfeatures = 2000)

# Scale data - default scales on only variable features and not all genes
f1b1_obj <- ScaleData(f1b1_obj)
f1b2_obj <- ScaleData(f1b2_obj)

# Do PCA for the set
f1b1_obj <- RunPCA(f1b1_obj)
f1b2_obj <- RunPCA(f1b2_obj)

# Do the jackstraw to check for dimensions to include
# f1b1_obj <- JackStraw(f1b1_obj, num.replicate = 100)
# f1b1_obj <- ScoreJackStraw(f1b1_obj, dims = 1:20)
# 
# f1b2_obj <- JackStraw(f1b2_obj, num.replicate = 100)
# f1b2_obj <- ScoreJackStraw(f1b2_obj, dims = 1:20)

# Plot the elbow plot to check the dims
# ep1 <- ElbowPlot(f1b1_obj)
# ep2 <- ElbowPlot(f1b2_obj)

# Find neighbors and clusters
f1b1_obj <- FindNeighbors(f1b1_obj, dims = 1:12)
f1b1_obj <- FindClusters(f1b1_obj, resolution = 0.5)

f1b2_obj <- FindNeighbors(f1b2_obj, dims = 1:12)
f1b2_obj <- FindClusters(f1b2_obj, resolution = 0.5)

# Run UMAP/tSNE
f1b1_obj <- RunTSNE(f1b1_obj, dims = 1:12)
f1b2_obj <- RunTSNE(f1b2_obj, dims = 1:12)

# get the final cluster plots
cp1 <- DimPlot(f1b1_obj, reduction = "tsne")
cp2 <- DimPlot(f1b2_obj, reduction = "tsne")

plot(cp1 + cp2)

# save the rds
saveRDS(f1b1_obj, "C:\\Users\\scsac\\Desktop\\GATech\\GhosnLab\\SachinKumar_GhosnLab Shared folder\\Sample_5\\F1B1.rds")
saveRDS(f1b2_obj, "C:\\Users\\scsac\\Desktop\\GATech\\GhosnLab\\SachinKumar_GhosnLab Shared folder\\Sample_5\\F1B2.rds")