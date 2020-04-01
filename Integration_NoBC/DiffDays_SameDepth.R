library(Seurat)
library(dplyr)
library(patchwork)
options(future.globals.maxSize = 16000 * 1024^2)

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
h1sets <- c(subset(h1, idents = "H1B1"), subset(h1, idents = "H1B2"))
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
h1sets <- lapply(h1sets, function(obj) ThreshFilter(obj))
# do sctransform on all objects
h1sets <- lapply(h1sets, function(obj) SCTransform(obj))

# find features for integration
h1.features <- SelectIntegrationFeatures(object.list = h1sets, nfeatures = 2000)
h1.list <- PrepSCTIntegration(object.list = h1sets, anchor.features = h1.features, 
                              verbose = FALSE)

# get anchors for integration
h1.anchors <- FindIntegrationAnchors(object.list = h1.list, anchor.features = h1.features, normalization.method = "SCT", dims = 1:30)
# integrate data
h1.combined <- IntegrateData(anchorset = h1.anchors, normalization.method = "SCT", dims = 1:30)
# change the default assay for downstream analysis to integrated uncorrected values
DefaultAssay(h1.combined) <- "RNA"
# variable features for the uncorrected integrated values
h1.combined <- FindVariableFeatures(h1.combined, selection.method = "vst", nfeatures = 2000)
# scale data
h1.combined <- ScaleData(h1.combined, verbose = FALSE)
# run pca
h1.combined <- RunPCA(h1.combined)
# find neighbors and clusters
h1.combined <- FindNeighbors(h1.combined, dims = 1:30)
h1.combined <- FindClusters(h1.combined, resolution = 0.5)
# t-SNE and Clustering
h1.combined <- RunTSNE(h1.combined, dims = 1:30)
p1 <- DimPlot(h1.combined, reduction = "tsne", split.by = "orig.ident")
plot(p1)
# save the rds
saveRDS(h1.combined, "C:\\Users\\scsac\\Desktop\\GATech\\GhosnLab\\SachinKumar_GhosnLab Shared folder\\Sample_5\\H1.rds")

# change the default assay for downstream analysis to integrated batch-corrected values
h1.combined <- readRDS("C:\\Users\\scsac\\Desktop\\GATech\\GhosnLab\\SachinKumar_GhosnLab Shared folder\\Sample_5\\H1.rds")
DefaultAssay(h1.combined) <- "integrated"
# run pca
h1.combined <- ScaleData(h1.combined, verbose = FALSE)
# run pca
h1.combined <- RunPCA(h1.combined)
# find neighbors and clusters
h1.combined <- FindNeighbors(h1.combined, dims = 1:30)
h1.combined <- FindClusters(h1.combined, resolution = 0.5)
# t-SNE and Clustering
h1.combined <- RunTSNE(h1.combined, dims = 1:30)
p2 <- DimPlot(h1.combined, reduction = "tsne", split.by = "orig.ident")
plot(p2)