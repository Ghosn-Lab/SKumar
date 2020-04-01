library(Seurat)
library(dplyr)
library(patchwork)
options(future.globals.maxSize = 16000 * 1024^2)
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

f1 <- merge(f1b1_obj, f1b2_obj, add.cell.ids = c("F1B1", "F1B2"), project = "F1") # make a set for running future functions
f1sets <- c(subset(f1, idents = "F1B1"), subset(f1, idents = "F1B2"))
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
f1sets <- lapply(f1sets, function(obj) ThreshFilter(obj))
# do sctransform on all objects
f1sets <- lapply(f1sets, function(obj) SCTransform(obj))

# find features for integration
f1.features <- SelectIntegrationFeatures(object.list = f1sets, nfeatures = 2000)
f1.list <- PrepSCTIntegration(object.list = f1sets, anchor.features = f1.features, 
                                   verbose = FALSE)

# get anchors for integration
f1.anchors <- FindIntegrationAnchors(object.list = f1.list, anchor.features = f1.features, normalization.method = "SCT", dims = 1:30)
# integrate data
f1.combined <- IntegrateData(anchorset = f1.anchors, normalization.method = "SCT", dims = 1:30)
# change the default assay for downstream analysis to integrated uncorrected values
DefaultAssay(f1.combined) <- "RNA"
# variable features for the uncorrected integrated values
f1.combined <- FindVariableFeatures(f1.combined, selection.method = "vst", nfeatures = 2000)
# scale the data
f1.combined <- ScaleData(f1.combined, verbose = FALSE)
# run pca
f1.combined <- RunPCA(f1.combined)
# find neighbors and clusters
f1.combined <- FindNeighbors(f1.combined, dims = 1:30)
f1.combined <- FindClusters(f1.combined, resolution = 0.5)
# t-SNE and Clustering
f1.combined <- RunTSNE(f1.combined, dims = 1:30)
p1 <- DimPlot(f1.combined, reduction = "tsne", split.by = "orig.ident")
plot(p1)
# save the rds
saveRDS(f1.combined, "C:\\Users\\scsac\\Desktop\\GATech\\GhosnLab\\SachinKumar_GhosnLab Shared folder\\Sample_5\\F1.rds")

# change the default assay for downstream analysis to integrated batch-corrected values
f1.combined <- readRDS("C:\\Users\\scsac\\Desktop\\GATech\\GhosnLab\\SachinKumar_GhosnLab Shared folder\\Sample_5\\F1.rds")
DefaultAssay(f1.combined) <- "integrated"
# scale the data
f1.combined <- ScaleData(f1.combined, verbose = FALSE)
# run pca
f1.combined <- RunPCA(f1.combined)
# find neighbors and clusters
f1.combined <- FindNeighbors(f1.combined, dims = 1:30)
f1.combined <- FindClusters(f1.combined, resolution = 0.5)
# t-SNE and Clustering
f1.combined <- RunTSNE(f1.combined, dims = 1:30)
p2 <- DimPlot(f1.combined, reduction = "tsne", split.by = "orig.ident")
plot(p2)