library(Seurat)
library(dplyr)
library(patchwork)
options(future.globals.maxSize = 16000 * 1024^2)

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
bm_pbmc_sets <- c(subset(bm_pbmc, idents = "BM"), subset(bm_pbmc, idents = "PBMC"))
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
bm_pbmc_sets <- lapply(bm_pbmc_sets, function(obj) ThreshFilter(obj))

# do sctransform on all objects
bm_pbmc_sets <- lapply(bm_pbmc_sets, function(obj) SCTransform(obj))

# find features for integration
bm_pbmc.features <- SelectIntegrationFeatures(object.list = bm_pbmc_sets, nfeatures = 2000)
bm_pbmc.list <- PrepSCTIntegration(object.list = bm_pbmc_sets, anchor.features = bm_pbmc.features, 
                              verbose = FALSE)

# get anchors for integration
bm_pbmc.anchors <- FindIntegrationAnchors(object.list = bm_pbmc.list, anchor.features = bm_pbmc.features, normalization.method = "SCT", dims = 1:30)
# integrate data
bm_pbmc.combined <- IntegrateData(anchorset = bm_pbmc.anchors, normalization.method = "SCT", dims = 1:30)
# change the default assay for downstream analysis to integrated uncorrected values
DefaultAssay(bm_pbmc.combined) <- "RNA"
# variable features for the uncorrected integrated values
bm_pbmc.combined <- FindVariableFeatures(bm_pbmc.combined, selection.method = "vst", nfeatures = 2000)
# scale data
bm_pbmc.combined <- ScaleData(bm_pbmc.combined, verbose = FALSE)
# run pca
bm_pbmc.combined <- RunPCA(bm_pbmc.combined)
# find neighbors and clusters
bm_pbmc.combined <- FindNeighbors(bm_pbmc.combined, dims = 1:30)
bm_pbmc.combined <- FindClusters(bm_pbmc.combined, resolution = 0.5)
# t-SNE and Clustering
bm_pbmc.combined <- RunTSNE(bm_pbmc.combined, dims = 1:30)
p1 <- DimPlot(bm_pbmc.combined, reduction = "tsne", split.by = "orig.ident")
plot(p1)
# save the rds
saveRDS(bm_pbmc.combined, "C:\\Users\\scsac\\Desktop\\GATech\\GhosnLab\\SachinKumar_GhosnLab Shared folder\\DifferentTissues_sameBatch\\BM_PBMC_Integrated.rds")

# change the default assay for downstream analysis to integrated batch-corrected values
bm_pbmc.combined <- readRDS("C:\\Users\\scsac\\Desktop\\GATech\\GhosnLab\\SachinKumar_GhosnLab Shared folder\\DifferentTissues_sameBatch\\BM_PBMC_Integrated.rds")
DefaultAssay(bm_pbmc.combined) <- "integrated"
# run pca
bm_pbmc.combined <- ScaleData(bm_pbmc.combined, verbose = FALSE)
# run pca
bm_pbmc.combined <- RunPCA(bm_pbmc.combined)
# find neighbors and clusters
bm_pbmc.combined <- FindNeighbors(bm_pbmc.combined, dims = 1:30)
bm_pbmc.combined <- FindClusters(bm_pbmc.combined, resolution = 0.5)
# t-SNE and Clustering
bm_pbmc.combined <- RunTSNE(bm_pbmc.combined, dims = 1:30)
p2 <- DimPlot(bm_pbmc.combined, reduction = "tsne", split.by = "orig.ident")
plot(p2)