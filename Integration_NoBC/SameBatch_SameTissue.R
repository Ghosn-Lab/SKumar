library(Seurat)
library(dplyr)
library(patchwork)
options(future.globals.maxSize = 16000 * 1024^2)

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
pbmc_sets <- c(subset(pbmc, idents = "PBMC1"), subset(pbmc, idents = "PBMC2"), subset(pbmc, idents = "PBMC3"))
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
pbmc_sets <- lapply(pbmc_sets, function(obj) ThreshFilter(obj))

# do sctransform on all objects
pbmc_sets <- lapply(pbmc_sets, function(obj) SCTransform(obj))

# find features for integration
pbmc.features <- SelectIntegrationFeatures(object.list = pbmc_sets, nfeatures = 2000)
pbmc.list <- PrepSCTIntegration(object.list = pbmc_sets, anchor.features = pbmc.features, 
                                   verbose = FALSE)

# get anchors for integration
pbmc.anchors <- FindIntegrationAnchors(object.list = pbmc.list, anchor.features = pbmc.features, normalization.method = "SCT", dims = 1:30)
# integrate data
pbmc.combined <- IntegrateData(anchorset = pbmc.anchors, normalization.method = "SCT", dims = 1:30)
# change the default assay for downstream analysis to integrated uncorrected values
DefaultAssay(pbmc.combined) <- "RNA"
# variable features for the uncorrected integrated values
pbmc.combined <- FindVariableFeatures(pbmc.combined, selection.method = "vst", nfeatures = 2000)
# scale data
pbmc.combined <- ScaleData(pbmc.combined, verbose = FALSE)
# run pca
pbmc.combined <- RunPCA(pbmc.combined)
# find neighbors and clusters
pbmc.combined <- FindNeighbors(pbmc.combined, dims = 1:30)
pbmc.combined <- FindClusters(pbmc.combined, resolution = 0.5)
# t-SNE and Clustering
pbmc.combined <- RunTSNE(pbmc.combined, dims = 1:30)
p1 <- DimPlot(pbmc.combined, reduction = "tsne", split.by = "orig.ident")
plot(p1)
# save the rds
saveRDS(pbmc.combined, "C:\\Users\\scsac\\Desktop\\GATech\\GhosnLab\\SachinKumar_GhosnLab Shared folder\\SameTissues_SameBatch\\PBMC_Integrated.rds")

# change the default assay for downstream analysis to integrated batch-corrected values
pbmc.combined <- readRDS("C:\\Users\\scsac\\Desktop\\GATech\\GhosnLab\\SachinKumar_GhosnLab Shared folder\\SameTissues_SameBatch\\PBMC_Integrated.rds")
DefaultAssay(pbmc.combined) <- "integrated"
# run pca
pbmc.combined <- ScaleData(pbmc.combined, verbose = FALSE)
# run pca
pbmc.combined <- RunPCA(pbmc.combined)
# find neighbors and clusters
pbmc.combined <- FindNeighbors(pbmc.combined, dims = 1:30)
pbmc.combined <- FindClusters(pbmc.combined, resolution = 0.5)
# t-SNE and Clustering
pbmc.combined <- RunTSNE(pbmc.combined, dims = 1:30)
p2 <- DimPlot(pbmc.combined, reduction = "tsne", split.by = "orig.ident")
plot(p2)