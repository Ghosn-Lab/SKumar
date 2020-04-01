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

# normalize - log normalize
bm_pbmc_sets <- lapply(bm_pbmc_sets, function(obj) NormalizeData(obj))

# Find variable features
bm_pbmc_sets <- lapply(bm_pbmc_sets, function(obj) FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000))

# Scale data - default scales on only variable features and not all genes
var.genes <- union(bm_pbmc_sets[[1]]@assays$RNA@var.features, bm_pbmc_sets[[2]]@assays$RNA@var.features)
bm_pbmc_sets <- lapply(bm_pbmc_sets, function(obj) ScaleData(obj, features = var.genes))
scale.data <- cbind(bm_pbmc_sets[[1]]@assays$RNA@scale.data[var.genes, ],
                    bm_pbmc_sets[[2]]@assays$RNA@scale.data[var.genes, ])
bm_pbmc <- subset(bm_pbmc, cells = c(Cells(bm_pbmc_sets[[1]]), Cells(bm_pbmc_sets[[2]])))
bm_pbmc@assays$RNA@var.features <- var.genes
bm_pbmc@assays$RNA@scale.data <- scale.data[var.genes, Cells(bm_pbmc)]

# Do PCA for the set
bm_pbmc <- RunPCA(bm_pbmc)

# Find neighbors and clusters
bm_pbmc <- FindNeighbors(bm_pbmc, dims = 1:12)
bm_pbmc <- FindClusters(bm_pbmc, resolution = 0.5)

# Run UMAP/tSNE
bm_pbmc <- RunTSNE(bm_pbmc, dims = 1:12)

# get the final cluster plots
cp1 <- DimPlot(bm_pbmc, reduction = "tsne", split.by = "orig.ident")

plot(cp1)

# save the rds
saveRDS(bm_pbmc, "C:\\Users\\scsac\\Desktop\\GATech\\GhosnLab\\SachinKumar_GhosnLab Shared folder\\DifferentTissues_sameBatch\\BM_PBMC.rds")