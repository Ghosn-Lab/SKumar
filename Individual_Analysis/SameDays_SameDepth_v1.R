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

# normalize - log normalize
pbmc_sets <- lapply(pbmc_sets, function(obj) NormalizeData(obj))

# Find variable features
pbmc_sets <- lapply(pbmc_sets, function(obj) FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000))

# Scale data - default scales on only variable features and not all genes
var.genes <- union(pbmc_sets[[1]]@assays$RNA@var.features, pbmc_sets[[2]]@assays$RNA@var.features)
var.genes <- union(var.genes, pbmc_sets[[3]]@assays$RNA@var.features)
pbmc_sets <- lapply(pbmc_sets, function(obj) ScaleData(obj, features = var.genes))
scale.data <- cbind(pbmc_sets[[1]]@assays$RNA@scale.data[var.genes, ],
                    pbmc_sets[[2]]@assays$RNA@scale.data[var.genes, ],
                    pbmc_sets[[3]]@assays$RNA@scale.data[var.genes, ])
pbmc <- subset(pbmc, cells = c(Cells(pbmc_sets[[1]]), Cells(pbmc_sets[[2]]), Cells(pbmc_sets[[3]])))
pbmc@assays$RNA@var.features <- var.genes
pbmc@assays$RNA@scale.data <- scale.data[var.genes, Cells(pbmc)]

# Do PCA for the set
pbmc <- RunPCA(pbmc)

# Find neighbors and clusters
pbmc <- FindNeighbors(pbmc, dims = 1:12)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Run UMAP/tSNE
pbmc <- RunTSNE(pbmc, dims = 1:12)

# get the final cluster plots
cp1 <- DimPlot(pbmc, reduction = "tsne", split.by = "orig.ident")

plot(cp1)

# save the rds
saveRDS(pbmc, "C:\\Users\\scsac\\Desktop\\GATech\\GhosnLab\\SachinKumar_GhosnLab Shared folder\\SameTissues_SameBatch\\PBMC.rds")