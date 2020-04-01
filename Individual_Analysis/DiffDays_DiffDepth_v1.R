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

# normalize - log normalize
f1sets <- lapply(f1sets, function(obj) NormalizeData(obj))

# Find variable features
f1sets <- lapply(f1sets, function(obj) FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000))

# Scale data - default scales on only variable features and not all genes
var.genes <- union(f1sets[[1]]@assays$RNA@var.features, f1sets[[2]]@assays$RNA@var.features)
f1sets <- lapply(f1sets, function(obj) ScaleData(obj, features = var.genes))
scale.data <- cbind(f1sets[[1]]@assays$RNA@scale.data[var.genes, ],
                    f1sets[[2]]@assays$RNA@scale.data[var.genes, ])
f1 <- subset(f1, cells = c(Cells(f1sets[[1]]), Cells(f1sets[[2]])))
f1@assays$RNA@var.features <- var.genes
f1@assays$RNA@scale.data <- scale.data[var.genes, Cells(f1)]

# Do PCA for the set
f1 <- RunPCA(f1)

# Find neighbors and clusters
f1 <- FindNeighbors(f1, dims = 1:12)
f1 <- FindClusters(f1, resolution = 0.5)

# Run UMAP/tSNE
f1 <- RunTSNE(f1, dims = 1:12)

# get the final cluster plots
cp1 <- DimPlot(f1, reduction = "tsne", split.by = "orig.ident")

plot(cp1)

# save the rds
saveRDS(f1, "C:\\Users\\scsac\\Desktop\\GATech\\GhosnLab\\SachinKumar_GhosnLab Shared folder\\Sample_5\\F1B1_2.rds")