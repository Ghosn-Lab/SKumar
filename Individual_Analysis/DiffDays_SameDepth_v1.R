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

# normalize - log normalize
h1sets <- lapply(h1sets, function(obj) NormalizeData(obj))

# Find variable features
h1sets <- lapply(h1sets, function(obj) FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000))

# Scale data - default scales on only variable features and not all genes
var.genes <- union(h1sets[[1]]@assays$RNA@var.features, h1sets[[2]]@assays$RNA@var.features)
h1sets <- lapply(h1sets, function(obj) ScaleData(obj, features = var.genes))
scale.data <- cbind(h1sets[[1]]@assays$RNA@scale.data[var.genes, ],
                    h1sets[[2]]@assays$RNA@scale.data[var.genes, ])
h1 <- subset(h1, cells = c(Cells(h1sets[[1]]), Cells(h1sets[[2]])))
h1@assays$RNA@var.features <- var.genes
h1@assays$RNA@scale.data <- scale.data[var.genes, Cells(h1)]

# Do PCA for the set
h1 <- RunPCA(h1)

# Find neighbors and clusters
h1 <- FindNeighbors(h1, dims = 1:12)
h1 <- FindClusters(h1, resolution = 0.5)

# Run UMAP/tSNE
h1 <- RunTSNE(h1, dims = 1:12)

# get the final cluster plots
cp1 <- DimPlot(h1, reduction = "tsne", split.by = "orig.ident")

plot(cp1)

# save the rds
saveRDS(h1, "C:\\Users\\scsac\\Desktop\\GATech\\GhosnLab\\SachinKumar_GhosnLab Shared folder\\Sample_5\\H1B1_2.rds")