library(harmony)
library(Seurat)
library(cowplot)
library(future)
library(Matrix)
memory.limit(size = 12000000000)

matrix1_dir = "C:\\Users\\scsac\\Desktop\\GATech\\GhosnLab\\SachinKumar_GhosnLab Shared folder\\DifferentTissues_sameBatch\\BM B cells\\"
matrix1.path <- paste0(matrix1_dir, "matrix.mtx.gz")
barcodes1.path <- paste0(matrix1_dir, "barcodes.tsv")
features1.path <- paste0(matrix1_dir, "features.tsv")
mat1 <- readMM(file = matrix1.path)
barcode1.names = read.delim(barcodes1.path, header = FALSE, stringsAsFactors = FALSE)
feature1.names = read.delim(features1.path, header = FALSE, stringsAsFactors = FALSE)
colnames(mat1) = barcode1.names$V1
rownames(mat1) = feature1.names$V1

matrix2_dir = "C:\\Users\\scsac\\Desktop\\GATech\\GhosnLab\\SachinKumar_GhosnLab Shared folder\\DifferentTissues_sameBatch\\FL B cells\\"
matrix2.path <- paste0(matrix2_dir, "matrix.mtx")
barcodes2.path <- paste0(matrix2_dir, "barcodes.tsv")
genes2.path <- paste0(matrix2_dir, "genes.tsv")
mat2 <- readMM(file = matrix2.path)
barcode2.names = read.delim(barcodes2.path, header = FALSE, stringsAsFactors = FALSE)
gene2.names = read.delim(genes2.path, header = FALSE, stringsAsFactors = FALSE)
colnames(mat2) = barcode2.names$V1
rownames(mat2) = gene2.names$V1

matrix3_dir = "C:\\Users\\scsac\\Desktop\\GATech\\GhosnLab\\SachinKumar_GhosnLab Shared folder\\DifferentTissues_sameBatch\\PBMC B cells\\"
matrix3.path <- paste0(matrix3_dir, "matrix.mtx.gz")
barcodes3.path <- paste0(matrix3_dir, "barcodes.tsv")
features3.path <- paste0(matrix3_dir, "features.tsv")
mat3 <- readMM(file = matrix3.path)
barcode3.names = read.delim(barcodes3.path, header = FALSE, stringsAsFactors = FALSE)
feature3.names = read.delim(features3.path, header = FALSE, stringsAsFactors = FALSE)
colnames(mat3) = barcode3.names$V1
rownames(mat3) = feature3.names$V1

rownames(barcode1.names) = barcode1.names$V1
rownames(barcode2.names) = barcode2.names$V1
rownames(barcode3.names) = barcode3.names$V1
meta.data <- rbind(barcode1.names, barcode2.names, barcode3.names)

dT <- CreateSeuratObject(counts = cbind(mat1, mat2, mat3), project = "DT", min.cells = 5) %>%
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(pc.genes = dT@var.genes, npcs = 20, verbose = FALSE)

dT@meta.data$data <- c(rep("BM", ncol(mat1)), rep("FL", ncol(mat2)), rep("PBMC", ncol(mat3)))
options(repr.plot.height = 2.5, repr.plot.width = 6)
dT <- dT %>% 
  RunHarmony("data", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(dT, 'harmony')
options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = dT, reduction = "harmony", pt.size = .1, group.by = "data", do.return = TRUE)
plot(p1)
dT <- dT %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()
options(repr.plot.height = 4, repr.plot.width = 10)
DimPlot(dT, reduction = "umap", group.by = "data", pt.size = .1, split.by = 'data')
options(repr.plot.height = 4, repr.plot.width = 6)
DimPlot(dT, reduction = "umap", label = TRUE, pt.size = .1)