library(harmony)
library(Seurat)
library(ggplot2)
library(sctransform)
library(future)
library(Matrix)
memory.limit(size = 12000000000)
matrix1_dir = "C:\\Users\\scsac\\Desktop\\GATech\\GhosnLab\\SachinKumar_GhosnLab Shared folder\\Sample_5\\5GEX_PBMC1_36yoF_F1\\Batch1\\"
matrix1.path <- paste0(matrix1_dir, "matrix.mtx.gz")
barcodes1.path <- paste0(matrix1_dir, "barcodes.tsv.gz")
features1.path <- paste0(matrix1_dir, "features.tsv.gz")
mat1 <- readMM(file = matrix1.path)
barcode1.names = read.delim(barcodes1.path, header = FALSE, stringsAsFactors = FALSE)
feature1.names = read.delim(features1.path, header = FALSE, stringsAsFactors = FALSE)
colnames(mat1) = barcode1.names$V1
rownames(mat1) = feature1.names$V1

matrix2_dir = "C:\\Users\\scsac\\Desktop\\GATech\\GhosnLab\\SachinKumar_GhosnLab Shared folder\\Sample_5\\5GEX_PBMC1_36yoF_F1\\Batch2\\"
matrix2.path <- paste0(matrix2_dir, "matrix.mtx.gz")
barcodes2.path <- paste0(matrix2_dir, "barcodes.tsv.gz")
features2.path <- paste0(matrix2_dir, "features.tsv.gz")
mat2 <- readMM(file = matrix2.path)
barcode2.names = read.delim(barcodes2.path, header = FALSE, stringsAsFactors = FALSE)
feature2.names = read.delim(features2.path, header = FALSE, stringsAsFactors = FALSE)
colnames(mat2) = barcode2.names$V1
rownames(mat2) = feature2.names$V1

matrix3_dir = "C:\\Users\\scsac\\Desktop\\GATech\\GhosnLab\\SachinKumar_GhosnLab Shared folder\\Sample_5\\5GEX_PBMC3_46yoM_H1\\Batch1\\"
matrix3.path <- paste0(matrix3_dir, "matrix.mtx.gz")
barcodes3.path <- paste0(matrix3_dir, "barcodes.tsv.gz")
features3.path <- paste0(matrix3_dir, "features.tsv.gz")
mat3 <- readMM(file = matrix3.path)
barcode3.names = read.delim(barcodes3.path, header = FALSE, stringsAsFactors = FALSE)
feature3.names = read.delim(features3.path, header = FALSE, stringsAsFactors = FALSE)
colnames(mat3) = barcode3.names$V1
rownames(mat3) = feature3.names$V1

matrix4_dir = "C:\\Users\\scsac\\Desktop\\GATech\\GhosnLab\\SachinKumar_GhosnLab Shared folder\\Sample_5\\5GEX_PBMC3_46yoM_H1\\Batch2\\"
matrix4.path <- paste0(matrix4_dir, "matrix.mtx.gz")
barcodes4.path <- paste0(matrix4_dir, "barcodes.tsv.gz")
features4.path <- paste0(matrix4_dir, "features.tsv.gz")
mat4 <- readMM(file = matrix4.path)
barcode4.names = read.delim(barcodes4.path, header = FALSE, stringsAsFactors = FALSE)
feature4.names = read.delim(features4.path, header = FALSE, stringsAsFactors = FALSE)
colnames(mat4) = barcode4.names$V1
rownames(mat4) = feature4.names$V1

rownames(barcode1.names) = barcode1.names$V1
rownames(barcode2.names) = barcode2.names$V1
rownames(barcode3.names) = barcode3.names$V1
rownames(barcode4.names) = barcode4.names$V1

dSD <- CreateSeuratObject(counts = cbind(mat1, mat2, mat3, mat4), project = "DSD", min.cells = 5) %>%
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(pc.genes = dSD@var.genes, npcs = 20, verbose = FALSE)

dSD@meta.data$data <- c(rep("F1B1", ncol(mat1)), rep("F1B2", ncol(mat2)), rep("H1B1", ncol(mat3)), rep("H1B2", ncol(mat4)))
options(repr.plot.height = 2.5, repr.plot.width = 6)
dSD <- dSD %>% 
  RunHarmony("data", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(dSD, 'harmony')
options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = dSD, reduction = "harmony", pt.size = .1, group.by = "data", do.return = TRUE)
plot(p1)
dSD <- dSD %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()
options(repr.plot.height = 4, repr.plot.width = 10)
DimPlot(dSD, reduction = "umap", group.by = "data", pt.size = .1, split.by = 'data')
options(repr.plot.height = 4, repr.plot.width = 6)
DimPlot(dSD, reduction = "umap", label = TRUE, pt.size = .1)