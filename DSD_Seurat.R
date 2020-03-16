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

f1b1 <- CreateSeuratObject(counts = mat1, meta.data = barcode1.names, project = "F1B1")
f1b1 <- SCTransform(f1b1, verbose = FALSE)
f1b1@meta.data$batch <- "F1B1"

f1b2 <- CreateSeuratObject(counts = mat2, meta.data = barcode2.names, project = "F1B2")
f1b2 <- SCTransform(f1b2, verbose = FALSE)
f1b2@meta.data$batch <- "F1B2"

h1b1 <- CreateSeuratObject(counts = mat3, meta.data = barcode3.names, project = "H1B1")
h1b1 <- SCTransform(h1b1, verbose = FALSE)
h1b1@meta.data$batch <- "H1B1"

h1b2 <- CreateSeuratObject(counts = mat4, meta.data = barcode4.names, project = "H1B2")
h1b2 <- SCTransform(h1b2, verbose = FALSE)
h1b2@meta.data$batch <- "H1B2"

options(future.globals.maxSize = 16000 * 1024^2)

batches.features <- SelectIntegrationFeatures(object.list = list(f1b1, f1b2, h1b1, h1b2), nfeatures = 3000)
batches.list <- PrepSCTIntegration(object.list = list(f1b1, f1b2, h1b1, h1b2), anchor.features = batches.features, 
                                 verbose = FALSE)

batches.anchors <- FindIntegrationAnchors(object.list = batches.list, anchor.features = batches.features, normalization.method = "SCT", dims = 1:30)
batches.combined <- IntegrateData(anchorset = batches.anchors, normalization.method = "SCT", dims = 1:30)

batches.combined <- ScaleData(batches.combined, verbose = FALSE)
batches.combined <- RunPCA(batches.combined, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
batches.combined <- RunUMAP(batches.combined, dims = 1:30)
p1 <- DimPlot(batches.combined, split.by = "batch")
plot(p1)
