install.packages("Seurat")
library(Seurat)
library(ggplot2)
library(sctransform)
library(future)
library(Matrix)
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

bm_b <- CreateSeuratObject(counts = mat1, meta.data = barcode1.names, project = "bm")
bm_b <- SCTransform(bm_b, verbose = FALSE)
bm_b@meta.data$ctype <- "BM"

fl_b <- CreateSeuratObject(counts = mat2, meta.data = barcode2.names, project = "fl")
fl_b <- SCTransform(fl_b, verbose = FALSE)
fl_b@meta.data$ctype <- "FL"

pbmc_b <- CreateSeuratObject(counts = mat3, meta.data = barcode3.names, project = "pbmc")
pbmc_b <- SCTransform(pbmc_b, verbose = FALSE)
pbmc_b@meta.data$ctype <- "PBMC"

options(future.globals.maxSize = 4000 * 1024^2)

b_cell.features <- SelectIntegrationFeatures(object.list = list(bm_b, fl_b, pbmc_b), nfeatures = 3000)
b_cell.list <- PrepSCTIntegration(object.list = list(bm_b, fl_b, pbmc_b), anchor.features = b_cell.features, 
                                     verbose = FALSE)
b_cell.anchors <- FindIntegrationAnchors(object.list = b_cell.list, normalization.method = "SCT", dims = 1:30)
b_cell.combined <- IntegrateData(anchorset = b_cell.anchors, normalization.method = "SCT", dims = 1:30)

b_cell.combined <- ScaleData(b_cell.combined, verbose = FALSE)
b_cell.combined <- RunPCA(b_cell.combined, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
b_cell.combined <- RunUMAP(b_cell.combined, dims = 1:30)

p1 <- DimPlot(b_cell.combined, pt.size = 0.1, split.by = "ctype")
plot(p1)
