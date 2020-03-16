library(liger)
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

colnames(mat1) <- paste("BM", colnames(mat1), sep = "_")
colnames(mat2) <- paste("FL", colnames(mat2), sep = "_")
colnames(mat3) <- paste("PBMC", colnames(mat3), sep = "_")

diffTissue.data = list(bm=mat1, fl=mat2, pbmc=mat3)
lObj.diffTissue <- createLiger(diffTissue.data)
lObj.diffTissue <- normalize(lObj.diffTissue)
lObj.diffTissue <- selectGenes(lObj.diffTissue, var.thresh = c(0.3, 0.875), do.plot = F)
lObj.diffTissue <- scaleNotCenter(lObj.diffTissue)
k.suggest <- suggestK(lObj.diffTissue, num.cores = 6, gen.new = T, plot.log2 = F,
                      nrep = 5)
lObj.diffTissue <- optimizeALS(lObj.diffTissue, k=25, thresh = 5e-5, nrep = 3)
lObj.diffTissue <- quantile_norm(lObj.diffTissue)
lObj.diffTissue <- runTSNE(lObj.diffTissue)
p_a <- plotByDatasetAndCluster(lObj.diffTissue, return.plots = T) 
# Modify plot output slightly
print(p_a[[1]])
print(p_a[[2]])