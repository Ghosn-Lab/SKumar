Generate clusters for individual data sets before integration.\
File names explain the data set being worked on - 
1. DiffDays_DiffDepth - F1 data set (PBMC)
2. DiffDays_SameDepth - H1 data set (PBMC)
3. SameBatch_DiffTiss - BM and PBMC sequenced in the same batch
4. SameBatch_SameTissue - 3 PBMC samples sequenced in the same batch

Notes: Harmony based workflow. Processing done on merged Seurat objects until PCA.
Then uses Harmony embeddings for Cluster and Nearest neighbors.
Changes of code to use SCTransform instead of LogNormalize:
1. "Seurat::NormalizeData(obj)" -> "Seurat::SCTransform(obj)"
2. "RunHarmony("orig.ident", plot_convergence = TRUE)" -> "RunHarmony("orig.ident", assay.use = "SCT", plot_convergence = TRUE)"
3. No need to use FindVariableFeatures when using SCTransform. It is a built-in function for that.
