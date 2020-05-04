library(Matrix)
library(data.table)

matrix_dir = "./data/10X/filtered_feature_bc_matrix/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1

df = as.data.frame(as.matrix(mat))
setDT(df, keep.rownames = TRUE)[]
setnames(df, old = c("rn"), new = c("id"))

dir_out = "./data/10X/cellranger"
file_out = file.path(dir_out, "readcounts.xls")
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
write.table(df, file_out, quote=FALSE, row.names=FALSE, sep="\t")



analysis_dir = "./data/10X/analysis/clustering/graphclust/"
cluster_file = paste0(analysis_dir, "clusters.csv")
df_cluster = read.table(cluster_file, sep=",", quote="",
                        header=TRUE, stringsAsFactors = FALSE)
colnames(df_cluster) = c("sample", "subgroup")
dir_out = "./data/design/10X/all"
file_out = file.path(dir_out, "design.tsv")
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
write.table(df_cluster, file_out, quote=FALSE, row.names=FALSE, sep="\t")
