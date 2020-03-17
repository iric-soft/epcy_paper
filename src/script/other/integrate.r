
library(data.table)
library(batchelor)
library(scater)
library(BiocSingular)
library(SingleCellExperiment)
library(scran)
library(Rtsne)
library(umap)

library(here)

#VARIABLE
script_dir = here()



file_leu = file.path(script_dir, "data", "leucegene_v2", "STAR", "tpm.xls")
leucegene = fread(file_leu, header = TRUE, stringsAsFactors=FALSE, sep="\t", quote = "")
setDF(leucegene)
rownames(leucegene) = leucegene$ID

file_leu = file.path(script_dir, "data", "design", "leucegene_v2", "all", "design.tsv")
leucegene_design = fread(file_leu, header = TRUE, stringsAsFactors=FALSE, sep="\t", quote = "")
setDF(leucegene_design)

file_tcga_aml = file.path(script_dir, "data", "TCGA_LAML", "htseq", "tpm.xls")
tcga_aml = fread(file_tcga_aml, header = TRUE, stringsAsFactors=FALSE, sep="\t", quote = "")
setDF(tcga_aml)
rownames(tcga_aml) = tcga_aml$ID

file_leu = file.path(script_dir, "data", "design", "TCGA_LAML", "all", "design.tsv")
tcga_aml_design = fread(file_leu, header = TRUE, stringsAsFactors=FALSE, sep="\t", quote = "")
setDF(tcga_aml_design)

leucegene = leucegene[which(leucegene$ID %in% tcga_aml$ID),]
tcga_aml = tcga_aml[which(tcga_aml$ID %in% leucegene$ID),]

leucegene = leucegene[, which(colnames(leucegene) %in% leucegene_design$sample)]
tcga_aml = tcga_aml[, which(colnames(tcga_aml) %in% tcga_aml_design$sample)]

mat_leucegene = as.matrix(leucegene)
mat_tcga_aml = as.matrix(tcga_aml)

mat_leucegene = log2(mat_leucegene + 1)
mat_tcga_aml = log2(mat_tcga_aml + 1)

out <- batchelor::fastMNN(mat_leucegene, mat_tcga_aml)

#out <- mnnCorrect(leucegene, tcga_aml)

sce <- out
assay(sce, "original") <- cbind(mat_leucegene, mat_tcga_aml)
sce$batch = replace(sce$batch, sce$batch == 1, "leucegene")
sce$batch = replace(sce$batch, sce$batch == 2, "TCGA-LAML")

set.seed(100)
osce <- runPCA(sce, exprs_values="original", ntop=Inf, BSPARAM=IrlbaParam())
osce <- runTSNE(osce, use_dimred="PCA")
ot <- plotTSNE(osce, colour_by="batch") + ggtitle("Original")

set.seed(100)
sce <- runTSNE(sce, use_dimred="corrected")
ct <- plotTSNE(sce, colour_by="batch") + ggtitle("Corrected")

multiplot(ot, ct, cols=2)

top=10
selected_ids = NULL
file="prediction_capability.xls"
for (str_design in c("6_NUP98NSD1", "13_EVI1", "18_t8_21", "18_Tri8", "23_Mono5", "30_t15_17", "31_inv16", "48_MLL", "87_Inter", "92_CK", "275_NK")) {
  file_comp = file.path(script_dir, "data", "design", "leucegene_v2", str_design, "STAR", "readcounts", file)
  data_comp = fread(file_comp, header = TRUE, stringsAsFactors=FALSE, sep="\t", quote = "")
  setDF(data_comp)
  data_comp = data_comp[order(data_comp$KERNEL_MCC, decreasing=TRUE), ]
  selected_ids = c(selected_ids, data_comp$ID[1:top])
}

leucegene_reduce = leucegene[which(rownames(leucegene) %in% selected_ids),]
tcga_aml_reduce = tcga_aml[which(rownames(tcga_aml) %in% selected_ids),]

mat_leucegene = as.matrix(leucegene_reduce)
mat_tcga_aml = as.matrix(tcga_aml_reduce)

mat_leucegene = log2(mat_leucegene + 1)
mat_tcga_aml = log2(mat_tcga_aml + 1)

out <- batchelor::fastMNN(mat_leucegene, mat_tcga_aml)

sce <- out
assay(sce, "original") <- cbind(mat_leucegene, mat_tcga_aml)
sce$batch = replace(sce$batch, sce$batch == 1, "leucegene")
sce$batch = replace(sce$batch, sce$batch == 2, "TCGA-LAML")
sce$subg = c(leucegene_design$subgroup, tcga_aml_design$sub_reduce)

set.seed(100)
osce <- runPCA(sce, exprs_values="original", ntop=Inf, BSPARAM=IrlbaParam())
osce <- runTSNE(osce, use_dimred="PCA")
ot <- plotTSNE(osce, colour_by="subg") + ggtitle("Original")

set.seed(100)
sce <- runTSNE(sce, use_dimred="corrected")
ct <- plotTSNE(sce, colour_by="subg") + ggtitle("Corrected")

multiplot(ot, ct, cols=2)




mat_all = cbind(mat_leucegene, mat_tcga_aml)
res_umap = umap(t(mat_leucegene), n_epochs=2000, random_state=42, min_dist=0.01)
df_res_umap = as.data.frame(res_umap$layout, stringsAsFactors=FALSE)
colnames(df_res_umap)=c("x", "y")
df_res_umap$subg = leucegene_design$subgroup
df_res_umap$cohort = c(rep("leucegene", length(leucegene_design$subgroup)))

ggplot(df_res_umap, aes(x=x, y=y, colour=cohort)) + geom_point(aes(shape=subg)) + scale_shape_manual(values=c(0:6,15,16,17,18,19))


mat_leucegene = as.matrix(leucegene_reduce)
mat_tcga_aml = as.matrix(tcga_aml_reduce)




mat_all = cbind(mat_leucegene, mat_tcga_aml)
res_umap = umap(t(mat_all), n_epochs=2000, random_state=42, min_dist=0.01)
df_res_umap = as.data.frame(res_umap$layout, stringsAsFactors=FALSE)
colnames(df_res_umap)=c("x", "y")
df_res_umap$subg = c(leucegene_design$subgroup, tcga_aml_design$sub_reduce)
df_res_umap$cohort = c(rep("leucegene", length(leucegene_design$subgroup)), rep("tcga_laml", length(tcga_aml_design$sub_reduce)))

ggplot(df_res_umap, aes(x=x, y=y, colour=cohort)) + geom_point(aes(shape=subg)) + scale_shape_manual(values=c(0:6,15,16,17,18,19))


out <- batchelor::fastMNN(mat_leucegene, mat_tcga_aml, k=40)
res_umap = umap(t(as.matrix(assay(out))), n_epochs=2000, random_state=42, min_dist=0.01)
df_res_umap = as.data.frame(res_umap$layout, stringsAsFactors=FALSE)
colnames(df_res_umap)=c("x", "y")
df_res_umap$subg = c(leucegene_design$subgroup, tcga_aml_design$sub_reduce)
df_res_umap$cohort = c(rep("leucegene", length(leucegene_design$subgroup)), rep("tcga_laml", length(tcga_aml_design$sub_reduce)))


ggplot(df_res_umap, aes(x=x, y=y, colour=cohort)) + geom_point(aes(shape=subg), size = 4) + scale_shape_manual(values=c(0:6,15,16,17,18,19))
