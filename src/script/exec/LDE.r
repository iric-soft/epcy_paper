#Rscript --vanilla path_design path_counts
#ex:
# - [path to design]
# - [path to count]
#Rscript --vanilla /u/eaudemard/project/epcy_paper//src/script/exec/LDE.r /u/eaudemard/project/epcy_paper//data/design/leucegene/28_inv16_vs_28/cv/foldloo/dataset_0/train /u/eaudemard/project/epcy_paper//data/leucegene/STAR/ /u/eaudemard/project/epcy_paper//data/design/leucegene/28_inv16_vs_28/cv/foldloo/dataset_0/train/STAR/readcounts 4

library(limma)
library(edgeR)
library(DESeq2)

args = commandArgs(trailingOnly=TRUE)

path_design=args[1]
path_counts=args[2]
path_output=args[3]

dir.create(path_output, recursive = TRUE, showWarnings = FALSE)


design_file = file.path(path_design, "design.tsv")
design = read.table(design_file, header = TRUE, stringsAsFactors=FALSE, sep="\t")
sampleTable = data.frame(condition = design$subgroup)


file_counts = file.path(path_counts, "readcounts.xls")
counts = read.delim(file_counts, header=TRUE, row.names=1, check.names=FALSE)
design$sample = gsub("-", ".", design$sample)
counts = counts[, which(colnames(counts) %in% design$sample)]
counts = counts[, match(design$sample, colnames(counts))]

run_limma <- function(counts, design, path_output) {
  dge <- DGEList(counts=counts, group=design$subgroup)

  #keep <- filterByExpr(dge)
  #dge <- dge[keep,,keep.lib.sizes=FALSE]
  dge <- calcNormFactors(dge)

  design = model.matrix(~condition, data = sampleTable)
  v <- voom(dge, design, plot=F)
  fit <- lmFit(v, design)
  fit <- eBayes(fit)

  #To give more weight to FC in the ranking
  #fit <- treat(fit, lfc=log2(1.2))
  #topTreat(fit, coef=ncol(design))
  res = topTable(fit, n=Inf, sort.by="p")
  res = cbind(data.frame(ID=row.names(res), stringsAsFactors=FALSE), res)

  file_out = file.path(path_output, "limma_voom_genes.xls")
  write.table(res, file_out, quote=FALSE, row.names=FALSE, sep="\t")
}

run_deseq2 <- function(counts, design, path_output, num_proc) {
  library("BiocParallel")

  counts = trunc(counts)
  #rownames(counts) = unlist(lapply(rownames(counts), function(x) unlist(strsplit(x, "[.]"))[1]))

  dds <- DESeqDataSetFromMatrix(counts, sampleTable, ~condition)

  ## wald
  dds_wald <- DESeq(dds, parallel=TRUE, BPPARAM=MulticoreParam(num_proc))
  res <- results(dds_wald, contrast=c("condition","Query","Ref"))

  res = cbind(data.frame(ID=row.names(res), stringsAsFactors=FALSE), res)

  file_out = file.path(path_output, "deseq2_genes.xls")
  write.table(res, file_out, quote=FALSE, row.names=FALSE, sep="\t")
}

run_edger <- function(counts, design, path_output) {
  y <- DGEList(counts=counts)
  #keep <- filterByExpr(y)
  #y <- y[keep,,keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)

  design = model.matrix(~design$subgroup)

  y = estimateDisp(y,design)
  fit = glmQLFit(y,design)
  qlf = glmQLFTest(fit,coef=2)
  res_qlf = topTags(qlf, n=Inf)
  res_qlf = cbind(data.frame(ID=row.names(res_qlf), stringsAsFactors=FALSE), res_qlf)

  file_out = file.path(path_output, "edger_genes.xls")
  write.table(res_qlf, file_out, quote=FALSE, row.names=FALSE, sep="\t")
}

run_limma(counts, design, path_output)
run_edger(counts, design, path_output)
run_deseq2(counts, design, path_output, args[4])
