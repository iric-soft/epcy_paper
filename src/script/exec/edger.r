#Rscript --vanilla path_design path_counts
#ex:
# - [path to design]
# - [path to count]

library(edgeR)

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
counts = counts[, which(colnames(counts) %in% design$sample)]
counts = counts[, match(design$sample, colnames(counts))]
#rownames(counts) = unlist(lapply(rownames(counts), function(x) unlist(strsplit(x, "[.]"))[1]))

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
