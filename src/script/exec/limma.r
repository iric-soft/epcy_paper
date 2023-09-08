#Rscript --vanilla path_design path_counts
#ex:
# - [path to design]
# - [path to count]


library(limma)
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
