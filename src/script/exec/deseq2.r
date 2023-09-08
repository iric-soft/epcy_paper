#Rscript --vanilla path_design path_counts num_proc
#ex:
# - [path to design]
# - [path to count]
# - num_proc

library(DESeq2)
library("BiocParallel")

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
counts = trunc(counts)
#rownames(counts) = unlist(lapply(rownames(counts), function(x) unlist(strsplit(x, "[.]"))[1]))

dds <- DESeqDataSetFromMatrix(counts, sampleTable, ~condition)

## wald
dds_wald <- DESeq(dds, parallel=TRUE, BPPARAM=MulticoreParam(args[4]))
res <- results(dds_wald, contrast=c("condition","Query","Ref"))

res = cbind(data.frame(ID=row.names(res), stringsAsFactors=FALSE), res)
file_out = file.path(path_output, "deseq2_genes.xls")
write.table(res, file_out, quote=FALSE, row.names=FALSE, sep="\t")
