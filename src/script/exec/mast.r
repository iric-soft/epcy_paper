#Rscript --vanilla path_design path_counts
#ex:
# - [path to design]
# - [path to count]
library(edgeR)
library(MAST)

args = commandArgs(trailingOnly=TRUE)

path_design=args[1]
path_counts=args[2]
path_output=args[3]

dir.create(path_output, recursive = TRUE, showWarnings = FALSE)

design_file = file.path(path_design, "design.tsv")
design = read.table(design_file, header = TRUE, stringsAsFactors=FALSE, sep="\t")

file_counts = file.path(path_counts, "readcounts.xls")
counts = read.delim(file_counts, header=TRUE, row.names=1, check.names=FALSE)
counts = counts[, which(colnames(counts) %in% design$sample)]
counts = counts[, match(design$sample, colnames(counts))]
#rownames(counts) = unlist(lapply(rownames(counts), function(x) unlist(strsplit(x, "[.]"))[1]))

y = DGEList(counts=counts)
y = edgeR::calcNormFactors(y)
cpms = edgeR::cpm(y)
sca = FromMatrix(exprsArray = log2(cpms + 1),
                  cData = data.frame(wellKey = design$sample,
                                     grp = design$subgroup))
zlmdata = zlm(~grp, sca)
mast = lrTest(zlmdata, "grp")

df_res = data.frame(
  ID = names(mast[, "hurdle", "Pr(>Chisq)"],
  pval = mast[, "hurdle", "Pr(>Chisq)"]
  )
)
file_out = file.path(path_output, "mast_genes.xls")
write.table(df_res, file_out, quote=FALSE, row.names=FALSE, sep="\t")
