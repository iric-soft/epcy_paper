#Rscript --vanilla path_design path_counts
#ex:
# - [path to design]
# - [path to count]


library(presto)

args = commandArgs(trailingOnly=TRUE)

path_design=args[1]
path_counts=args[2]
path_output=args[3]

dir.create(path_output, recursive = TRUE, showWarnings = FALSE)


design_file = file.path(path_design, "design.tsv")
design = read.table(design_file, header = TRUE, stringsAsFactors=FALSE, sep="\t")
design = design[order(design$subgroup), ]
design$subgroup = factor(design$subgroup,levels=c("Ref", "Query"))

sampleTable = data.frame(condition = design$subgroup)


file_counts = file.path(path_counts, "readcounts.xls")
counts = read.delim(file_counts, header=TRUE, row.names=1, check.names=FALSE)
counts = counts[, which(colnames(counts) %in% design$sample)]
counts = counts[, match(design$sample, colnames(counts))]

res = presto::wilcoxauc(counts, design$subgroup)
colnames(res)[1] = "ID"
res_query = res[res$group=="Query",]

file_out = file.path(path_output, "presto_genes.xls")
write.table(res_query, file_out, quote=FALSE, row.names=FALSE, sep="\t")

#run_wilcox <- function(row, design) {
#    x=row[design$subgroup == 'Query']
#    y=row[design$subgroup=='Ref']

#    tmp = wilcox.test(x, y, alternative="two.sided")
#    return(tmp$p.value)
#}

#res_wilcox = apply(counts,1,run_wilcox, design=design)
