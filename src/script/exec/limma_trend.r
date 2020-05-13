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
counts = read.delim(file_counts, header=TRUE, row.names=1)
colnames(counts) = gsub("^X", "", colnames(counts))
design$sample = gsub("-", ".", design$sample)
counts = counts[, which(colnames(counts) %in% design$sample)]
counts = counts[, match(design$sample, colnames(counts))]
#rownames(counts) = unlist(lapply(rownames(counts), function(x) unlist(strsplit(x, "[.]"))[1]))

dge <- DGEList(counts=counts, group=design$subgroup)

#keep <- filterByExpr(dge)
#dge <- dge[keep,,keep.lib.sizes=FALSE]
#dge <- calcNormFactors(dge)
logCPM <- cpm(dge, log=TRUE, prior.count=3)

design = model.matrix(~condition, data = sampleTable)
fit <- lmFit(logCPM, design)
fit <- eBayes(fit, trend=TRUE)
topTable(fit, coef=ncol(design), n=Inf, sort.by="p")

#To give more weight to FC in the ranking
#fit <- treat(fit, lfc=log2(1.2))
#topTreat(fit, coef=ncol(design))
res = topTable(fit, n=Inf, sort.by="p")
res = cbind(data.frame(ID=row.names(res), stringsAsFactors=FALSE), res)
file_out = file.path(path_output, "limma_trend_genes.xls")
write.table(res, file_out, quote=FALSE, row.names=FALSE, sep="\t")

############################################
############################################
############################################

# Function to run limma trend from powsimR
.run.limma.trend <- function(normData, countData, DEOpts, verbose) {

  # calculate normalisation factors
  sf <- normData$size.factors
  sf[sf<0] <- min(sf[sf > 0])
  nsf <- log(sf/colSums(countData))
  nsf <- exp(nsf - mean(nsf, na.rm=T))

  if (attr(normData, 'normFramework')  %in% c('sctransform')) {
    norm.cnts <- normData$RoundNormCounts
    ixx.valid <- rownames(countData) %in% rownames(norm.cnts)
    countData[ixx.valid, ] <- norm.cnts
  }

  # construct input object
  dge <- edgeR::DGEList(counts = countData,
                        lib.size = colSums(countData),
                        norm.factors = nsf,
                        group = factor(DEOpts$designs),
                        remove.zeros = FALSE)
  # run DE testing
  p.DE <- DEOpts$p.DE
  design.mat <- stats::model.matrix( ~ DEOpts$designs)
  y <- new("EList")
  y$E <- edgeR::cpm(dge, log = TRUE, prior.count = 3)

  if (attr(normData, 'normFramework') %in% c('SCnorm', 'Linnorm')) {
    scale.facts <- normData$scale.factors
    ixx.valid <- rownames(countData) %in% rownames(scale.facts)
    wgenes <- countData
    wgenes[1:nrow(wgenes), 1:ncol(wgenes)] <- NA
    wgenes[ixx.valid, ] <- scale.facts
    fit <- limma::lmFit(object = y, design = design.mat, weights = wgenes)
  }
  if (!attr(normData, 'normFramework') %in% c('SCnorm', 'Linnorm')) {
    fit <- limma::lmFit(object = y, design = design.mat)
  }

  fit <- limma::eBayes(fit, trend=TRUE, proportion=p.DE, robust=TRUE)
  resT <- limma::topTable(fit=fit, coef=2, number=Inf, adjust.method = "BH", sort.by = "none")

  # construct results
  result <- data.frame(geneIndex=rownames(resT),
                       pval=resT$P.Value,
                       fdr=rep(NA, nrow(resT)),
                       lfc=resT$logFC,
                       stringsAsFactors = F)
  return(result)
}
