
##Bioconductor: https://bioconductor.org/
repo <- "http://cran.us.r-project.org"

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos=repo)

#BiocManager::install(version = "3.15")

BiocManager::install("GenomicDataCommons", ask=F)
BiocManager::install("SummarizedExperiment", ask=F)
BiocManager::install("DESeq2", ask=F)
BiocManager::install("BiocParallel", ask=F)
BiocManager::install("edgeR", ask=F)
BiocManager::install("limma", ask=F)
BiocManager::install("biomaRt", ask=F)
BiocManager::install("batchelor", ask=F)
BiocManager::install("scater", ask=F)
BiocManager::install("BiocSingular", ask=F)
BiocManager::install("scran", ask=F)
BiocManager::install("SingleCellExperiment", ask=F)
BiocManager::install("MAST", ask=F)

##Other

install.packages("data.table", repos=repo)
install.packages("dplyr", repos=repo)
install.packages("here", repos=repo)
install.packages("xlsx", repos=repo)

install.packages("devtools", repos=repo)
devtools::install_github("immunogenomics/presto")
