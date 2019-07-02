
##Bioconductor: https://bioconductor.org/

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("GenomicDataCommons")
BiocManager::install("SummarizedExperiment")
BiocManager::install("DESeq2")
BiocManager::install("edgeR")
BiocManager::install("limma")
BiocManager::install("biomaRt")

##Other

install.packages("data.table")
install.packages("dplyr")
