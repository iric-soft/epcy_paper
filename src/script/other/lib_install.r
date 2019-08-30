
##Bioconductor: https://bioconductor.org/

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("GenomicDataCommons", ask=F)
BiocManager::install("SummarizedExperiment", ask=F)
BiocManager::install("DESeq2", ask=F)
BiocManager::install("edgeR", ask=F)
BiocManager::install("limma", ask=F)
BiocManager::install("biomaRt", ask=F)
#For genes length 2nd solution 
#BiocManager::install("EnsDb.Hsapiens.v86", ask=F)

##Other

install.packages("data.table")
install.packages("dplyr")
install.packages("here")
