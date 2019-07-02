library(data.table)
library(dplyr)
library(here)

#VARIABLE
script_dir = here()

# for GRCh38_84 host="http://mar2016.archive.ensembl.org"
# check here https://useast.ensembl.org/Help/ArchiveList to find host for an other version
mart = biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host="http://mar2016.archive.ensembl.org")
e2g = biomaRt::getBM(attributes = c("ensembl_exon_id", "ensembl_gene_id", "version"), mart=mart)

exon_db = biomaRt::getBM(attributes = c("ensembl_exon_id", "exon_chrom_start", "exon_chrom_end"), mart = mart)

exon_db = merge(exon_db, e2g, by="ensembl_exon_id")
exon_db$ensembl_gene_id = paste(exon_db$ensembl_gene_id, exon_db$version, sep=".")

genes_length = lapply(
  unique(exon_db$ensembl_gene_id),
  function(x) {
    exon_selected = exon_db[which(exon_db$ensembl_gene_id == x),]
    res = exon_selected %>%
      arrange(exon_chrom_start) %>%
      group_by(g = cumsum(cummax(lag(exon_chrom_end, default = first(exon_chrom_end))) < exon_chrom_start)) %>%
      summarise(start = first(exon_chrom_start), end = max(exon_chrom_end))
    res$length = abs(res$start - res$end) + 1
    res = data.frame(ID=x, length=sum(res$length), stringsAsFactors=F)

    return(res)
  }
)

genes_length = rbindlist(genes_length)

dir_out = file.path(script_dir, "data", "other")
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
file_out = file.path(dir_out, "GRCh38_84_genes_length.tsv")
write.table(genes_length, file_out, quote=FALSE, row.names=FALSE, sep="\t")

library(EnsDb.Hsapiens.v86)
edb <- EnsDb.Hsapiens.v86
gene.lengths <- as.data.frame(lengthOf(edb, of = "gene"))
