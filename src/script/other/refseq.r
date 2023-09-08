library(here)

#VARIABLE
script_dir = here()

# for GRCh38_84 host="http://mar2016.archive.ensembl.org"
# check here https://useast.ensembl.org/Help/ArchiveList to find host for an other version
mart = biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host="http://mar2017.archive.ensembl.org")
genes = biomaRt::getBM(attributes = c("ensembl_gene_id",
    "external_gene_name", "gene_biotype", "description"),
  mart=mart
)

path_design = "/u/eaudemard/project/epcy_paper/data/design/mader/"

for (design in c("ER", "HER2", "TNBC", "relapse", "relapse_ER", "relapse_HER2", "relapse_TNBC")) {
  file_epcy = file.path(path_design, design, "readcounts", "predictive_capability.xls")
  epcy = read.table(file_epcy, header=T, sep="\t", stringsAsFactors=F, check.names=F)

  epcy = merge(epcy, genes, by.x="id", by.y="ensembl_gene_id", all.x=T)
  write.table(epcy, file_epcy, quote=F, sep = "\t", row.names = FALSE, col.names = TRUE)
}
