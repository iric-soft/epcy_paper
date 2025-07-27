library(rtracklayer)

library(here)



#VARIABLE
script_dir = here()
#time to DL the gtf
options(timeout=600)

df <- import("https://ftp.ensembl.org/pub/release-84/gtf/homo_sapiens/Homo_sapiens.GRCh38.84.gtf.gz")
df_gene <- as.data.frame(df[,c("gene_id", "gene_version", "gene_name", "gene_biotype")])
df_gene <- as.data.frame(df_gene[,c("gene_id", "gene_version", "gene_name", "gene_biotype")])
df_gene <- unique(df_gene)
df_gene$ensembl_gene_id = paste(df_gene$gene_id, df_gene$gene_version, sep=".")
df_gene = df_gene[,c(-1,-2)]
df_gene = df_gene[,c(3,1,2)]
colnames(df_gene) = c("ensembl_gene_id", "external_gene_name", "gene_biotype")

dir_out = file.path(script_dir, "data", "other")
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
file_out = file.path(dir_out, "GRCh38_84_genes_biotype.tsv")
write.table(df_gene, file_out, quote=FALSE, row.names=FALSE, sep="\t")


df <- import("https://ftp.ensembl.org/pub/release-98/gtf/homo_sapiens/Homo_sapiens.GRCh38.98.gtf.gz")
df_gene <- as.data.frame(df[,c("gene_id", "gene_version", "gene_name", "gene_biotype")])
df_gene <- as.data.frame(df_gene[,c("gene_id", "gene_version", "gene_name", "gene_biotype")])
df_gene <- unique(df_gene)
df_gene$ensembl_gene_id = paste(df_gene$gene_id, df_gene$gene_version, sep=".")
df_gene = df_gene[,c(-1,-2)]
df_gene = df_gene[,c(3,1,2)]
colnames(df_gene) = c("ensembl_gene_id", "external_gene_name", "gene_biotype")

dir_out = file.path(script_dir, "data", "other")
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
file_out = file.path(dir_out, "GRCh38_98_genes_biotype.tsv")
write.table(df_gene, file_out, quote=FALSE, row.names=FALSE, sep="\t")


# Old version, ensembl server are not available

# for GRCh38_84 host="http://mar2016.archive.ensembl.org"
# check here https://useast.ensembl.org/Help/ArchiveList to find host for an other version
#mart = biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host="https://mar2016.archive.ensembl.org")
#genes = biomaRt::getBM(attributes = c("ensembl_gene_id", "version", "external_gene_name", "gene_biotype"), mart=mart)
#genes$ensembl_gene_id = paste(genes$ensembl_gene_id, genes$version, sep=".")
#genes = genes[,-2]

#dir_out = file.path(script_dir, "data", "other")
#dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
#file_out = file.path(dir_out, "GRCh38_84_genes_biotype.tsv")
#write.table(genes, file_out, quote=FALSE, row.names=FALSE, sep="\t")

# for GRCh38_98 host="http://mar2016.archive.ensembl.org"
# check here https://useast.ensembl.org/Help/ArchiveList to find host for an other version
#mart = biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host="https://sep2019.archive.ensembl.org/")
#genes = biomaRt::getBM(attributes = c("ensembl_gene_id", "version", "external_gene_name", "gene_biotype"), mart=mart)
#genes$ensembl_gene_id = paste(genes$ensembl_gene_id, genes$version, sep=".")
#genes = genes[,-2]

#dir_out = file.path(script_dir, "data", "other")
#dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
#file_out = file.path(dir_out, "GRCh38_98_genes_biotype.tsv")
#write.table(genes, file_out, quote=FALSE, row.names=FALSE, sep="\t")
