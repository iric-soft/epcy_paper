library(data.table)
library(here)

##VARIABLE
script_dir = here()

set.seed(42)

dir_design = file.path(script_dir, "data", "design", "STAG2", "4458_ok")
file_design = file.path(dir_design, "design.tsv")

design = fread(file_design, header = TRUE, stringsAsFactors=FALSE, sep="\t", quote = "")

uniq_subg = unique(design$subgroup)

for (i in c(1:20)) {
  selected_subgroup = sample(1:length(uniq_subg), 4424, replace=T)
  tab_ss = as.matrix(table(selected_subgroup))

  samples_query = NULL
  for (j in 1:length(tab_ss)) {
    selected_design = design[which(design$subgroup == uniq_subg[as.numeric(rownames(tab_ss)[j])])]
    id_samples = sample(1:length(selected_design$sample), tab_ss[j], replace=F)
    samples_query = c(samples_query, selected_design$sample[id_samples])
  }

  design_rand = design
  design_rand$subgroup[which(design$sample %in% samples_query)] = "Query"
  design_rand$subgroup[which(!design$sample %in% samples_query)] = "Ref"

  dir_rand = file.path(script_dir, "data", "design", "STAG2_random", i)
  file_rand = file.path(dir_rand, "design.tsv")
  dir.create(dir_rand, recursive = TRUE, showWarnings = FALSE)
  write.table(design_rand, file_rand, quote=FALSE, row.names=FALSE, sep="\t")
}
