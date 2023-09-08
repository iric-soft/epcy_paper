#Rscript --vanilla path_design subgroup
#ex:
# - [path to design]
# - [subgroup]: column name used to select subgroups


#library(here)
#VARIABLE
#script_dir = here()

args = commandArgs(trailingOnly=TRUE)
path_design=args[1]
subgroup_collumn=args[2]
num_sample = TRUE
if (length(args) == 3) {
  num_sample = args[3]
}

dir_design = file.path(path_design, "all")
file_in = file.path(dir_design, "design.tsv")
df_all = read.table(file_in, header=TRUE, stringsAsFactors=FALSE, sep="\t", quote = "")

for (subg in unique(df_all[[subgroup_collumn]])) {
  query_ids = which(df_all[[subgroup_collumn]] == subg)
  ref_ids = which(df_all[[subgroup_collumn]] != subg)

  dir_out = file.path(path_design, subg, sep="")
  if (num_sample) {
    dir_out = file.path(path_design, paste(length(query_ids), "_", subg, sep=""))
  }
  dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
  df_res = df_all
  if (length(df_all) > 5) {
    df_res = df_all[, c(1,2,3,4, 5)]
  }

  df_res[query_ids, 2] = "Query"
  df_res[ref_ids, 2] = "Ref"

  file_out = file.path(dir_out, "design.tsv")
  write.table(df_res, file_out, row.names = FALSE, col.names = TRUE, sep="\t", quote=FALSE)
}
