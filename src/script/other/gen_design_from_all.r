# ------------------------------------------------------------------------------
# Script Name: gen_design_from_all.r
#
# Description:
#   This script generates subgroup-specific design files from a master design file.
#   For each unique value in a specified subgroup column, it creates a new directory
#   and writes a design.tsv file where samples are labeled as "Query" (for the subgroup)
#   or "Ref" (for all others). Optionally, output directories can be named with the
#   number of samples in the subgroup.
#
# Usage:
#   Rscript --vanilla gen_design_from_all.r [path_to_design] [subgroup_column] [num_sample (optional)]
#
#   - path_to_design: Path to the directory containing the "all/design.tsv" file.
#   - subgroup_column: Column name in design.tsv used to define subgroups.
#   - num_sample (optional): If provided, output directories are named with the sample count.
#
# Example:
#   Rscript --vanilla gen_design_from_all.r /path/to/design Group TRUE
#
# Output:
#   For each subgroup, a directory is created with a design.tsv file labeling samples as "Query" or "Ref".
#
# Author: eric.audemard@umontreal.ca
# ------------------------------------------------------------------------------


args = commandArgs(trailingOnly=TRUE)
path_design=args[1]
subgroup_collumn=args[2]
num_sample = TRUE
if (length(args) == 3) {
  num_sample = args[3]
}

cat("Starting gen_design_from_all.r\n")
cat("Input design directory:", path_design, "\n")
cat("Subgroup column:", subgroup_collumn, "\n")
cat("Use sample count in directory name:", num_sample, "\n")

dir_design = file.path(path_design, "all")
file_in = file.path(dir_design, "design.tsv")
cat("Reading design file:", file_in, "\n")
df_all = read.table(file_in, header=TRUE, stringsAsFactors=FALSE, sep="\t", quote = "")


cat("Unique subgroups found:\n")
print(unique(df_all[[subgroup_collumn]]))

for (subg in unique(df_all[[subgroup_collumn]])) {
  if (subg == "Other") {
    cat("Skipping subgroup: Other\n")
    next
  }
  cat("\nProcessing subgroup:", subg, "\n")
  query_ids = which(df_all[[subgroup_collumn]] == subg)
  ref_ids = which(df_all[[subgroup_collumn]] != subg)

  dir_out = file.path(path_design, subg, sep="")
  if (num_sample) {
    dir_out = file.path(path_design, paste(length(query_ids), "_", subg, sep=""))
  }
  cat("Creating output directory:", dir_out, "\n")
  dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
  df_res = df_all
  if (length(df_all) > 5) {
    df_res = df_all[, c(1,2,3,4, 5)]
  }

  df_res[query_ids, 2] = "Query"
  df_res[ref_ids, 2] = "Ref"

  file_out = file.path(dir_out, "design.tsv")
  cat("Writing design file:", file_out, "\n")
  write.table(df_res, file_out, row.names = FALSE, col.names = TRUE, sep="\t", quote=FALSE)
}

cat("\nAll subgroups processed. Script finished.\n")