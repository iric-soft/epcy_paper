# ------------------------------------------------------------------------------
# Script Name: gen_rep_leucegene_subsampling.r
#
# Description:
#   This script generates subsampled replicate design files for Leucegene data.
#   For a given design, it creates multiple subsampled versions at different
#   proportions (p_subs) and for multiple iterations (nums). Each subsample
#   is saved as a new design.tsv file. The script also replicates the readcount
#   matrix accordingly and saves it.
#
# Usage:
#   Rscript --vanilla gen_rep_leucegene_subsampling.r
#
# Input:
#   - Design files: data/design/leucegene3/[design]/design.tsv
#   - Readcount file: data/leucegene3/STAR_RSEM/readcounts.xls
#
# Output:
#   - Subsampled design files: data/design/leucegene3_rep_ss/[design]/[p]/[num]/design.tsv
#   - Replicated readcount file: data/leucegene3_rep/STAR_RSEM/readcounts.xls
#
# Author: eric.audemard@umontreal.ca
# ------------------------------------------------------------------------------

library(data.table)
library(here)

##VARIABLE
script_dir = here()

set.seed(42)

nums = c(1:10)
#designs = c("30_t15_17", "30_inv16")
design = "30_t15_17"
p_subs = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)

cat("Starting gen_rep_leucegene_subsampling.r\n")
cat("Script directory:", script_dir, "\n")
cat("Design:", design, "\n")
cat("Proportions to subsample:", paste(p_subs, collapse = ", "), "\n")
cat("Number of iterations:", length(nums), "\n")

sub_sampling <- function(p, design, name_design, num, script_dir, samples_rep) {
  cat("  Subsampling: p =", p, ", iteration =", num, "\n")
  num_subg = as.data.frame(table(design$subgroup), stringsAsFactors=FALSE)
  num_subg$Freq = 3 + round((num_subg$Freq - 3) * p)

  design$ntimes = 20
  design_rep = design[rep(seq_len(nrow(design)), design$ntimes),]
  design_rep$sample = samples_rep

  selected_sample = NULL
  for (subg in num_subg$Var1) {
    tmp = design_rep[which(design_rep$subgroup == subg),]
    cat("    Subgroup:", subg, "- selecting", num_subg$Freq[which(num_subg$Var == subg)], "samples\n")
    selected_sample = c(
      selected_sample,
      tmp$sample[
        sample(
          1:length(tmp$subgroup),
          num_subg$Freq[which(num_subg$Var == subg)]
        )
      ]
    )
  }
  design_ss = design_rep[which(design_rep$sample %in% selected_sample), ]
  design_ss = design_ss[,-3]

  dir_ss = file.path(script_dir, "data", "design", "leucegene3_rep_ss", name_design, p, num)
  file_ss = file.path(dir_ss, "design.tsv")
  cat("    Writing subsampled design to:", file_ss, "\n")
  dir.create(dir_ss, recursive = TRUE, showWarnings = FALSE)
  write.table(design_ss, file_ss, quote=FALSE, row.names=FALSE, sep="\t")
  return(design_ss$sample)
}

foreach_num <- function(num, design, name_design, p_subs, script_dir, samples_rep) {
  cat(" Iteration:", num, "\n")
  return(lapply(p_subs, function(x) sub_sampling(x, design, name_design, num, script_dir, samples_rep)))
}

foreach_design <- function(name_design, p_subs, script_dir, nums, samples_rep) {
  dir_design = file.path(script_dir, "data", "design", "leucegene3", name_design)
  file_design = file.path(dir_design, "design.tsv")
  cat("Reading design file:", file_design, "\n")
  design = fread(file_design, header = TRUE, stringsAsFactors=FALSE, sep="\t", quote = "")

  return(lapply(nums, function(x) foreach_num(x, design, name_design, p_subs, script_dir, samples_rep)))
}

file_mat = file.path(script_dir, "data", "leucegene3", "STAR_RSEM", "readcounts.xls")
cat("Reading readcount file:", file_mat, "\n")
df_mat = fread(file_mat, header = TRUE, stringsAsFactors=FALSE, sep="\t", quote = "")
setDF(df_mat)
cat("Replicating readcount matrix columns...\n")
df_mat_rep = df_mat[c(1,rep(2:ncol(df_mat), each = 20))]

#all_sample_used = lapply(designs,  function(x) foreach_design(x, p_subs, script_dir, nums, colnames(df_mat_rep)[-1])
cat("Generating subsampled designs...\n")
all_sample_used = foreach_design(design, p_subs, script_dir, nums, colnames(df_mat_rep)[-1])

dir_nat_rep = file.path(script_dir, "data", "leucegene3_rep", "STAR_RSEM")
file_rep = file.path(dir_nat_rep, "readcounts.xls")
cat("Writing replicated readcount matrix to:", file_rep, "\n")
dir.create(dir_nat_rep, recursive = TRUE, showWarnings = FALSE)
write.table(df_mat_rep, file_rep, quote=FALSE, row.names=FALSE, sep="\t")

cat("All subsampling completed.\n")