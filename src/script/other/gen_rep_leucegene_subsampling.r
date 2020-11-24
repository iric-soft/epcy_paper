library(data.table)
library(here)

##VARIABLE
script_dir = here()


set.seed(42)

nums = c(1:10)
designs = c("30_t15_17", "30_inv16")
p_subs = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)


sub_sampling <- function(p, design, name_design, num, script_dir, samples_rep) {

  num_subg = as.data.frame(table(design$subgroup), stringsAsFactors=FALSE)
  num_subg$Freq = 3 + round((num_subg$Freq - 3) * p)

  design$ntimes = 20
  design_rep = design[rep(seq_len(nrow(design)), design$ntimes),]
  design_rep$sample = samples_rep

  selected_sample = NULL
  for (subg in num_subg$Var1) {
    tmp = design_rep[which(design_rep$subgroup == subg),]
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
  dir.create(dir_ss, recursive = TRUE, showWarnings = FALSE)
  write.table(design_ss, file_ss, quote=FALSE, row.names=FALSE, sep="\t")
  return(design_ss$sample)
}

foreach_num <- function(num, design, name_design, p_subs, script_dir, samples_rep) {

  return(lapply(p_subs, function(x) sub_sampling(x, design, name_design, num, script_dir, samples_rep)))
}

foreach_design <- function(name_design, p_subs, script_dir, nums, samples_rep) {
  dir_design = file.path(script_dir, "data", "design", "leucegene3", designs[1])
  file_design = file.path(dir_design, "design.tsv")
  design = fread(file_design, header = TRUE, stringsAsFactors=FALSE, sep="\t", quote = "")

  return(lapply(nums, function(x) foreach_num(x, design, name_design, p_subs, script_dir, samples_rep)))
}







file_mat = file.path(script_dir, "data", "leucegene3", "STAR_RSEM", "readcounts.xls")
df_mat = fread(file_mat, header = TRUE, stringsAsFactors=FALSE, sep="\t", quote = "")
setDF(df_mat)
df_mat_rep = df_mat[c(1,rep(2:ncol(df_mat), each = 20))]

all_sample_used = lapply(designs, function(x) foreach_design(x, p_subs, script_dir, nums, colnames(df_mat_rep)))

dir_nat_rep = file.path(script_dir, "data", "leucegene3_rep", "STAR_RSEM")
file_rep = file.path(dir_nat_rep, "readcounts.xls")
dir.create(dir_nat_rep, recursive = TRUE, showWarnings = FALSE)
write.table(df_mat_rep, file_rep, quote=FALSE, row.names=FALSE, sep="\t")
