library(data.table)
library(here)

##VARIABLE
script_dir = here()

set.seed(42)

designs = c("30_t15_17")

p_subs = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)


sub_sampling <- function(p, design, name_design) {
  
  num_subg = as.data.frame(table(design$subgroup), stringsAsFactors=FALSE)
  num_subg$Freq = 3 + round((num_subg$Freq - 3) * p)
  
  selected_sample = NULL
  for (subg in num_subg$Var1) {
    tmp = design[which(design$subgroup == subg),]
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
  
  design_ss = design[which(design$sample %in% sample_test), ]
  
  dir_ss = file.path(script_dir, "data", "design", "leucegene_ss", name_design, p)
  file_ss = file.path(dir_ss, "design.tsv")
  dir.create(dir_ss, recursive = TRUE, showWarnings = FALSE)
  write.table(design_ss, file_ss, quote=FALSE, row.names=FALSE, sep="\t")
}


foreach_design <- function(name_design, p_subs, script_dir) {
  dir_design = file.path(script_dir, "data", "design", "leucegene", designs[1])
  file_design = file.path(dir_design, "design.tsv")
  design = fread(file_design, header = TRUE, stringsAsFactors=FALSE, sep="\t", quote = "")
  
  lapply(p_subs, function(x) sub_sampling(x, design, name_design))
}

lapply(designs, function(x) foreach_design(x, p_subs, script_dir))

