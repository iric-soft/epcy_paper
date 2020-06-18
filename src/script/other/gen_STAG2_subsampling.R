library(data.table)
library(here)

##VARIABLE
script_dir = here()

set.seed(42)

nums = c(1:10)
p_subs = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)


sub_sampling <- function(p, design, num, script_dir) {

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

  design_ss = design[which(design$sample %in% selected_sample), ]

  dir_ss = file.path(script_dir, "data", "design", "STAG2_ss", "2662_ko", p, num)
  file_ss = file.path(dir_ss, "design.tsv")
  dir.create(dir_ss, recursive = TRUE, showWarnings = FALSE)
  write.table(design_ss, file_ss, quote=FALSE, row.names=FALSE, sep="\t")
}

foreach_num <- function(num, design, p_subs, script_dir) {

  lapply(p_subs, function(x) sub_sampling(x, design, num, script_dir))
}

dir_design = file.path(script_dir, "data", "design", "STAG2", "2662_ko")
file_design = file.path(dir_design, "design.tsv")
design = fread(file_design, header = TRUE, stringsAsFactors=FALSE, sep="\t", quote = "")

lapply(nums, function(x) foreach_num(x, design, p_subs, script_dir))
