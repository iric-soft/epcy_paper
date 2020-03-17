library(data.table)
library(here)

##VARIABLE
script_dir = here()

set.seed(42)

dir_design = file.path(script_dir, "data", "design", "leucegene_v2", "all")
file_design = file.path(dir_design, "design.tsv")

design = fread(file_design, header = TRUE, stringsAsFactors=FALSE, sep="\t", quote = "")

num_subg = as.data.frame(table(design$subgroup), stringsAsFactors=FALSE)
num_subg$Freq = round(num_subg$Freq/3)
print(num_subg)

sample_test = NULL
for (subg in num_subg$Var1) {
  tmp = design[which(design$subgroup == subg),]
  sample_test = c(
    sample_test,
    tmp$sample[
      sample(
        1:length(tmp$subgroup),
        num_subg$Freq[which(num_subg$Var == subg)]
      )
    ]
  )
}


design_comp = design[which(design$sample %in% sample_test), ]
design_train_test = design[which(!design$sample %in% sample_test), ]

dir_train_test = file.path(script_dir, "data", "design", "leucegene_learning_curve", "test_train", "all")
file_train_test = file.path(dir_train_test, "design.tsv")
dir.create(dir_train_test, recursive = TRUE, showWarnings = FALSE)
write.table(design_train_test, file_train_test, quote=FALSE, row.names=FALSE, sep="\t")

dir_comp = file.path(script_dir, "data", "design", "leucegene_learning_curve", "comp", "all")
file_comp = file.path(dir_comp, "design.tsv")
dir.create(dir_comp, recursive = TRUE, showWarnings = FALSE)
write.table(design_comp, file_comp, quote=FALSE, row.names=FALSE, sep="\t")
