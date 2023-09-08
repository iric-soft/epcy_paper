
library(ggplot2)
library(here)

#VARIABLE
script_dir = here()

subg = "inv16"
quant = "kallisto"
dir_design = file.path(script_dir, "data", "design", "leucegene_comp_tt")

subset = "comp"
dir_a =  file.path(dir_design, subset, subg, quant, "tpm")
file_a = file.path(dir_a, "prediction_capability.xls")
tab_a = read.table(file_a, header=T, stringsAsFactors=F)
tab_a = tab_a[order(tab_a$kernel_mcc, abs(tab_a$l2fc), decreasing=TRUE),]
tab_a$order = 1:length(tab_a$kernel_mcc)

subset = "test_train"
dir_b =  file.path(dir_design, subset, subg, quant, "tpm")
file_b = file.path(dir_b, "prediction_capability.xls")
tab_b = read.table(file_b, header=T, stringsAsFactors=F)
tab_b = tab_b[order(tab_b$kernel_mcc, abs(tab_b$l2fc), decreasing=TRUE),]
tab_b$order = 1:length(tab_b$kernel_mcc)

ids_selected = unique(c(tab_a$id[1:100], tab_b$id[1:100]))


tab_merge = merge(tab_a, tab_b, by="id", all=T)
tab_merge = tab_merge[which(tab_merge$id %in% ids_selected),]


gg = ggplot(tab_merge, aes(x=kernel_mcc.x, y=kernel_mcc.y)) + geom_point(alpha=0.5)
gg

gg = ggplot(tab_merge, aes(x=order.x, y=order.y)) + geom_point(alpha=1)
gg



subset = "comp"
dir_a =  file.path(dir_design, subset, subg, quant, "tpm10_20")
file_a = file.path(dir_a, "prediction_capability.xls")
tab_a = read.table(file_a, header=T, stringsAsFactors=F)
tab_a = tab_a[order(tab_a$kernel_mcc, abs(tab_a$l2fc), decreasing=TRUE),]
tab_a = tab_a[1:1000,]
tab_a$order = 1:length(tab_a$kernel_mcc)

subset = "test_train"
dir_b =  file.path(dir_design, subset, subg, quant, "tpm10_20")
file_b = file.path(dir_b, "prediction_capability.xls")
tab_b = read.table(file_b, header=T, stringsAsFactors=F)
tab_b = tab_b[order(tab_b$kernel_mcc, abs(tab_b$l2fc), decreasing=TRUE),]
tab_b = tab_b[1:1000,]
tab_b$order = 1:length(tab_b$kernel_mcc)

tab_merge = merge(tab_a, tab_b, by="id", all=T)
tab_merge[is.na(tab_merge)] = 1001
tab_merge = tab_merge[order(tab_merge$kernel_mcc.x, abs(tab_merge$l2fc.x), decreasing=TRUE),]


gg = ggplot(tab_merge, aes(x=kernel_mcc.x, y=kernel_mcc.y)) + geom_point(alpha=0.3)
gg

gg = ggplot(tab_merge, aes(x=order.x, y=order.y)) + geom_point(alpha=1)
gg




quant = "STAR"
subset = "comp"
dir_a =  file.path(dir_design, subset, subg, quant, "readcounts")
file_a = file.path(dir_a, "prediction_capability.xls")
tab_a = read.table(file_a, header=T, stringsAsFactors=F)
tab_a = tab_a[order(tab_a$kernel_mcc, abs(tab_a$l2fc), decreasing=TRUE),]
tab_a$order = 1:length(tab_a$kernel_mcc)

subset = "test_train"
dir_b =  file.path(dir_design, subset, subg, quant, "readcounts")
file_b = file.path(dir_b, "prediction_capability.xls")
tab_b = read.table(file_b, header=T, stringsAsFactors=F)
tab_b = tab_b[order(tab_b$kernel_mcc, abs(tab_b$l2fc), decreasing=TRUE),]
tab_b$order = 1:length(tab_b$kernel_mcc)

ids_selected = unique(c(tab_a$id[1:100], tab_b$id[1:100]))


tab_merge = merge(tab_a, tab_b, by="id", all=T)
tab_merge = tab_merge[which(tab_merge$id %in% ids_selected),]


gg = ggplot(tab_merge, aes(x=kernel_mcc.x, y=kernel_mcc.y)) + geom_point(alpha=0.8)
gg

gg = ggplot(tab_merge, aes(x=order.x, y=order.y)) + geom_point(alpha=1)
gg






run_comp <- function(subg, quant, dir_design, epcy_folder) {
  #### EPCY
  subset = "comp"
  #dir_epcy = file.path(dir_design, subset, subg, quant, "readcounts_bagging")
  dir_epcy = file.path(dir_design, subset, subg, quant, epcy_folder)
  file_epcy = file.path(dir_epcy, "prediction_capability.xls")
  tab_epcy = read.table(file_epcy, header=T, stringsAsFactors=F)
  tab_epcy$id = unlist(lapply(tab_epcy$id, function(x) unlist(strsplit(x, "[.]"))[1]))
  subset = "test_train"
  #dir_ref_epcy = file.path(dir_design, subset, subg, quant, "readcounts_bagging")
  dir_ref_epcy = file.path(dir_design, subset, subg, quant, epcy_folder)
  file_ref_epcy = file.path(dir_ref_epcy, "prediction_capability.xls")
  tab_ref_epcy = read.table(file_ref_epcy, header=T, stringsAsFactors=F)
  tab_ref_epcy$id = unlist(lapply(tab_ref_epcy$id, function(x) unlist(strsplit(x, "[.]"))[1]))

  #tab_ref_epcy = tab_ref_epcy[which(tab_ref_epcy$mean_query + tab_ref_epcy$mean_ref >=0.5),]
  #tab_epcy = tab_epcy[which(tab_epcy$mean_query + tab_epcy$mean_ref >=0.5),]

  tab_epcy = tab_epcy[order(tab_epcy$kernel_mcc, abs(tab_epcy$l2fc), decreasing=TRUE),]
  tab_epcy$order = 1:length(tab_epcy$kernel_mcc)
  tab_ref_epcy = tab_ref_epcy[order(tab_ref_epcy$kernel_mcc, abs(tab_ref_epcy$l2fc), decreasing=TRUE),]
  tab_ref_epcy$order = 1:length(tab_ref_epcy$kernel_mcc)
  #epcy_merge = epcy_merge[order(epcy_merge$kernel_mcc.x, abs(epcy_merge$l2fc.x), decreasing=TRUE),]

  #tab_epcy_filtred = tab_epcy[1:5,]
  #tab_ref_epcy_filtred = tab_ref_epcy[1:5,]
  #epcy_merge2 = merge(tab_epcy_filtred, tab_ref_epcy_filtred, by="id", all=F)
  #epcy_merge3 = merge(tab_epcy_filtred, tab_ref_epcy_filtred, by="id", all=T)
  #length(epcy_merge2$id) / length(epcy_merge3$id)


  #### Limma
  quant = "STAR"
  subset = "comp"
  dir_limma = file.path(dir_design, subset, subg, quant, "readcounts")
  file_limma = file.path(dir_limma, "limma_voom_genes.xls")
  tab_limma = read.table(file_limma, header=T, stringsAsFactors=F)
  tab_limma$ID = unlist(lapply(tab_limma$ID, function(x) unlist(strsplit(x, "[.]"))[1]))
  tab_limma = tab_limma[which(tab_limma$P.Value<=0.05), ]
  tab_limma = tab_limma[order(abs(tab_limma$logFC), decreasing=TRUE),]
  tab_limma$order = 1:length(tab_limma$ID)

  #subset = "test_train"
  #dir_ref_limma = file.path(dir_design, subset, subg, quant, "readcounts")
  #file_ref_limma  = file.path(dir_ref_limma , "limma_voom_genes.xls")
  #tab_ref_limma = read.table(file_ref_limma, header=T, stringsAsFactors=F)
  #tab_ref_limma$ID = unlist(lapply(tab_ref_limma$ID, function(x) unlist(strsplit(x, "[.]"))[1]))
  #limma_merge = merge(tab_ref_limma, tab_limma, by="ID")
  #tab_ref_limma = tab_ref_limma[which(tab_ref_limma$P.Value<=0.05), ]
  #tab_ref_limma = tab_ref_limma[order(abs(tab_ref_limma$logFC), decreasing=TRUE),]
  #tab_ref_limma$order = 1:length(tab_ref_limma$ID)
  #limma_merge2 = merge(tab_ref_limma, tab_limma, by="ID", all=F)
  #limma_merge3 = merge(tab_ref_limma, tab_limma, by="ID", all=T)
  #length(limma_merge2$ID) / length(limma_merge3$ID)

  #### DESEQ2
  subset = "comp"
  dir_deseq2 = file.path(dir_design, subset, subg, quant, "readcounts")
  file_deseq2 = file.path(dir_deseq2, "deseq2_genes.xls")
  tab_deseq2 = read.table(file_deseq2, header=T, stringsAsFactors=F)
  tab_deseq2$ID = unlist(lapply(tab_deseq2$ID, function(x) unlist(strsplit(x, "[.]"))[1]))
  tab_deseq2 = tab_deseq2[which(tab_deseq2$pvalue<=0.05), ]
  tab_deseq2 = tab_deseq2[order(abs(tab_deseq2$log2FoldChange), decreasing=TRUE),]
  tab_deseq2$order = 1:length(tab_deseq2$ID)

  #subset = "test_train"
  #dir_ref_deseq2 = file.path(dir_design, subset, subg, quant, "readcounts")
  #file_ref_deseq2 = file.path(dir_ref_deseq2, "deseq2_genes.xls")
  #tab_ref_deseq2 = read.table(file_ref_deseq2, header=T, stringsAsFactors=F)
  #tab_ref_deseq2$ID = unlist(lapply(tab_ref_deseq2$ID, function(x) unlist(strsplit(x, "[.]"))[1]))
  #deseq2_merge = merge(tab_ref_deseq2, tab_deseq2, by="ID")
  #tab_ref_deseq2 = tab_ref_deseq2[which(tab_ref_deseq2$pvalue<=0.05), ]
  #tab_ref_deseq2 = tab_ref_deseq2[order(abs(tab_ref_deseq2$log2FoldChange), decreasing=TRUE),]
  #tab_ref_deseq2$order = 1:length(tab_ref_deseq2$ID)
  #deseq2_merge2 = merge(tab_ref_deseq2, tab_deseq2, by="ID", all=F)
  #deseq2_merge3 = merge(tab_ref_deseq2, tab_deseq2, by="ID", all=T)
  #length(deseq2_merge2$ID) / length(deseq2_merge3$ID)

  #### EDGER
  subset = "comp"
  dir_edger = file.path(dir_design, subset, subg, quant, "readcounts")
  file_edger = file.path(dir_edger, "edger_genes.xls")
  tab_edger = read.table(file_edger, header=T, stringsAsFactors=F)
  tab_edger = tab_edger[which(tab_edger$PValue<=0.05), ]
  tab_edger = tab_edger[order(abs(tab_edger$logFC), decreasing=TRUE),]
  tab_edger$ID = unlist(lapply(tab_edger$ID, function(x) unlist(strsplit(x, "[.]"))[1]))
  tab_edger$order = 1:length(tab_edger$ID)

  #subset = "test_train"
  #dir_ref_edger = file.path(dir_design, subset, subg, quant, "readcounts")
  #file_ref_edger = file.path(dir_ref_edger, "edger_genes.xls")
  #tab_ref_edger = read.table(file_ref_edger, header=T, stringsAsFactors=F)
  #edger_merge = merge(tab_ref_edger, tab_edger, by="ID")
  #tab_ref_edger = tab_ref_edger[which(tab_ref_edger$PValue<=0.05), ]
  #tab_ref_edger = tab_ref_edger[order(abs(tab_ref_edger$logFC), decreasing=TRUE),]
  #tab_ref_edger$ID = unlist(lapply(tab_ref_edger$ID, function(x) unlist(strsplit(x, "[.]"))[1]))
  #tab_ref_edger$order = 1:length(tab_ref_edger$ID)
  #edger_merge2 = merge(tab_ref_edger, tab_edger, by="ID", all=F)
  #edger_merge3 = merge(tab_ref_edger, tab_edger, by="ID", all=T)
  #length(edger_merge2$ID) / length(edger_merge3$ID)

  top=50
  tab_ref_keep = tab_ref_epcy[1:top,]
  res_top = data.frame(method="epcy", order=tab_epcy$order[which(tab_epcy$id %in% tab_ref_keep$id)], stringsAsFactors=F)
  res_top = rbind(res_top, data.frame(method="limma", order=tab_limma$order[which(tab_limma$ID %in% tab_ref_keep$id)], stringsAsFactors=F))
  res_top = rbind(res_top, data.frame(method="deseq2", order=tab_deseq2$order[which(tab_deseq2$ID %in% tab_ref_keep$id)], stringsAsFactors=F))
  res_top = rbind(res_top, data.frame(method="edger", order=tab_edger$order[which(tab_edger$ID %in% tab_ref_keep$id)], stringsAsFactors=F))

  res_perf = NULL
  for (top in c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 30, 50, 100)) {
    res_perf = rbind(res_perf, data.frame(method="epcy", topX_gene=top, mean_mcc=sum(tab_ref_epcy$kernel_mcc[which(tab_ref_epcy$id %in% tab_epcy$id[1:top])])/top, stringsAsFactors=F))
    res_perf = rbind(res_perf, data.frame(method="limma", topX_gene=top, mean_mcc=sum(tab_ref_epcy$kernel_mcc[which(tab_ref_epcy$id %in% tab_limma$ID[1:top])])/top, stringsAsFactors=F))
    res_perf = rbind(res_perf, data.frame(method="deseq2", topX_gene=top, mean_mcc=sum(tab_ref_epcy$kernel_mcc[which(tab_ref_epcy$id %in% tab_deseq2$ID[1:top])])/top, stringsAsFactors=F))
    res_perf = rbind(res_perf, data.frame(method="edger", topX_gene=top, mean_mcc=sum(tab_ref_epcy$kernel_mcc[which(tab_ref_epcy$id %in% tab_edger$ID[1:top])])/top, stringsAsFactors=F))
    res_perf = rbind(res_perf, data.frame(method="max_mean_mcc", topX_gene=top, mean_mcc=sum(tab_ref_epcy$kernel_mcc[1:top])/top, stringsAsFactors=F))
  }

  return(list(res_top, res_perf))
}




quant = "STAR"
epcy_folder = "readcounts"
res_top = NULL
res_perf = NULL
for (subg in c("inv16", "MLL", "t8_21", "t15_17", "CK", "NK")) {
  for (replicat in c(1:20)) {
    dir_design = file.path(script_dir, "data", "design", "leucegene_comp_tt_v2", replicat)
    res = run_comp(subg, quant, dir_design, epcy_folder)
    tmp_top = as.data.frame(res[[1]])
    tmp_perf = as.data.frame(res[[2]])

    tmp_top$subg = subg
    tmp_top$replicat = replicat

    tmp_perf$subg = subg
    tmp_perf$replicat = replicat

    res_top = rbind(res_top, tmp_top)
    res_perf = rbind(res_perf, tmp_perf)
  }
}



#gg = ggplot(res_top, aes(x=method, y=log10(order))) + geom_boxplot() #+ geom_jitter()
gg = ggplot(res_top, aes(x=method, y=order)) + geom_boxplot() + geom_jitter()
gg = gg + facet_grid(subg ~ replicat)
gg

png(file="/u/eaudemard/project/epcy_paper/data/fig/res_fig1.png", width = 620, height = 480)
gg
dev.off()


gg = ggplot(res_perf, aes(x=topX_gene, y=mean_mcc, color=method)) + geom_point() + geom_line()
gg = gg + facet_grid(subg ~ replicat)
png(file="/u/eaudemard/project/epcy_paper/data/fig/res_comp_detail.png", width = 620, height = 480)
gg
dev.off()

library(dplyr)

res_perf_mean <- res_perf %>%
    group_by(method, subg, topX_gene) %>%
    summarise(n = n(),
           mean = mean(mean_mcc),
         median = median(mean_mcc),
             sd = sd(mean_mcc)) %>%
    mutate(sem = sd / sqrt(n - 1),
      CI_lower = mean + qt((1-0.8)/2, n - 1) * sem,
      CI_upper = mean - qt((1-0.8)/2, n - 1) * sem)

gg = ggplot(res_perf_mean, aes(x=topX_gene, y=mean, color = method))
gg = gg + geom_line(aes(x=topX_gene, y=mean, color=method))
gg = gg + geom_ribbon(aes(ymin=CI_lower,ymax=CI_upper,fill=method),alpha=0.4)
gg = gg + facet_grid(subg ~ .)
png(file="/u/eaudemard/project/epcy_paper/data/fig/res_comp.png", width=480, height=620)
gg
dev.off()


tab_keep = tab_ref_epcy$id[1:top]
tab_keep = unique(c(tab_keep, tab_epcy$id[1:top]))

tab_ref_keep = tab_ref_epcy[which(tab_ref_epcy$id %in% tab_keep), ]
tab_epcy_keep = tab_epcy[which(tab_epcy$id %in% tab_keep), ]

tab_merge_keep = merge(tab_epcy_keep, tab_ref_keep, by="id", all=T)
gg = ggplot(tab_merge_keep, aes(x=kernel_mcc.x, y=kernel_mcc.y)) + geom_point()
gg


gg = ggplot(epcy_merge, aes(x=l2fc.x, y=l2fc.y)) + geom_point()
gg = gg + ggtitle(paste("pearson: ", cor(epcy_merge$l2fc.x, epcy_merge$l2fc.y), sep=" "))
gg
#gg = ggplot(epcy_merge[1:10,], aes(x=kernel_mcc.x, y=kernel_mcc.y)) + geom_point()
gg = ggplot(epcy_merge, aes(x=kernel_mcc.x, y=kernel_mcc.y)) + geom_point()
gg = gg + ggtitle(paste("pearson: ", cor(epcy_merge$kernel_mcc.x, epcy_merge$kernel_mcc.y), sep=" "))
gg

gg = ggplot(limma_merge, aes(x=logFC.x, y=logFC.y)) + geom_point()
gg = gg + ggtitle(paste("pearson: ", cor(limma_merge$logFC.x, limma_merge$logFC.y), sep=" "))
gg
gg = ggplot(limma_merge, aes(x=-log10(P.Value.x), y=-log10(P.Value.y))) + geom_point()
gg = gg + ggtitle(paste("pearson: ", cor(-log10(limma_merge$P.Value.x), -log10(limma_merge$P.Value.y)), sep=" "))
gg

gg = ggplot(deseq2_merge, aes(x=log2FoldChange.x, y=log2FoldChange.y)) + geom_point()
gg = gg + ggtitle(paste("pearson: ", cor(deseq2_merge$log2FoldChange.x, deseq2_merge$log2FoldChange.y), sep=" "))
gg
gg = ggplot(deseq2_merge, aes(x=-log10(pvalue.x), y=-log10(pvalue.y))) + geom_point()
gg = gg + ggtitle(paste("pearson: ", cor(-log10(deseq2_merge$pvalue.x), -log10(deseq2_merge$pvalue.y)), sep=" "))
gg

gg = ggplot(edger_merge, aes(x=logFC.x, y=logFC.y)) + geom_point()
gg = gg + ggtitle(paste("pearson: ", cor(edger_merge$logFC.x, edger_merge$logFC.y), sep=" "))
gg
gg = ggplot(edger_merge, aes(x=-log10(PValue.x), y=-log10(PValue.y))) + geom_point()
gg = gg + ggtitle(paste("pearson: ", cor(-log10(edger_merge$PValue.x), -log10(edger_merge$PValue.y)), sep=" "))
gg
