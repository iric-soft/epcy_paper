
devtools::install_github("ericloud/powsimR", build_vignettes = FALSE, dependencies = FALSE)

setwd("~/dev/epcy_paper/")

library("powsimR")
#source("./src/script/other/update_powsimr_epcy.r")

InputData <- readRDS(file = "/Users/ricou/dev/scRNA-seq-pipelines/simulation/SCRBseq.rds")
estparam.gene <- estimateParam(countData = InputData$GeneUMIs, 
                               readData = InputData$GeneReads, 
                               Lengths = InputData$GeneLengths,
                               RNAseq = "singlecell", Protocol = "UMI", Distribution = "NB", 
                               Normalisation = "scran", 
                               GeneFilter = 0.1, SampleFilter = 3)

estparam.spike <- estimateSpike(spikeData = InputData$SpikeCounts, spikeInfo = InputData$SpikeInfo, 
                                Normalisation = "depth")
plotParam(estparam.gene)

ngenes = 10000
n1 = 50
n2 = 200

lfc.gamma.narrow.asym = function(x) sample(c(-1,1), size=x, prob = c(0.25,0.75),replace=T)*rgamma(x, shape = 1, rate = 2)
setup <- Setup(ngenes = ngenes, nsims = 1, 
               p.DE = 0.2, pLFC = lfc.gamma.narrow.asym,
               n1 = c(n1), n2 = c(n2), 
               LibSize = 'equal', DropGenes = TRUE, 
               estParamRes = estparam.gene, estSpikeRes = estparam.spike)


simres <- simulateDE(SetupRes = setup,
                     Normalisation = "scran", Label = "clustering", 
                     DEmethod = "limma-trend",
                     path_count = "./data/simu/",
                     Counts=T,
                     CountsNorm=T)
df = as.data.frame(simres$CountsNorm$`50vs200`[[1]])
df <- data.frame(id=rownames(df), data.frame(df, row.names=NULL))
write.table(df, "./data/simu/simu1/50vs200/norm/readcounts.tsv", row.names = FALSE, sep="\t", quote=FALSE)
df_var = data.frame(id = df$id, var_query=rowVars(log2(as.matrix(df[,2:51])+1)), var_ref=rowVars(log2(as.matrix(df[,52:251])+1)), stringsAsFactors = F)


df = as.data.frame(simres$Counts$`50vs200`[[1]])
df <- data.frame(id=rownames(df), data.frame(df, row.names=NULL))
write.table(df, "./data/simu/simu1/50vs200/raw/readcounts.tsv", row.names = FALSE, sep="\t", quote=FALSE)

eval.de <- evaluateDE(simRes = simres, alpha.type = "adjusted", MTC = "BY")
eval.de$TP.marginal
eval.de$TN.marginal
eval.de$FN.marginal
eval.de$FP.marginal
plotEvalDE(eval.de, rate = "marginal")


pval = simres$SimulateRes$pvalue[,1,1]
pval = pval[which(!is.na(pval))]
meanexpr = log1p(simres$SimulateRes$mu[,1,1])
padj = stats::p.adjust(pval, method = "BY")
padj[is.na(padj)] = 1

tab_epcy = read.table("./data/simu/simu1/50vs200/norm/predictive_capability.xls", header=TRUE, sep="\t", stringsAsFactors=F)
tab_epcy$kernel_mcc[is.na(tab_epcy$kernel_mcc)] = 0
tab_epcy$diff = FALSE

ids_de = rep(0, ngenes)
ids_de[simres$DESetup$DEid[1][[1]]] = 1
ids_de = ids_de[which(!is.na(simres$SimulateRes$mu))]
tab_epcy$diff[which(ids_de == 1)] = TRUE
tab_epcy$padj = padj
tab_epcy$trend_padj = FALSE
tab_epcy$trend_padj[which(padj <= 0.1 )] = TRUE
tab_epcy$trend_pvalue = FALSE
tab_epcy$trend_pvalue[which(pval <= 0.05 )] = TRUE

library(ggplot2)
breaks = c(1.5,0.1,-0.1, -1.5)
gg = ggplot(tab_epcy, aes(x=diff, y=kernel_mcc, shape=trend_padj, color=-log10(padj)))
gg = gg + geom_point(alpha=0.5, position=position_jitter(width=0.3))
gg = gg + geom_hline(yintercept=0.1, linetype="dashed", color="black")
gg = gg + scale_colour_gradientn(colours=rainbow(4))
gg


sum(tab_epcy$kernel_mcc>=0.1 & tab_epcy$diff==T)
sum(tab_epcy$kernel_mcc<0.1 & tab_epcy$diff==F)
sum(tab_epcy$kernel_mcc>=0.1 & tab_epcy$diff==F)
sum(tab_epcy$kernel_mcc<0.1 & tab_epcy$diff==T)


sum(tab_epcy$padj<=0.1 & tab_epcy$diff==T)
sum(tab_epcy$padj>0.1 & tab_epcy$diff==F)
sum(tab_epcy$padj<=0.1 & tab_epcy$diff==F)
sum(tab_epcy$padj>0.1 & tab_epcy$diff==T)

#ENSMUSG00000107096.3_81 ENSMUSG00000074415.12_1835 ENSMUSG00000056486.17_2289 ENSMUSG00000017405.14_2641 ENSMUSG00000113450.1_3346 ENSMUSG00000052727.6_3952 ENSMUSG00000081440.2_9067 ENSMUSG00000038280.12_9622
tab_epcy$id[which(tab_epcy$padj<=0.1 & tab_epcy$diff==F)]
#ENSMUSG00000055341.10_2 ENSMUSG00000037110.19_63 ENSMUSG00000048174.2_129
tab_epcy$id[which(tab_epcy$padj<=0.1 & tab_epcy$diff==T & tab_epcy$kernel_mcc<0.1)][1:3]

#ENSMUSG00000066613.14_229 ENSMUSG00000039530.17_1036 ENSMUSG00000021127.7_8438
tab_epcy$id[which(-log10(tab_epcy$padj)<60 & tab_epcy$diff==T & tab_epcy$kernel_mcc>0.95)]

tab_epcy$id[which(-log10(tab_epcy$padj)>150 & tab_epcy$diff==T & tab_epcy$kernel_mcc>0.8 & tab_epcy$kernel_mcc < 0.95)]


tab_epcy_var = merge(tab_epcy, df_var, by="id")
tab_epcy_var$var = pmax(tab_epcy_var$var_query, tab_epcy_var$var_ref)

num_diff_epcy = sum(tab_epcy_var$kernel_mcc>= 0.1)
tab_epcy_var$padj[order(tab_epcy_var$padj, decreasing = F)][num_diff_epcy]


gg = ggplot(tab_epcy_var, aes(x=diff, y=kernel_mcc, shape=trend_padj, color=var))
gg = gg + geom_point(alpha=0.5, position=position_jitter(width=0.3))
gg = gg + geom_hline(yintercept=0.1, linetype="dashed", color="black")
gg = gg + scale_colour_gradientn(colours=rainbow(4))
gg

gg = ggplot(tab_epcy_var, aes(x=diff, y=-log10(padj), shape=kernel_mcc>=0.1, color=var))
gg = gg + geom_point(alpha=0.5, position=position_jitter(width=0.3))
gg = gg + geom_hline(yintercept=0.1, linetype="dotted", color="black")
gg = gg + geom_hline(yintercept=-log10(tab_epcy_var$padj[order(tab_epcy_var$padj, decreasing = F)][num_diff_epcy]), linetype="dashed", color="black")
gg = gg + scale_colour_gradientn(colours=rainbow(4))
gg


gg = ggplot(tab_epcy_var, aes(x=var, y=kernel_mcc, shape=trend_padj, color=-log10(padj)))
gg = gg + geom_point(alpha=0.5, position=position_jitter(width=0.3))
gg = gg + scale_colour_gradientn(colours=rainbow(4))
gg

