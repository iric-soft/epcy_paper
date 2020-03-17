.libPaths("/u/gendrop/R/x86_64-pc-linux-gnu-library/3.6.2")
library("powsimR")

source("./src/script/other/update_powsimr2.r")


data("SmartSeq2_Gene_Read_Counts")
Batches = data.frame(Batch = sapply(strsplit(colnames(SmartSeq2_Gene_Read_Counts), "_"), "[[", 1),
                     stringsAsFactors = F,
                     row.names = colnames(SmartSeq2_Gene_Read_Counts))
data("GeneLengths_mm10")
estparam_gene <- estimateParam(countData = SmartSeq2_Gene_Read_Counts,
                               readData = NULL,
                              batchData = Batches,
                               spikeData = NULL, spikeInfo = NULL,
                               Lengths = GeneLengths_mm10, MeanFragLengths = NULL,
                               RNAseq = 'singlecell', Protocol = 'Read',
                               Distribution = 'ZINB', Normalisation = "scran",
                               GeneFilter = 0.1, SampleFilter = 3,
                               sigma = 1.96, NCores = NULL, verbose = TRUE)
# estimate spike parameters
data("SmartSeq2_SpikeIns_Read_Counts")
data("SmartSeq2_SpikeInfo")
Batches = data.frame(Batch = sapply(strsplit(colnames(SmartSeq2_SpikeIns_Read_Counts), "_"), "[[", 1),
                     stringsAsFactors = F,
                     row.names = colnames(SmartSeq2_SpikeIns_Read_Counts))
estparam_spike <- estimateSpike(spikeData = SmartSeq2_SpikeIns_Read_Counts,
                                spikeInfo = SmartSeq2_SpikeInfo,
                                MeanFragLength = NULL,
                                batchData = Batches,
                                Normalisation = 'depth')
# define log2 fold change
p.lfc <- function(x) sample(c(-1,1), size=x,replace=T)*rgamma(x, shape = 1, rate = 2)
# set up simulations
setupres <- Setup(ngenes = 10000, nsims = 10,
                  p.DE = 0.1, pLFC = p.lfc,
                  n1 = c(20,50,100), n2 = c(30,60,120),
                  Thinning = c(1,0.9,0.8), LibSize = 'given',
                  estParamRes = estparam_gene,
                  estSpikeRes = estparam_spike,
                  DropGenes = FALSE,
                  sim.seed = 52679, verbose = TRUE)
# run simulation
simres <- simulateDE(SetupRes = setupres,
                     Prefilter = "FreqFilter", Imputation = NULL,
                     Normalisation = 'scran', Label = 'none',
                     DEmethod = "limma-trend", DEFilter = FALSE,
                     NCores = NULL, verbose = TRUE,
                     path_count="./data/powsimr/test")
