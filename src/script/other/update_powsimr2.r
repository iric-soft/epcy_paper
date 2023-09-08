simulateDE <- function (SetupRes,
                        Prefilter = NULL,
                        Imputation = NULL,
                        Normalisation = c("TMM", "MR", "PosCounts", "UQ",
                                          "scran", "Linnorm", "sctransform",
                                          "SCnorm", "Census", "depth"),
                        Label = "none",
                        DEmethod = c("T-Test", "edgeR-LRT", "edgeR-QL",
                                     "edgeR-zingeR", "edgeR-ZINB-WaVE",
                                     "limma-voom", "limma-trend",
                                     "DESeq2", "DESeq2-zingeR", "DESeq2-ZINB-WaVE",
                                     "ROTS", "baySeq", "NOISeq", "EBSeq",
                                     "MAST", "BPSC", "scDD", "DECENT"),
                        DEFilter = FALSE,
                        NCores = NULL,
                        verbose = TRUE,
                        path_count="")
{
    if (!is.null(NCores) && DEmethod %in% c("edgeR-LRT", "edgeR-QL",
        "edgeR-zingeR", "DESeq2-zingeR", "limma-voom", "limma-trend",
        "NOISeq", "EBSeq", "ROTS")) {
        if (verbose) {
            message(paste0(DEmethod, " has no parallel computation option!"))
        }
    }
    if (attr(SetupRes, "RNAseq") == "singlecell" && DEmethod %in%
        c("edgeR-LRT", "edgeR-QL", "limma-voom", "limma-trend",
            "DESeq2", "baySeq", "NOISeq", "EBSeq")) {
        if (verbose) {
            message(paste0(DEmethod, " is developed for bulk RNA-seq experiments."))
        }
    }
    if (attr(SetupRes, "RNAseq") == "bulk" && DEmethod %in% c("MAST",
        "BPSC", "edgeR-zingeR", "DESeq2-zingeR", "edgeR-ZINB-WaVE",
        "DESeq2-ZINB-WaVE")) {
        if (verbose) {
            message(paste0(DEmethod, " is developed for single cell RNA-seq experiments."))
        }
    }
    if (attr(SetupRes, "RNAseq") == "bulk" && DEmethod %in% c("scDD",
        "DECENT")) {
        stop(message(paste0(DEmethod, " is only developed and implemented for single cell RNA-seq experiments.")))
    }
    if (DEmethod == "DECENT") {
        if (verbose) {
            message(paste0(DEmethod, " does not require additional normalisation nor imputation."))
        }
        Normalisation = "none"
        Prefilter = NULL
        Imputation = NULL
    }
    if (c(all(is.null(Imputation), is.null(Prefilter)) && isTRUE(DEFilter))) {
        stop(message(paste0("You wish to use imputed/filtered gene expression values for DE testing but you did not specify the imputation/filtering method. Aborting.")))
    }
    if (!is.null(Imputation) && attr(SetupRes, "RNAseq") == "bulk") {
        message(paste0("You wish to use imputation but powsimR has only methods implemented for single cell RNA-seq and in most cases imputation is not needed for bulk RNA-seq. Setting Imputation to NULL."))
        Imputation = NULL
    }
    start.time.pipe = proc.time()
    max.n = max(SetupRes$SimSetup$n1, SetupRes$SimSetup$n2)
    min.n = min(SetupRes$SimSetup$n1, SetupRes$SimSetup$n2)
    if (Label == "clustering") {
        PreclustNumber <- min.n
    }
    if (!Label == "clustering") {
        PreclustNumber <- NULL
    }
    Pipeline <- list(Prefilter = Prefilter, Imputation = Imputation,
        Normalisation = Normalisation, Label = Label, DEmethod = DEmethod,
        DEFilter = DEFilter, NCores = NCores, clustNumber = ifelse(SetupRes$DESetup$design ==
            "2grp", 2, NULL), PreclustNumber = PreclustNumber)
    SetupRes <- c(SetupRes, Pipeline = list(Pipeline))
    if (verbose) {
        message(paste0("Preparing output arrays."))
    }
    my.names = paste0(SetupRes$SimSetup$n1, "vs", SetupRes$SimSetup$n2)
    pvalues = fdrs = elfcs = rlfcs = mus = true.mus = disps = true.disps = dropouts = true.drops = array(NA,
        dim = c(SetupRes$DESetup$ngenes, length(SetupRes$SimSetup$n1),
            SetupRes$DESetup$nsims))
    tstep <- c("Simulation", "Preprocessing", "Normalisation",
        "DE", "Moments", "Total")
    tmoment <- c("User", "System", "Elapsed")
    tname <- paste(rep(tstep, each = 3), rep(tmoment, times = 6),
        sep = "_")
    true.sf <- lapply(1:length(my.names), function(x) {
        matrix(NA, nrow = SetupRes$DESetup$nsims, ncol = SetupRes$SimSetup$n1[x] +
            SetupRes$SimSetup$n2[x])
    })
    est.sf <- lapply(1:length(my.names), function(x) {
        matrix(NA, nrow = SetupRes$DESetup$nsims, ncol = SetupRes$SimSetup$n1[x] +
            SetupRes$SimSetup$n2[x])
    })
    true.depth <- lapply(1:length(my.names), function(x) {
        matrix(NA, nrow = SetupRes$DESetup$nsims, ncol = SetupRes$SimSetup$n1[x] +
            SetupRes$SimSetup$n2[x])
    })
    est.depth <- lapply(1:length(my.names), function(x) {
        matrix(NA, nrow = SetupRes$DESetup$nsims, ncol = SetupRes$SimSetup$n1[x] +
            SetupRes$SimSetup$n2[x])
    })
    time.taken <- lapply(1:length(my.names), function(x) {
        matrix(NA, nrow = SetupRes$DESetup$nsims, ncol = length(tname),
            dimnames = list(NULL, tname))
    })
    true.designs <- lapply(1:length(my.names), function(x) {
        matrix(NA, nrow = SetupRes$DESetup$nsims, ncol = SetupRes$SimSetup$n1[x] +
            SetupRes$SimSetup$n2[x])
    })
    names(true.sf) = names(est.sf) = names(true.depth) = names(est.depth) = names(time.taken) = names(true.designs) = my.names
    if (SetupRes$Pipeline$Normalisation %in% c("SCnorm", "Linnorm",
        "sctransform", "bayNorm")) {
        est.gsf <- lapply(1:length(my.names), function(x) {
            array(NA, dim = c(SetupRes$DESetup$nsims, SetupRes$DESetup$ngenes,
                SetupRes$SimSetup$n1[x] + SetupRes$SimSetup$n2[x]))
        })
        names(est.gsf) = my.names
    }
    if (!SetupRes$Pipeline$Normalisation %in% c("SCnorm", "Linnorm",
        "sctransform", "bayNorm")) {
        est.gsf = NULL
    }
    for (i in 1:SetupRes$DESetup$nsims) {
        if (verbose) {
            message(paste0("\n  SIMULATION   NUMBER   ", i, "\n"))
        }
        start.time.sim1 <- proc.time()
        tmp.simOpts = SetupRes
        tmp.simOpts$DESetup$DEid = SetupRes$DESetup$DEid[[i]]
        tmp.simOpts$DESetup$pLFC = SetupRes$DESetup$pLFC[[i]]
        tmp.simOpts$DESetup$Bid = SetupRes$DESetup$Bid[[i]]
        tmp.simOpts$DESetup$bLFC = SetupRes$DESetup$bLFC[[i]]
        tmp.simOpts$DESetup$sim.seed = SetupRes$DESetup$sim.seed[[i]]
        if (verbose) {
            message(paste0("Generating gene expression."))
        }
        gene.data = powsimR::.simRNAseq.2grp(simOptions = tmp.simOpts,
            n1 = max.n, n2 = max.n, verbose = verbose)
        if (isTRUE(tmp.simOpts$SimSetup$DropGenes)) {
            gene.data = .dropGene(simOptions = tmp.simOpts, simData = gene.data)
        }
        if (isTRUE(tmp.simOpts$SimSetup$spikeIns)) {
            if (verbose) {
                message(paste0("Generating spike-in expression."))
            }
            spike.data = .simSpike(SpikeOptions = tmp.simOpts$estSpikeRes,
                n1 = max.n, n2 = max.n, sf = gene.data$sf)
            spike.info = tmp.simOpts$estSpikeRes$FilteredInput$spikeInfo
        }
        if (!isTRUE(tmp.simOpts$SimSetup$spikeIns)) {
            spike.data = NULL
            spike.info = NULL
        }
        if (!is.null(tmp.simOpts$estParamRes$MeanFragLengths)) {
            if (verbose) {
                message(paste0("Sampling from observed mean fragment lengths"))
            }
            MeanFrag.data = sample(tmp.simOpts$estParamRes$MeanFragLengths,
                max.n + max.n, replace = TRUE)
            names(MeanFrag.data) = colnames(gene.data$counts)
        }
        if (is.null(tmp.simOpts$estParamRes$MeanFragLengths)) {
            MeanFrag.data = NULL
        }
        if (!is.null(tmp.simOpts$estParamRes$Lengths)) {
            gene.id = sub("_([^_]*)$", "", rownames(gene.data$counts))
            Length.data = tmp.simOpts$estParamRes$Lengths
            Length.data = Length.data[match(gene.id, names(Length.data))]
        }
        if (is.null(tmp.simOpts$estParamRes$Lengths)) {
            Length.data = NULL
        }
        end.time.sim1 <- proc.time()
        for (j in seq(along = tmp.simOpts$SimSetup$n1)) {
            start.time.sim2 <- proc.time()
            Nrep1 = tmp.simOpts$SimSetup$n1[j]
            Nrep2 = tmp.simOpts$SimSetup$n2[j]
            Thin = tmp.simOpts$SimSetup$Thinning[j]
            if (verbose) {
                message(paste0(Nrep1, " vs. ", Nrep2))
            }
            idx = c(1:Nrep1, max.n + (1:Nrep2))
            true.design = gene.data$designs[idx]
            sim.cnts = gene.data$counts[, idx]
            gene.sf = gene.data$sf[idx]
            if (!is.null(Thin) && !Thin == 1) {
                if (verbose) {
                  message(paste0("Reduce the size of gene counts to ",
                    Thin, " using binomial thinning."))
                }
                sim.cnts = .run.thin(countData = sim.cnts, Thin = Thin,
                  simOptions = tmp.simOpts)
            }
            ix.valid = rowSums(sim.cnts) > 0
            count.data = sim.cnts[ix.valid, , drop = FALSE]
            sim.depth = colSums(count.data)
            if (!is.null(Length.data)) {
                if (verbose) {
                  message(paste0("Associating gene lengths with gene expression"))
                }
                gene.id = sub("_([^_]*)$", "", rownames(count.data))
                length.data = Length.data
                length.data = length.data[match(gene.id, names(length.data))]
            }
            if (is.null(Length.data)) {
                length.data = NULL
            }
            if (!is.null(spike.data) && !is.null(spike.info)) {
                sim.spike <- spike.data$counts
                spike.valid = rowSums(sim.spike) > 0
                count.spike = sim.spike[spike.valid, idx, drop = FALSE]
                info.spike <- tmp.simOpts$estSpikeRes$FilteredInput$spikeInfo
                info.spike <- info.spike[rownames(info.spike) %in%
                  rownames(count.spike), , drop = FALSE]
                info.spike <- info.spike[match(rownames(count.spike),
                  rownames(info.spike)), , drop = FALSE]
            }
            if (c(!is.null(spike.data) && !is.null(spike.info) &&
                !is.null(Thin) && !Thin == 1 && isTRUE(tmp.simOpts$SimSetup$thinSpike))) {
                if (verbose) {
                  message(paste0("Reduce the size of spike-in counts to ",
                    Thin, " using binomial thinning."))
                }
                count.spike = .run.thin(countData = count.spike,
                  Thin = Thin, simOptions = tmp.simOpts)
                count.spike = count.spike[rowSums(count.spike) >
                  0, ]
                info.spike <- tmp.simOpts$estSpikeRes$FilteredInput$spikeInfo
                info.spike <- info.spike[rownames(info.spike) %in%
                  rownames(count.spike), , drop = FALSE]
                info.spike <- info.spike[match(rownames(count.spike),
                  rownames(info.spike)), , drop = FALSE]
            }
            if (is.null(spike.data) && is.null(spike.info)) {
                count.spike = NULL
                info.spike = NULL
            }
            if (!is.null(MeanFrag.data)) {
                meanfrag.data = MeanFrag.data[idx]
            }
            if (is.null(MeanFrag.data)) {
                meanfrag.data = NULL
            }
            def.design <- true.design
            end.time.sim2 <- proc.time()
            start.time.preprocess <- proc.time()
            if (!is.null(tmp.simOpts$Pipeline$Prefilter)) {
                if (verbose) {
                  message(paste0("Applying ", tmp.simOpts$Pipeline$Prefilter,
                    " prefiltering"))
                }
                filter.data <- .prefilter.calc(Prefilter = tmp.simOpts$Pipeline$Prefilter,
                  countData = count.data, NCores = tmp.simOpts$Pipeline$NCores)
                filter.count.data <- filter.data
                if (!is.null(Length.data)) {
                  gene.id = sub("_([^_]*)$", "", rownames(filter.count.data))
                  length.data = Length.data
                  length.data = length.data[match(gene.id, names(length.data))]
                }
                if (is.null(Length.data)) {
                  length.data = NULL
                }
            }
            if (is.null(tmp.simOpts$Pipeline$Prefilter)) {
                filter.count.data <- count.data
            }
            if (!is.null(tmp.simOpts$Pipeline$Imputation)) {
                if (verbose) {
                  message(paste0("Applying ", tmp.simOpts$Pipeline$Imputation,
                    " imputation"))
                }
                impute.data <- .impute.calc(Imputation = tmp.simOpts$Pipeline$Imputation,
                  countData = filter.count.data, spikeData = count.spike,
                  batchData = def.design, clustNumber = tmp.simOpts$Pipeline$clustNumber,
                  Lengths = length.data, MeanFragLengths = meanfrag.data,
                  NCores = tmp.simOpts$Pipeline$NCores, verbose = verbose)
                fornorm.count.data <- impute.data
                if (!is.null(Length.data)) {
                  gene.id = sub("_([^_]*)$", "", rownames(fornorm.count.data))
                  length.data = Length.data
                  length.data = length.data[match(gene.id, names(length.data))]
                }
                if (is.null(Length.data)) {
                  length.data = NULL
                }
            }
            if (is.null(Imputation)) {
                fornorm.count.data <- filter.count.data
            }
            end.time.preprocess <- proc.time()
            start.time.norm <- proc.time()
            if (verbose) {
                message(paste0("Applying ", tmp.simOpts$Pipeline$Normalisation,
                  " normalisation"))
            }
            if (tmp.simOpts$Pipeline$Normalisation == "sctransform") {
                def.design <- NULL
            }
            norm.data <- .norm.calc(Normalisation = tmp.simOpts$Pipeline$Normalisation,
                sf = gene.sf, countData = fornorm.count.data,
                spikeData = count.spike, spikeInfo = info.spike,
                batchData = def.design, Lengths = length.data,
                MeanFragLengths = meanfrag.data, PreclustNumber = tmp.simOpts$Pipeline$PreclustNumber,
                Step = "Simulation", Protocol = attr(tmp.simOpts$estParamRes,
                  "Protocol"), Label = tmp.simOpts$Pipeline$Label,
                NCores = tmp.simOpts$Pipeline$NCores, verbose = verbose)
            end.time.norm <- proc.time()
            DEOpts <- list(designs = def.design, p.DE = tmp.simOpts$DESetup$p.DE)
            if (!is.null(Length.data)) {
                if (verbose) {
                  message(paste0("Reassociate gene lengths with gene expression"))
                }
                gene.id = sub("_([^_]*)$", "", rownames(count.data))
                length.data = Length.data
                length.data = length.data[match(gene.id, names(length.data))]
            }
            if (is.null(Length.data)) {
                length.data = NULL
            }

            ## save simulate readcount
            if(path_count != "") {
              num_sample = paste0(Nrep1, " vs. ", Nrep2)
              num_simu = paste0("simu", i)
              path_out = file.path(path_count, num_simu, num_sample)
              file_out = file.path(path_out, "readcounts.tsv")
              dir.create(path_out, recursive = TRUE, showWarnings = FALSE)
              write.table(count.data, file_out, row.names = FALSE)
            }

            start.time.DE <- proc.time()
            if (!isTRUE(tmp.simOpts$Pipeline$DEFilter)) {
                if (verbose) {
                  message(paste0("Applying ", tmp.simOpts$Pipeline$DEmethod,
                    " for DE analysis on raw count data."))
                }
                res.de = .de.calc(DEmethod = tmp.simOpts$Pipeline$DEmethod,
                  normData = norm.data, countData = count.data,
                  Lengths = length.data, MeanFragLengths = meanfrag.data,
                  DEOpts = DEOpts, spikeData = count.spike, spikeInfo = info.spike,
                  NCores = tmp.simOpts$Pipeline$NCores, verbose = verbose)
            }
            if (isTRUE(tmp.simOpts$Pipeline$DEFilter)) {
                if (verbose) {
                  message(paste0("Applying ", tmp.simOpts$Pipeline$DEmethod,
                    " for DE analysis on imputed/filtered count data."))
                }
                res.de = .de.calc(DEmethod = tmp.simOpts$Pipeline$DEmethod,
                  normData = norm.data, countData = fornorm.count.data,
                  Lengths = length.data, MeanFragLengths = meanfrag.data,
                  DEOpts = DEOpts, spikeData = count.spike, spikeInfo = info.spike,
                  NCores = tmp.simOpts$Pipeline$NCores, verbose = verbose)
            }
            end.time.DE <- proc.time()
            if (tmp.simOpts$Pipeline$DEmethod == "DECENT") {
                res.de <- res.de[["DEresults"]]
                norm.data <- res.de[["NormData"]]
            }
            if (attr(norm.data, "normFramework") %in% c("SCnorm",
                "Linnorm", "scTransform")) {
                allgenes <- rownames(sim.cnts)
                estgenes <- rownames(norm.data$scale.factors)
                ixx.valid <- allgenes %in% estgenes
                est.gsf[[j]][i, ixx.valid, ] = norm.data$scale.factors
                norm.data$scale.factors = est.gsf[[j]][i, , ]
            }
            start.time.moments <- proc.time()
            if (isTRUE(tmp.simOpts$Pipeline$DEFilter)) {
                if (verbose) {
                  message(paste0("Estimating moments of imputed/filtered count data."))
                }
                res.params <- .run.params(countData = fornorm.count.data,
                  normData = norm.data, group = DEOpts$designs)
                est.depths <- colSums(fornorm.count.data)
            }
            if (!isTRUE(tmp.simOpts$Pipeline$DEFilter)) {
                if (verbose) {
                  message(paste0("Estimating moments of raw count data."))
                }
                res.params <- .run.params(countData = count.data,
                  normData = norm.data, group = DEOpts$designs)
                est.depths <- colSums(count.data)
            }
            end.time.moments <- proc.time()
            pval = fdr = est.lfc = raw.lfc = mu.tmp = true.mu.tmp = disp.tmp = true.disp.tmp = p0.tmp = true.p0.tmp = rep(NA,
                nrow(sim.cnts))
            allgenes <- rownames(sim.cnts)
            testedgenes <- res.de$geneIndex
            ixx.de.valid <- allgenes %in% testedgenes
            pval[ixx.de.valid] = res.de$pval
            fdr[ixx.de.valid] = res.de$fdr
            est.lfc[ixx.de.valid] = res.de$lfc
            paramgenes <- res.params$geneIndex
            ixx.param.valid <- allgenes %in% paramgenes
            raw.lfc[ixx.param.valid] = res.params$lfc
            mu.tmp[ixx.param.valid] = res.params$means
            disp.tmp[ixx.param.valid] = res.params$dispersion
            p0.tmp[ixx.param.valid] = res.params$dropout
            true.mu.tmp = gene.data$mus
            true.disp.tmp = gene.data$disps
            true.p0.tmp = gene.data$drops
            pvalues[, j, i] = pval
            fdrs[, j, i] = fdr
            elfcs[, j, i] = est.lfc
            rlfcs[, j, i] = raw.lfc
            mus[, j, i] = mu.tmp
            true.mus[, j, i] = true.mu.tmp
            disps[, j, i] = disp.tmp
            true.disps[, j, i] = true.disp.tmp
            dropouts[, j, i] = p0.tmp
            true.drops[, j, i] = true.p0.tmp
            true.sf[[j]][i, ] = gene.sf
            est.sf[[j]][i, ] = norm.data$size.factors
            true.depth[[j]][i, ] = sim.depth
            est.depth[[j]][i, ] = est.depths
            true.designs[[j]][i, ] = true.design
            end.time.pipe = proc.time()
            time.taken.sim <- (end.time.sim1 - start.time.sim1) +
                (end.time.sim2 - start.time.sim2)
            if (any(c(!is.null(tmp.simOpts$Pipeline$Prefilter),
                !is.null(tmp.simOpts$Pipeline$Imputation)))) {
                time.taken.preprocess <- end.time.preprocess -
                  start.time.preprocess
            }
            if (all(c(is.null(tmp.simOpts$Pipeline$Prefilter),
                is.null(tmp.simOpts$Pipeline$Imputation)))) {
                time.taken.preprocess = c(NA, NA, NA)
            }
            time.taken.norm <- end.time.norm - start.time.norm
            time.taken.DE <- end.time.DE - start.time.DE
            time.taken.moments <- end.time.moments - start.time.moments
            time.taken.total <- end.time.pipe - start.time.pipe
            timing <- c(time.taken.sim, time.taken.preprocess,
                time.taken.norm, time.taken.DE, time.taken.moments,
                time.taken.total)
            timing <- timing[!grepl(pattern = "child", names(timing))]
            time.taken[[j]][i, ] = timing
        }
    }
    Simulate <- list(pvalue = pvalues, fdr = fdrs, elfc = elfcs,
        rlfc = rlfcs, mu = mus, true.mu = true.mus, disp = disps,
        true.disp = true.disps, dropout = dropouts, true.dropout = true.drops,
        true.sf = true.sf, true.depth = true.depth, est.depth = est.depth,
        est.sf = est.sf, est.gsf = est.gsf, true.designs = true.designs,
        time.taken = time.taken)
    attr(Simulate, "Simulation") <- "DE"
    res.out <- c(SetupRes, SimulateRes = list(Simulate))
    return(res.out)
}
