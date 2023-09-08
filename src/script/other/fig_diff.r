library(ggplot2)


## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    library(plyr)

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
      .fun = function(xx, col) {
        c(N    = length2(xx[[col]], na.rm=na.rm),
          mean = mean   (xx[[col]], na.rm=na.rm),
          sd   = sd     (xx[[col]], na.rm=na.rm)
        )
      },
      measurevar
    )

    # Rename the "mean" column
    datac <- rename(datac, c("mean" = measurevar))

    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval:
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
}




df_solo = data.frame(
  quant_value=c(1,1,1,1,3,6),
  class=c("A", "A", "A", "B", "B", "B"),
  label=c("Case A", "Case B", "Case C")
)

gg_fig1 = ggplot(df_solo, aes(x=class, y=quant_value)) + geom_point()
gg_fig1 = gg_fig1 + facet_wrap(~ label, ncol=3)
gg_fig1 = gg_fig1 + theme_bw()
gg_fig1 = gg_fig1 + theme(plot.title = element_text(hjust = 0.5))
gg_fig1 = gg_fig1 + ggtitle("Which case is different ?")
gg_fig1

png(file="/u/eaudemard/project/epcy_paper/data/fig/fig1.png", width = 620, height = 480)
gg_fig1
dev.off()



df_sd = data.frame(
  quant_value = c(
    rnorm(500, 1, 1),
    rnorm(500, 6, 1),
    rnorm(500, 1, 2),
    rnorm(500, 6, 2),
    rnorm(500, 1, 4),
    rnorm(500, 6, 4)
  ),
  class = c(
    rep("A", 500),
    rep("B", 500),
    rep("A", 500),
    rep("B", 500),
    rep("A", 500),
    rep("B", 500)
  ),
  label = c(
    rep("Case A", 1000),
    rep("Case B", 1000),
    rep("Case C", 1000)
  )
)

df_sd_se <- summarySE(df_sd, measurevar="quant_value", groupvars="class")

gg_fig2 = ggplot() + geom_jitter(
                      aes(x=class, y=quant_value, color=I("blue")), data = df_sd,
                      position = position_jitter(width = 0.05))
gg_fig2 = gg_fig2 + geom_crossbar(
                      data=df_sd_se,
                      aes(x=class, ymin=quant_value,
                          ymax=quant_value,y=quant_value,
                          group=class),
                      width = 0.5)
gg_fig2 = gg_fig2 + facet_wrap(~ label, ncol=3)
gg_fig2 = gg_fig2 + theme_bw()
gg_fig2 = gg_fig2 + theme(plot.title = element_text(hjust = 0.5))
gg_fig2 = gg_fig2 + ggtitle("Which case have the most different groups ?")

gg_fig2

png(file="/u/eaudemard/project/epcy_paper/data/fig/fig2.png", width = 620, height = 480)
gg_fig2
dev.off()



df_stat = NULL
df_title = NULL
mean_a = 1
mean_b = 6
for (n in c(5, 10, 100)) {
  for (label in c("Case A", "Case B", "Case C")) {
    if (label == "Case A") {
      a = rnorm(n, mean_a, 1)
      b = rnorm(n, mean_b, 1)
    } else if (label == "Case B") {
      a = rnorm(n, mean_a, 3)
      b = rnorm(n, mean_b, 3)
    } else {
      a = rnorm(n, mean_a, 6)
      b = rnorm(n, mean_b, 6)
    }

    res_stat = signif(wilcox.test(a, b)$p.value,3)
    res_l2fc = signif(log2(mean(b) + 1) - log2(mean(a) + 1), 3)
    df_stat = rbind(df_stat, data.frame(
                               quant_value = c(a,b),
                               class = c(rep("A", n), rep("B", n)),
                               label = rep(label, 2*n),
                               n = rep(n, 2*n)
                             )
                   )
    df_title = rbind(df_title, data.frame(
                                 title = paste("pv=", res_stat, " | fc=", res_l2fc, sep=""),
                                 pvalue = res_stat,
                                 l2fc = res_l2fc,
                                 label = label,
                                 n=n,
                                 class="A",
                                 quant_value=27
                                )
                    )
  }
}

df_stat$n <- factor(df_stat$n, levels = c("5", "10", "100"))
df_title$n <- factor(df_title$n, levels = c("5", "10", "100"))

gg_fig3 = ggplot() + geom_jitter(
                      aes(x=class, y=quant_value, color=I("blue")), data = df_stat,
                      position = position_jitter(width = 0.05))
gg_fig3 = gg_fig3 + geom_text(data=df_title, aes(x=class, y=quant_value, label=title))
gg_fig3 = gg_fig3 + facet_grid(n ~ label)
gg_fig3 = gg_fig3 + theme_bw()
gg_fig3 = gg_fig3 + ylim(c(-20,30))
gg_fig3 = gg_fig3 + theme(plot.title = element_text(hjust = 0.5))
gg_fig3 = gg_fig3 + ggtitle("Which case have the most different groups ?")

png(file="/u/eaudemard/project/epcy_paper/data/fig/fig3.png", width = 680, height = 480)
gg_fig3
dev.off()


gg_fig3b = ggplot(df_title, aes(x=l2fc, y=-log10(pvalue), color=label)) + geom_point()
gg_fig3b = gg_fig3b + facet_grid(n ~ .)
gg_fig3b = gg_fig3b + theme_bw()
gg_fig3b = gg_fig3b + xlim(c(-3,3))
gg_fig3b = gg_fig3b + geom_hline(yintercept=range(1.30103))
gg_fig3b = gg_fig3b + geom_vline(xintercept=range(0))
gg_fig3b = gg_fig3b + theme(plot.title = element_text(hjust = 0.5))
gg_fig3b = gg_fig3b + ggtitle("Volcano plot")
png(file="/u/eaudemard/project/epcy_paper/data/fig/fig3b.png", width = 480, height = 680)
gg_fig3b
dev.off()



df_simu = NULL
for (i in 1:1000) {
  for (n in c(5, 10, 100)) {
    for (label in c("Case A", "Case B", "Case C")) {
      if (label == "Case A") {
        a = rnorm(n, mean_a, 1)
        b = rnorm(n, mean_b, 1)
      } else if (label == "Case B") {
        a = rnorm(n, mean_a, 3)
        b = rnorm(n, mean_b, 3)
      } else {
        a = rnorm(n, mean_a, 6)
        b = rnorm(n, mean_b, 6)
      }

      res_stat = signif(wilcox.test(a, b)$p.value,3)
      res_l2fc = signif(log2(mean(b) + 1) - log2(mean(a) + 1), 3)
      df_simu = rbind(df_simu, data.frame(
                          title = paste("pv=", res_stat, " | fc=", res_l2fc, sep=""),
                          pvalue = res_stat,
                          l2fc = res_l2fc,
                          label = label,
                          n=n,
                          simu=i
                         )
             )
    }
  }
}
df_simu$significant = df_simu$pvalue <= 0.05


gg_fig4 = ggplot(df_simu, aes(x=significant)) + geom_bar(stat="count")
gg_fig4 = gg_fig4 + facet_grid(n ~ label)
gg_fig4 = gg_fig4 + theme_bw()
gg_fig4 = gg_fig4 + theme(plot.title = element_text(hjust = 0.5))
gg_fig4 = gg_fig4 + ggtitle("Significant for 1000 samples of case A, B, C")
png(file="/u/eaudemard/project/epcy_paper/data/fig/fig4.png", width = 680, height = 480)
gg_fig4
dev.off()

gg_fig4b = ggplot(df_simu, aes(x=l2fc)) + geom_bar(aes(fill=significant), width=.5)
gg_fig4b = gg_fig4b + geom_vline(xintercept=range(1.81))
gg_fig4b = gg_fig4b + facet_grid(n ~ label)
gg_fig4b = gg_fig4b + theme_bw()
gg_fig4b = gg_fig4b + theme(plot.title = element_text(hjust = 0.5))
gg_fig4b = gg_fig4b + ggtitle("L2FC for 1000 samples of case A, B, C")
png(file="/u/eaudemard/project/epcy_paper/data/fig/fig4b.png", width = 680, height = 480)
gg_fig4b
dev.off()



a = rnorm(200, 64, 5)
b = rnorm(200, 65, 5)

all_value = c(a,b)
mean_all = mean(all_value)
sd_all = sd(all_value)

n = 10000
a = rnorm(n, 64, 5)
b = rnorm(n, 65, 5)
df_dist = data.frame(
  quant_value = c(
    a, b
  ),
  class = c(
    rep("A", n),
    rep("B", n)
  )
)


gg_fig3 = ggplot(df_dist, aes(x=quant_value))
gg_fig3 = gg_fig3 + ggtitle(wilcox.test(a, b)$p.value)
gg_fig3 = gg_fig3 + geom_density(fill="bleu")
gg_fig3 = gg_fig3 + stat_function(
                      fun = dnorm, n = 20000, colour = "red",
                      args = list(mean = mean_all, sd = sd_all))
gg_fig3

gg_fig3 = ggplot(df_dist, aes(x=quant_value)) + geom_histogram(binwidth=0.5,aes(y=..density.., fill=class))
gg_fig3
