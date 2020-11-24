import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import umap

from scipy import stats
from collections import defaultdict

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'DejaVu Sans'
mpl.rcParams['pdf.fonttype'] = 42

col_pal = plt.get_cmap("tab10")

df_design = pd.read_csv('./data/design/leucegene3/30_inv16/design.tsv', sep="\t")
df_design_query = df_design[df_design.subgroup == "Query"]
df_design_ref = df_design[df_design.subgroup != "Query"]

path_design = os.path.join("./data/design/leucegene3/", "30_inv16", "STAR_RSEM", "readcounts")
epcy = pd.read_csv(os.path.join(path_design, "predictive_capability.xls"), index_col=0, sep="\t")
deseq = pd.read_csv(os.path.join(path_design, "deseq2_genes.xls"), index_col=0, sep="\t")
deseq["padj"] = -np.log10(deseq["padj"].values)
voom = pd.read_csv(os.path.join(path_design, "limma_voom_genes.xls"), index_col=0, sep="\t")
voom["adj.P.Val"] = -np.log10(voom["adj.P.Val"].values)
edger = pd.read_csv(os.path.join(path_design, "edger_genes.xls"), index_col=0, sep="\t")
edger["FDR"] = -np.log10(edger["FDR"].values)

all = epcy.merge(deseq, left_index=True, right_index=True, how='outer')
all = all.merge(voom, left_index=True, right_index=True, how='outer')
all = all.merge(edger, left_index=True, right_index=True, how='outer')

all.kernel_mcc = round(all.kernel_mcc, 4)
all.kernel_ppv = round(all.kernel_ppv, 4)
all.kernel_npv = round(all.kernel_npv, 4)
all.PValue = round(all.PValue, 4)
all["adj.P.Val"] = round(all["adj.P.Val"], 4)



path_res = os.path.join("./data/res/leucegene3/")
if not os.path.exists(path_res):
    os.makedirs(path_res)

design = "30_inv16"
path_res_design = os.path.join(path_res, design)
if not os.path.exists(path_res_design):
    os.makedirs(path_res_design)

top = 20
all.sort_values(by='kernel_mcc', ascending=False, inplace=True)
d = {
    'method': ["EPCY"] * top,
    'log2FC': all.l2fc[0:top].values,
    'top': [top] * top
}
df_res = pd.DataFrame(data=d)

all.sort_values(by='padj', ascending=False, inplace=True)
d = {
    'method': ["DESeq2"] * top,
    'log2FC': all.l2fc[0:top].values,
    'top': [top] * top
}
df_tmp = pd.DataFrame(data=d)
df_res = df_res.append(df_tmp)

all.sort_values(by='adj.P.Val', ascending=False, inplace=True)
d = {
    'method': ["Voom"] * top,
    'log2FC': all.l2fc[0:top].values,
    'top': [top] * top
}
df_tmp = pd.DataFrame(data=d)
df_res = df_res.append(df_tmp)

col_pal_method = sns.color_palette("colorblind")
col_pal_method = [col_pal_method[i] for i in [4, 1, 2, 9]]
col_pal_method = [col_pal_method[3], col_pal_method[0], col_pal_method[2]]

epcy_deseq = stats.mannwhitneyu(
    df_res.log2FC.loc[df_res.method == "EPCY"],
    df_res.log2FC.loc[df_res.method == "DESeq2"]
)[1]
epcy_voom = stats.mannwhitneyu(
    df_res.log2FC.loc[df_res.method == "EPCY"],
    df_res.log2FC.loc[df_res.method == "Voom"]
)[1]
deseq_voom = stats.mannwhitneyu(
    df_res.log2FC.loc[df_res.method == "DESeq2"],
    df_res.log2FC.loc[df_res.method == "Voom"]
)[1]

sns_plot = sns.boxplot(x="method", y="log2FC",
                hue="method", linewidth=1,
                palette=col_pal_method,
                data=df_res)
sns_plot.set_title("log2FC distribution of top " + str(top) + "\nepcy_deseq=" + str(epcy_deseq) + "\nepcy_voom=" + str(epcy_voom) + "\ndeseq_voom=" + str(deseq_voom))
#sns_plot.set(ylim=(-0.05, 1.2))

fig_out = os.path.join(path_res_design, "log2FC_top" + str(top) + ".pdf")
sns_plot.figure.savefig(fig_out)
plt.close('all')


tmp = all.loc[all["kernel_mcc"] > 0.9]
print(tmp)

tmp = all.loc[(all["kernel_mcc"] > 0.6) & (all["adj.P.Val"] > 100) & (all["padj"] < 30)]
print(tmp)

tmp = all.loc[(all["kernel_mcc"] > 0.2) & (all["adj.P.Val"] < 15) & (all["padj"] > 75)]
print(tmp)

tmp = all.loc[(all["kernel_mcc"] < 0.2) & (all["adj.P.Val"] < 15) & (all["padj"] > 120)]
print(tmp)
