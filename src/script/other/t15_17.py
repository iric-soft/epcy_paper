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

df_design = pd.read_csv('./data/design/leucegene3/30_t15_17/design.tsv', sep="\t")
df_design_query = df_design[df_design.subgroup == "Query"]
df_design_ref = df_design[df_design.subgroup != "Query"]

path_design = os.path.join("./data/design/leucegene3/", "30_t15_17", "STAR_RSEM", "readcounts")
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

design = "30_t15_17"
path_res_design = os.path.join(path_res, design)
if not os.path.exists(path_res_design):
    os.makedirs(path_res_design)


d = defaultdict(list)
all.sort_values(by='kernel_mcc', ascending=False, inplace=True)
for top in range(1, 51, 1):
    num_pos = np.sum(all.l2fc[0:top] >= 0)
    num_neg = top - num_pos
    d['top'].append(top)
    d['method'].append('EPCY')
    d['type'].append('highly')
    d['value'].append(num_pos)
    d['top'].append(top)
    d['method'].append('EPCY')
    d['type'].append('lowly')
    d['value'].append(num_neg)

all.sort_values(by='padj', ascending=False, inplace=True)
for top in range(1, 51, 1):
    num_pos = np.sum(all.l2fc[0:top] >= 0)
    num_neg = top - num_pos
    d['top'].append(top)
    d['method'].append('DESeq2')
    d['type'].append('highly')
    d['value'].append(num_pos)
    d['top'].append(top)
    d['method'].append('DESeq2')
    d['type'].append('lowly')
    d['value'].append(num_neg)

all.sort_values(by='adj.P.Val', ascending=False, inplace=True)
for top in range(1, 51, 1):
    num_pos = np.sum(all.l2fc[0:top] >= 0)
    num_neg = top - num_pos
    d['top'].append(top)
    d['method'].append('voom')
    d['type'].append('highly')
    d['value'].append(num_pos)
    d['top'].append(top)
    d['method'].append('voom')
    d['type'].append('lowly')
    d['value'].append(num_neg)

df_pos_neg = pd.DataFrame(d)

sns_plot = sns.relplot(
    data=df_pos_neg, x="top", y="value",
    col="type", hue="method",
    kind="line"
)
fig_out = os.path.join(path_res_design, "num_pos_neg.pdf")
sns_plot.savefig(fig_out)
plt.close()

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






top = 100
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










all.sort_values(by='kernel_mcc', ascending=False, inplace=True)
d = {
    'type': ["ppv"] * top,
    'method': ["epcy"] * top,
    'value': all.kernel_ppv[0:top].values,
    'top': [top] * top
}
#df_tmp = pd.DataFrame(data=d)
#df_res = df_res.append(df_tmp)
df_res = pd.DataFrame(data=d)
d = {
    'type': ["npv"] * top,
    'method': ["epcy"] * top,
    'value': all.kernel_npv[0:top].values,
    'top': [top] * top
}
df_tmp = pd.DataFrame(data=d)
df_res = df_res.append(df_tmp)

all.sort_values(by='padj', ascending=False, inplace=True)
d = {
    'type': ["ppv"] * top,
    'method': ["deseq2"] * top,
    'value': all.kernel_ppv[0:top].values,
    'top': [top] * top
}
df_tmp = pd.DataFrame(data=d)
df_res = df_res.append(df_tmp)
d = {
    'type': ["npv"] * top,
    'method': ["deseq2"] * top,
    'value': all.kernel_npv[0:top].values,
    'top': [top] * top
}
df_tmp = pd.DataFrame(data=d)
df_res = df_res.append(df_tmp)

all.sort_values(by='adj.P.Val', ascending=False, inplace=True)
d = {
    'type': ["ppv"] * top,
    'method': ["voom"] * top,
    'value': all.kernel_ppv[0:top].values,
    'top': [top] * top
}
df_tmp = pd.DataFrame(data=d)
df_res = df_res.append(df_tmp)
d = {
    'type': ["npv"] * top,
    'method': ["voom"] * top,
    'value': all.kernel_npv[0:top].values,
    'top': [top] * top
}
df_tmp = pd.DataFrame(data=d)
df_res = df_res.append(df_tmp)

all.sort_values(by='FDR', ascending=False, inplace=True)
d = {
    'type': ["ppv"] * top,
    'method': ["edger"] * top,
    'value': all.kernel_ppv[0:top].values,
    'top': [top] * top
}
df_tmp = pd.DataFrame(data=d)
df_res = df_res.append(df_tmp)
d = {
    'type': ["npv"] * top,
    'method': ["edger"] * top,
    'value': all.kernel_npv[0:top].values,
    'top': [top] * top
}
df_tmp = pd.DataFrame(data=d)
df_res = df_res.append(df_tmp)

df_res2 = df_res.loc[df_res.method!="edger"]

col_pal_method = sns.color_palette("colorblind")
col_pal_method = [col_pal_method[i] for i in [4, 1, 2, 9]]
col_pal_method = [col_pal_method[3], col_pal_method[0], col_pal_method[2]]

df_res_npv = df_res2.loc[df_res2.type == "npv"]
npv_epcy_deseq = stats.mannwhitneyu(
    df_res_npv.value.loc[df_res_npv.method == "epcy"],
    df_res_npv.value.loc[df_res_npv.method == "deseq2"]
)[1]
npv_epcy_voom = stats.mannwhitneyu(
    df_res_npv.value.loc[df_res_npv.method == "epcy"],
    df_res_npv.value.loc[df_res_npv.method == "voom"]
)[1]
npv_deseq_voom = stats.mannwhitneyu(
    df_res_npv.value.loc[df_res_npv.method == "deseq2"],
    df_res_npv.value.loc[df_res_npv.method == "voom"]
)[1]

df_res_ppv = df_res2.loc[df_res2.type == "ppv"]
ppv_epcy_deseq = stats.mannwhitneyu(
    df_res_ppv.value.loc[df_res_ppv.method == "epcy"],
    df_res_ppv.value.loc[df_res_ppv.method == "deseq2"]
)[1]
ppv_epcy_voom = stats.mannwhitneyu(
    df_res_ppv.value.loc[df_res_ppv.method == "epcy"],
    df_res_ppv.value.loc[df_res_ppv.method == "voom"]
)[1]
ppv_deseq_voom = stats.mannwhitneyu(
    df_res_ppv.value.loc[df_res_ppv.method == "deseq2"],
    df_res_ppv.value.loc[df_res_ppv.method == "voom"]
)[1]

fig = plt.figure(figsize=(10, 7))
gs = plt.GridSpec(1, 2)

ax_ppv = fig.add_subplot(gs[0, 0])
ax_npv = fig.add_subplot(gs[0, 1])

sns_plot = sns.boxplot(x="method", y="value",
                hue="method", ax=ax_ppv, linewidth=1,
                palette=col_pal_method,
                data=df_res2.loc[df_res2.type=="ppv"])
sns_plot.set_title("PPV distribution of top " + str(top) + "\nepcy_deseq=" + str(ppv_epcy_deseq) + "\nepcy_voom=" + str(ppv_epcy_voom) + "\ndeseq_voom=" + str(ppv_deseq_voom))
sns_plot.set(ylim=(-0.05, 1.2))

sns_plot = sns.boxplot(x="method", y="value",
                hue="method", ax=ax_npv, linewidth=1,
                palette=col_pal_method,
                data=df_res2.loc[df_res2.type=="npv"])
sns_plot.set_title("NPV distribution of top " + str(top) + "\nepcy_deseq=" + str(npv_epcy_deseq) + "\nepcy_voom=" + str(npv_epcy_voom) + "\ndeseq_voom=" + str(npv_deseq_voom))
sns_plot.set(ylim=(0.843, 1.030))

fig_file = os.path.join(path_res, design, "ppv_npv_" + str(top)+ ".pdf")
sns_plot.figure.savefig(fig_file)
plt.close('all')







tmp = all.loc[(all["kernel_mcc"] < 0.3) & (all["padj"] > 60) & (all["padj"] < 80)]
print(tmp)
