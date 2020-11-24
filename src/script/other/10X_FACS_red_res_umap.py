import os
import h5py as h5
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import umap

from scipy import stats

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'DejaVu Sans'
mpl.rcParams['pdf.fonttype'] = 42

col_pal = plt.get_cmap("tab10")
biomart = pd.read_csv('./data/other/ensembl_gene.txt', sep="\t")
biomart.shape
biomart = biomart[biomart["Gene.names"].str.match('^(RPL)') == False]
biomart = biomart[biomart["Gene.names"].str.match('^(RPS)') == False]

design = pd.read_csv('./data/design/10X_FACS_reduce/all/design.tsv', sep="\t")
matrix_num_samples = design.shape[0]

counts = pd.read_csv('./data/10X_FACS_reduce/cellranger/readcounts.xls', index_col=0, sep="\t" )
counts = counts.T
total_mass = np.nansum(counts, axis=1)
med_count = np.median(total_mass)
counts = counts.div(total_mass, axis=0) * med_count
counts = np.log2(counts + 1)

path_res = os.path.join("./data/res/10X_FACS_reduce_umap/")
if not os.path.exists(path_res):
    os.makedirs(path_res)

standard_embedding = umap.UMAP(random_state=42).fit_transform(counts)

counts["x"] = standard_embedding[:, 0]
counts["y"] = standard_embedding[:, 1]
counts = counts.reindex(design['sample'].values)
counts["type"] = design.subgroup.values

fig, ax = plt.subplots()
for type, color in zip(design.subgroup.unique(), col_pal.colors):
    df_tmp = counts[counts['type'] == type]
    sns_plot = ax.scatter(
        df_tmp.x, df_tmp.y, label=type,
        alpha=.8, s=1, c=np.array([color])
    )

ax.legend()
fig_file = os.path.join(path_res, "umap_test.pdf")
sns_plot.figure.savefig(fig_file)
plt.close()

design = "948_cd34"
path_res_design = os.path.join(path_res, design)
if not os.path.exists(path_res_design):
    os.makedirs(path_res_design)

path_design = os.path.join("./data/design/10X_FACS_reduce/", design, "cellranger", "readcounts")
epcy = pd.read_csv(os.path.join(path_design, "predictive_capability.xls"), index_col=0, sep="\t")
mast = pd.read_csv(os.path.join(path_design, "mast_genes.xls"), index_col=0, sep="\t")
mast["pval"] = -np.log10(mast["pval"].values)
mast = mast.replace([np.inf, -np.inf], 340)
trend = pd.read_csv(os.path.join(path_design, "limma_trend_genes.xls"), index_col=0, sep="\t")
trend["adj.P.Val"] = -np.log10(trend["adj.P.Val"].values)
trend = trend.replace([np.inf, -np.inf], 340)
all = epcy.merge(mast, left_index=True, right_index=True, how='outer')
all = all.merge(trend, left_index=True, right_index=True, how='outer')
all = all.merge(biomart, left_index=True, right_on='ens_gene', how='inner')
all.kernel_mcc = round(all.kernel_mcc, 4)
all.kernel_ppv = round(all.kernel_ppv, 4)
all.kernel_npv = round(all.kernel_npv, 4)
all.kernel_npv = round(all.kernel_npv, 4)
all.pval = round(all.pval, 4)
all["adj.P.Val"] = round(all["adj.P.Val"], 4)


top = 50
all.sort_values(by='kernel_ppv', ascending=False, inplace=True)
d = {
    'type': ["ppv"] * top,
    'method': ["max"] * top,
    'value': all.kernel_ppv[0:top].values,
    'top': [top] * top
}
df_res = pd.DataFrame(data=d)

all.sort_values(by='kernel_npv', ascending=False, inplace=True)
d = {
    'type': ["npv"] * top,
    'method': ["max"] * top,
    'value': all.kernel_npv[0:top].values,
    'top': [top] * top
}
df_tmp = pd.DataFrame(data=d)
df_res = df_res.append(df_tmp)

all.sort_values(by='kernel_mcc', ascending=False, inplace=True)
d = {
    'type': ["ppv"] * top,
    'method': ["epcy"] * top,
    'value': all.kernel_ppv[0:top].values,
    'top': [top] * top
}
df_tmp = pd.DataFrame(data=d)
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

all.sort_values(by='pval', ascending=False, inplace=True)
d = {
    'type': ["ppv"] * top,
    'method': ["mast"] * top,
    'value': all.kernel_ppv[0:top].values,
    'top': [top] * top
}
df_tmp = pd.DataFrame(data=d)
df_res = df_res.append(df_tmp)
d = {
    'type': ["npv"] * top,
    'method': ["mast"] * top,
    'value': all.kernel_npv[0:top].values,
    'top': [top] * top
}
df_tmp = pd.DataFrame(data=d)
df_res = df_res.append(df_tmp)

all.sort_values(by='adj.P.Val', ascending=False, inplace=True)
d = {
    'type': ["ppv"] * top,
    'method': ["trend"] * top,
    'value': all.kernel_ppv[0:top].values,
    'top': [top] * top
}
df_tmp = pd.DataFrame(data=d)
df_res = df_res.append(df_tmp)
d = {
    'type': ["npv"] * top,
    'method': ["trend"] * top,
    'value': all.kernel_npv[0:top].values,
    'top': [top] * top
}
df_tmp = pd.DataFrame(data=d)
df_res = df_res.append(df_tmp)

col_pal_method = sns.color_palette("colorblind")
col_pal_method = [col_pal_method[i] for i in [4, 1, 2, 9]]
col_pal_method = [col_pal_method[3], col_pal_method[0], col_pal_method[2]]

df_res_npv = df_res.loc[df_res.type == "npv"]
npv_epcy_trend = stats.mannwhitneyu(
    df_res_npv.value.loc[df_res_npv.method == "epcy"],
    df_res_npv.value.loc[df_res_npv.method == "trend"]
)[1]
npv_epcy_mast = stats.mannwhitneyu(
    df_res_npv.value.loc[df_res_npv.method == "epcy"],
    df_res_npv.value.loc[df_res_npv.method == "mast"]
)[1]
npv_mast_trend = stats.mannwhitneyu(
    df_res_npv.value.loc[df_res_npv.method == "mast"],
    df_res_npv.value.loc[df_res_npv.method == "trend"]
)[1]

df_res_ppv = df_res.loc[df_res.type == "ppv"]
ppv_epcy_trend = stats.mannwhitneyu(
    df_res_ppv.value.loc[df_res_ppv.method == "epcy"],
    df_res_ppv.value.loc[df_res_ppv.method == "trend"]
)[1]
ppv_epcy_mast = stats.mannwhitneyu(
    df_res_ppv.value.loc[df_res_ppv.method == "epcy"],
    df_res_ppv.value.loc[df_res_ppv.method == "mast"]
)[1]
ppv_mast_trend = stats.mannwhitneyu(
    df_res_ppv.value.loc[df_res_ppv.method == "mast"],
    df_res_ppv.value.loc[df_res_ppv.method == "trend"]
)[1]

fig = plt.figure(figsize=(10, 7))
gs = plt.GridSpec(1, 2)

ax_ppv = fig.add_subplot(gs[0, 0])
ax_npv = fig.add_subplot(gs[0, 1])

sns_plot = sns.boxplot(x="method", y="value",
                       hue="method", ax=ax_ppv, linewidth=1,
                       palette=col_pal_method,
                       data=df_res.loc[df_res.type == "ppv"])
sns_plot.set_title("PPV distribution of top " + str(top) + "\nepcy_trend=" + str(ppv_epcy_trend) + "\nepcy_mast=" + str(ppv_epcy_mast) + "\ntrend_mast=" + str(ppv_mast_trend))
sns_plot.set(ylim=(0, 1))

sns_plot = sns.boxplot(x="method", y="value",
                       hue="method", ax=ax_npv, linewidth=1,
                       palette=col_pal_method,
                       data=df_res.loc[df_res.type == "npv"])
sns_plot.set_title("NPV distribution of top " + str(top) + "\nepcy_trend=" + str(npv_epcy_trend) + "\nepcy_mast=" + str(npv_epcy_mast) + "\ntrend_mast=" + str(npv_mast_trend))
sns_plot.set(ylim=(0, 1))

fig_file = os.path.join(path_res, design, "ppv_npv.pdf")
sns_plot.figure.savefig(fig_file)
plt.close('all')


def gene_umap(df, path_res, counts):
    if not os.path.exists(path_res):
        os.makedirs(path_res)
    for i in range(df.shape[0]):
        counts.sort_values(by=[df["ens_gene"].values[i]], inplace=True)
        fig_file = os.path.join(path_res, "umap_" + df["ens_gene"].values[i] + "_" + df["Gene.names"].values[i] + ".pdf")
        sns_plot = plt.scatter(
            counts.x, counts.y,
            alpha=.8, c=counts[df["ens_gene"].values[i]], s=1,
            cmap='jet'
        )
        plt.title(df["ens_gene"].values[i] + " " + df["Gene.names"].values[i] + "\n MCC=" + str(df["kernel_mcc"].values[i]) + ", mast=" + str(df["pval"].values[i]) + ", trend=" + str(df["adj.P.Val"].values[i]) + "\n PPV=" + str(df["kernel_ppv"].values[i])+ ", NPV=" + str(df["kernel_npv"].values[i]))
        cbar = plt.colorbar()
        sns_plot.figure.savefig(fig_file)
        plt.close()


##############
tmp = all.loc[(all["kernel_mcc"] < 0.15) & (all["pval"] > 300)]
print(tmp)
path_res_gene = os.path.join(path_res, design, "mast_p_300_mcc_m_15")
gene_umap(tmp, path_res_gene, counts)

##############
tmp = all.loc[(all["kernel_mcc"] > 0.3) & (all["kernel_mcc"] < 0.4) & (all["adj.P.Val"] < 100)]
print(tmp)
path_res_gene = os.path.join(path_res, design, "trend_m_100_mcc_30_40")
gene_umap(tmp, path_res_gene, counts)

##############
tmp = all.loc[(all["kernel_mcc"] > 0.4) & (all["adj.P.Val"] < 150)]
print(tmp)
path_res_gene = os.path.join(path_res, design, "trend_m_150_mcc_p_40")
gene_umap(tmp, path_res_gene, counts)


##############
tmp = all.loc[(all["kernel_mcc"] < 0.15) & (all["adj.P.Val"] > 300)]
print(tmp)
path_res_gene = os.path.join(path_res, design, "trend_p_300_mcc_m_15")
gene_umap(tmp, path_res_gene, counts)

##############
tmp = all.loc[(all["kernel_mcc"] < 0.2) & (all["adj.P.Val"] == 340)]
print(tmp)
path_res_gene = os.path.join(path_res, design, "trend_340_mcc_m_20")
gene_umap(tmp, path_res_gene, counts)


##############
tmp = all.loc[(all["kernel_mcc"] > 0.6)]
print(tmp)
path_res_gene = os.path.join(path_res, design, "mcc_p_60")
gene_umap(tmp, path_res_gene, counts)

##############
tmp = all.loc[(all["adj.P.Val"] == 340) | (all["pval"] == 340)]
print(tmp)
path_res_gene = os.path.join(path_res, design, "340")
gene_umap(tmp, path_res_gene, counts)


##############
tmp = all.loc[(all["kernel_mcc"] > 0.25) & (all["kernel_mcc"] < 0.35) & (all["adj.P.Val"] == 340) & (all["pval"] == 340)]
print(tmp)
path_res_gene = os.path.join(path_res, design, "340_mcc_25_35")
gene_umap(tmp, path_res_gene, counts)

##############
all.sort_values(by='kernel_ppv', ascending=False, inplace=True)
tmp = all.loc[(all["kernel_ppv"] >= all.kernel_ppv.values[9])]
print(tmp)
path_res_gene = os.path.join(path_res, design, "top10_ppv")
gene_umap(tmp, path_res_gene, counts)

##############
all.sort_values(by='kernel_npv', ascending=False, inplace=True)
tmp = all.loc[(all["kernel_npv"] >= all.kernel_npv.values[9])]
print(tmp)
path_res_gene = os.path.join(path_res, design, "top10_npv")
gene_umap(tmp, path_res_gene, counts)
