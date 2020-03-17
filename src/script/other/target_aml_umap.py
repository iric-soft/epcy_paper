
import csv
import os
import umap
import hdbscan
import math

import numpy as np
import pandas as pd

import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm

from mnist import MNIST
from collections import defaultdict

from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score



design = os.path.join("data", "design", "TARGET_AML", "all", "design.tsv")
df_design = pd.read_csv(design, sep="\t")
df_design["subgroup_str"] = df_design["subgroup"]
df_design["subgroup"].replace(
    ["inv16", "MLL", "Other", "Normal", "t8_21", "Unknown"],
    [0, 1, 2, 3, 4, 5],
    inplace=True
)
df_design = df_design.sort_values(["subgroup"], ascending=True)

dir_exp = os.path.join("data", "TARGET_AML", "htseq")
df_exp = pd.read_csv(os.path.join(dir_exp, "tpm.xls"), sep="\t", index_col=0 )
df_exp = df_exp[df_design["sample"]]
df_exp = np.log2(df_exp + 1)


target_aml_embedding = umap.UMAP(
    n_neighbors=15,
    min_dist=0.2,
    n_components=2,
    random_state=42,
    n_epochs = 2000
).fit_transform(df_exp.T)

df_design["x"] = target_aml_embedding[:, 0]
df_design["y"] = target_aml_embedding[:, 1]


dir_out = os.path.join("data", "res", "TARGET_AML")
if not os.path.exists(dir_out):
    os.makedirs(dir_out)


sns.set(rc={'figure.figsize':(11.7,8.27)})
filled_markers = ('o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X', 'o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X')
color_palette=sns.color_palette("Paired", df_design["subgroup_str"].nunique())

sns_plot = sns.scatterplot(
    x="x",
    y="y",
    hue="subgroup_str",
    palette=color_palette,
    data=df_design,
    style="subgroup_str",
    markers=filled_markers
)
box = sns_plot.get_position() # get position of figure
sns_plot.set_position([box.x0, box.y0, box.width * 0.75, box.height]) # resize position
sns_plot.legend(loc='center right', bbox_to_anchor=(1.45, 0.5), ncol=1)


file_out = os.path.join(dir_out, "umap_all_genes.pdf")
sns_plot = sns_plot.get_figure()
sns_plot.savefig(file_out)
plt.close()



def read_diff_table(file_name, path_dir):
    path_file = path_dir
    path_file = os.path.join(path_file, "STAR", "tpm", file_name)
    df_diff = pd.read_csv(path_file, sep="\t")
    df_diff.rename(str.upper, axis='columns', inplace=True)
    df_diff["ABS_L2FC"] = df_diff.L2FC.abs()
    df_diff = df_diff.sort_values(["KERNEL_MCC", "ABS_L2FC"], ascending=False)

    return(df_diff)

top = 20
features_ids_keep = []
path_design="/u/eaudemard/project/epcy_paper/data/design/leucegene_v2"
for design in ["31_inv16", "48_MLL", "275_NK", "18_t8_21"]:
    print(design)
    dir_design = os.path.join(path_design, design)
    df_diff = read_diff_table("prediction_capability.xls", dir_design)
    features_ids_keep += df_diff["ID"][:top].tolist()

dir_exp = os.path.join("data", "TARGET_AML", "htseq")
df_exp = pd.read_csv(os.path.join(dir_exp, "tpm.xls"), sep="\t", index_col=0 )
df_exp = df_exp.reindex(features_ids_keep, method=None)
df_exp.dropna(inplace=True)
df_exp = df_exp[df_design["sample"]]
df_exp = np.log2(df_exp + 1)


target_aml_embedding = umap.UMAP(
    n_neighbors=15,
    min_dist=0.2,
    n_components=2,
    random_state=42,
    n_epochs = 2000
).fit_transform(df_exp.T)

df_design["x"] = target_aml_embedding[:, 0]
df_design["y"] = target_aml_embedding[:, 1]

sns.set(rc={'figure.figsize':(11.7,8.27)})
filled_markers = ('o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X', 'o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X')
color_palette=sns.color_palette("Paired", df_design["subgroup_str"].nunique())

sns_plot = sns.scatterplot(
    x="x",
    y="y",
    hue="subgroup_str",
    palette=color_palette,
    data=df_design,
    style="subgroup_str",
    markers=filled_markers
)
box = sns_plot.get_position() # get position of figure
sns_plot.set_position([box.x0, box.y0, box.width * 0.75, box.height]) # resize position
sns_plot.legend(loc='center right', bbox_to_anchor=(1.45, 0.5), ncol=1)


file_out = os.path.join(dir_out, "umap_leucegene_genes.pdf")
sns_plot = sns_plot.get_figure()
sns_plot.savefig(file_out)
plt.close()
