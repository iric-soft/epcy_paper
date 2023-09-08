
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


dir_out = os.path.join("data", "leucegene_v2", "STAR")
train = pd.read_csv(os.path.join(dir_out, "tpm.xls"), sep="\t", index_col=0 )
dir_out = os.path.join("data", "leucegene_v2", "STAR")
test = pd.read_csv(os.path.join(dir_out, "tpm.xls"), sep="\t", index_col=0)

design = os.path.join("data", "design", "leucegene_train", "all", "design.tsv")
train_labels = pd.read_csv(design, sep="\t")
train_labels["subgroup_str"] = train_labels["subgroup"]
train_labels["subgroup"].replace(
    ["inv16", "MLL", "CK", "EVI1", "Inter", "Mono5", "NK", "Tri8", "NUP98NSD1", "t8_21", "t15_17"],
    [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
    inplace=True
)
train_labels = train_labels.sort_values(["subgroup"], ascending=True)

design = os.path.join("data", "design", "leucegene_test", "all", "design.tsv")
test_labels = pd.read_csv(design, sep="\t")
test_labels["subgroup_str"] = test_labels["subgroup"]
test_labels["subgroup"].replace(
    ["inv16", "MLL", "CK", "EVI1", "Inter", "Mono5", "NK", "Tri8", "NUP98NSD1", "t8_21", "t15_17"],
    [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
    inplace=True
)
test_labels = test_labels.sort_values(["subgroup"], ascending=True)

test = test[test_labels["sample"]]
train = train[train_labels["sample"]]

test = np.log2(test + 1)
train = np.log2(train + 1)


def run_train_label(train, train_labels, test, test_labels, features_ids_keep=None, top=None, mcc=None, nneighbors=15):
    if features_ids_keep is None:
        dir_out = os.path.join("data", "res", "leucegene_tt", "all")
    else:
        if top is not None:
            dir_out = os.path.join("data", "res", "leucegene_tt", "top"+str(top))
            train = train.loc[features_ids_keep]
            test = test.loc[features_ids_keep]
        else:
            dir_out = os.path.join("data", "res", "leucegene_tt", "mcc"+str(mcc))
            train = train.loc[features_ids_keep]
            test = test.loc[features_ids_keep]


    if not os.path.exists(dir_out):
        os.makedirs(dir_out)

    train_embedding = umap.UMAP(
        n_neighbors=nneighbors,
        min_dist=0.2,
        n_components=2,
        random_state=42,
        n_epochs = 2000
    ).fit_transform(train.T, y=train_labels["subgroup"])

    train_labels["x"] = train_embedding[:, 0]
    train_labels["y"] = train_embedding[:, 1]
    train_labels["type"] = "train"

    labels = hdbscan.HDBSCAN(
        min_samples=int(math.floor(20/4)),
        min_cluster_size=int(math.floor(20/2))
    ).fit_predict(train_embedding)
    clustered = (labels >= 0)
    p_clustered = round(np.sum(clustered) / train_labels.shape[0], 2)
    full_ARS = round(adjusted_rand_score(train_labels["subgroup"], labels), 2)
    full_MIS = round(adjusted_mutual_info_score(train_labels["subgroup"], labels, average_method='arithmetic'), 2)
    ARS = round(adjusted_rand_score(train_labels["subgroup"][clustered], labels[clustered]), 2)
    MIS = round(adjusted_mutual_info_score(train_labels["subgroup"][clustered], labels[clustered], average_method='arithmetic'), 2)

    sns.set(rc={'figure.figsize':(11.7,8.27)})
    filled_markers = ('o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X')
    color_palette=sns.color_palette("Paired", train_labels["subgroup_str"].nunique())

    sns_plot = sns.scatterplot(
        x="x",
        y="y",
        style="subgroup_str",
        hue="subgroup_str",
        palette=color_palette,
        markers=filled_markers,
        data=train_labels
    )
    box = sns_plot.get_position() # get position of figure
    sns_plot.set_position([box.x0, box.y0, box.width * 0.75, box.height]) # resize position
    sns_plot.legend(loc='center right', bbox_to_anchor=(1.45, 0.5), ncol=2)
    sns_plot.set_title("ARS=" + str(full_ARS) + ", MIS=" + str(full_MIS) + ", p=" + str(p_clustered) + ", n=" + str(train.shape[0]))

    file_out = os.path.join(dir_out, "train.pdf")
    sns_plot = sns_plot.get_figure()
    sns_plot.savefig(file_out)
    plt.close()

    mapper = umap.UMAP(
        n_neighbors=nneighbors,
        min_dist=0.2,
        n_components=2,
        random_state=42,
        n_epochs = 2000
    ).fit(train.T, np.array(train_labels["subgroup"]))

    test_embedding = mapper.transform(test.T)

    test_labels["x"] = test_embedding[:, 0]
    test_labels["y"] = test_embedding[:, 1]
    test_labels["pred"] = test_embedding[:, 0]
    test_labels["type"] = "test"

    #both_labels = pd.concat([train_labels, test_labels])

    labels = hdbscan.HDBSCAN(
        min_samples=int(math.floor(10/4)),
        min_cluster_size=int(math.floor(10/2))
    ).fit_predict(test_embedding)
    clustered = (labels >= 0)
    p_clustered = round(np.sum(clustered) / test_labels.shape[0], 2)
    full_ARS = round(adjusted_rand_score(test_labels["subgroup"], labels), 2)
    full_MIS = round(adjusted_mutual_info_score(test_labels["subgroup"], labels, average_method='arithmetic'), 2)
    ARS = round(adjusted_rand_score(test_labels["subgroup"][clustered], labels[clustered]), 2)
    MIS = round(adjusted_mutual_info_score(test_labels["subgroup"][clustered], labels[clustered], average_method='arithmetic'), 2)


    from sklearn.metrics.pairwise import euclidean_distances

    all_dist= euclidean_distances(test_embedding, train_embedding)
    print(all_dist)
    for i in range(0,len(test_labels["pred"])):
        selected_subg = train_labels["subgroup_str"].iloc[np.argpartition(all_dist[i,:], nneighbors)[:nneighbors]]
        selected_subg = selected_subg.tolist()
        #if test_labels["subgroup_str"][i] == "t8_21":
        #    print(test_labels["subgroup_str"][i] + ": " +  max(selected_subg,key=selected_subg.count))
        #    print(selected_subg)
        test_labels.loc[test_labels.iloc[i].name, "pred"] = max(selected_subg,key=selected_subg.count)

    #print(test_labels)

    dict_acc = defaultdict(list)
    sub_title = ""
    all_tp = 0
    all_tn = 0
    all_fn = 0
    all_fp = 0
    all_sum = 0
    cpt = 0
    for subg in set(test_labels["subgroup_str"]):
        tp = sum(test_labels["pred"][test_labels["subgroup_str"] == subg] == subg)
        fp = sum(test_labels["pred"][test_labels["subgroup_str"] != subg] == subg)
        tn = sum(test_labels["pred"][test_labels["subgroup_str"] != subg] != subg)
        fn = sum(test_labels["pred"][test_labels["subgroup_str"] == subg] != subg)

        acc = (tp + tn) / len(test_labels["pred"])

        d1 = tp + fp
        d2 = tp + fn
        d3 = tn + fp
        d4 = tn + fn
        n1 = tp * tn
        n2 = fp * fn
        mcc = ((n1 - n2) / math.sqrt(d1 * d2 * d3 * d4) if d1!=0 and d2!=0 and d3!=0 and d4!=0 else 0)

        all_tp += tp
        all_tn += tn
        all_fn += fn
        all_fp += fp
        all_sum += len(test_labels["pred"])

        dict_acc['tp'].append(tp)
        dict_acc['tn'].append(tn)
        dict_acc['acc'].append(acc)
        dict_acc['mcc'].append(mcc)
        dict_acc['subg'].append(subg)
        if cpt % 4 == 0:
            sub_title += "\n"
        sub_title += subg + ": " + str(round(mcc, 2)) + ", "

        cpt += 1

    all_acc = round((all_tp + all_tn) / all_sum, 2)
    d1 = all_tp + all_fp
    d2 = all_tp + all_fn
    d3 = all_tn + all_fp
    d4 = all_tn + all_fn
    n1 = all_tp * all_tn
    n2 = all_fp * all_fn
    all_mcc = ((n1 - n2) / math.sqrt(d1 * d2 * d3 * d4) if d1!=0 and d2!=0 and d3!=0 and d4!=0 else 0)
    all_mcc = round(all_mcc, 2)

    sns_plot = sns.scatterplot(
        x="x",
        y="y",
        style="subgroup_str",
        hue="subgroup_str",
        palette=color_palette,
        markers=filled_markers,
        data=test_labels
    )
    box = sns_plot.get_position() # get position of figure
    sns_plot.set_position([box.x0, box.y0, box.width * 0.75, box.height]) # resize position
    sns_plot.legend(loc='center right', bbox_to_anchor=(1.45, 0.5), ncol=2)
    sns_plot.set_title(
        "ACC=" + str(all_acc) +
        ", MCC=" + str(all_mcc) +
        ", ARS=" + str(full_ARS) +
        ", MIS=" + str(full_MIS) +
        ", p=" + str(p_clustered) +
        ", n=" + str(test.shape[0]) +
        "\n" + sub_title
    )

    file_out = os.path.join(dir_out, "test.pdf")
    sns_plot = sns_plot.get_figure()
    sns_plot.savefig(file_out)
    plt.close()





def read_diff_table(file_name, path_dir):
    path_file = path_dir
    path_file = os.path.join(path_file, "STAR", "tpm", file_name)
    df_diff = pd.read_csv(path_file, sep="\t")
    df_diff["abs_l2fc"] = df_diff.l2fc.abs()
    df_diff = df_diff.sort_values(["kernel_mcc", "abs_l2fc"], ascending=False)

    return(df_diff)



run_train_label(train, train_labels, test, test_labels)

mcc = None
top = None
for top in [10,20,30,40,50,60,70]:
#for mcc in [0.9,0.8,0.7,0.6,0.5,0.4,0.3]:
    features_ids_keep = []
    path_design="/u/eaudemard/project/epcy_paper/data/design/leucegene_train"
    for design in ["32_MLL"]: #["21_inv16", "32_MLL", "61_CK", "9_EVI1", "58_Inter", "15_Mono5", "183_NK", "20_t15_17", "12_t8_21", "12_Tri8", "4_NUP98NSD1"]:
        print(design)
        dir_design = os.path.join(path_design, design)
        df_diff = read_diff_table("prediction_capability.xls", dir_design)
        if mcc is not None:
            df_diff_selected = df_diff.loc[df_diff["kernel_mcc"] >= mcc]
            features_ids_keep += df_diff_selected["id"].tolist()
        else:
            features_ids_keep += df_diff["id"][:top].tolist()

    run_train_label(train, train_labels, test, test_labels, features_ids_keep, top, mcc)
















features_ids_keep = []
path_design="/u/eaudemard/project/epcy_paper/data/design/leucegene_train"
for design in ["21_inv16", "32_MLL", "61_CK", "9_EVI1", "58_Inter", "15_Mono5", "183_NK", "20_t15_17", "12_t8_21", "12_Tri8", "4_NUP98NSD1"]:
    print(design)
    dir_design = os.path.join(path_design, design)
    df_diff = read_diff_table("prediction_capability.xls", dir_design)
    features_ids_keep += df_diff["id"][:10].tolist()

dir_out = os.path.join("data", "TCGA_LAML", "htseq")
tcga_laml = pd.read_csv(os.path.join(dir_out, "tpm.xls"), sep="\t", index_col=0 )
tcga_laml_top = tcga_laml.loc[features_ids_keep]
tcga_laml_top = tcga_laml_top[np.isfinite(tcga_laml_top['TCGA_AB_3001'])]

design = os.path.join("data", "design", "TCGA_LAML", "all", "design.tsv")
tcga_laml_labels = pd.read_csv(design, sep="\t")
tcga_laml_labels["subgroup_str"] = tcga_laml_labels["subgroup"]
#tcga_laml_labels["subgroup"].replace(
#    ["inv16", "MLL", "CK", "EVI1", "Inter", "Mono5", "NK", "Tri8", "NUP98NSD1", "t8_21", "t15_17"],
#    [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
#    inplace=True
#)
#tcga_laml_labels = tcga_laml_labels.sort_values(["subgroup"], ascending=True)


tcga_laml_top = tcga_laml_top[tcga_laml_labels["sample"]]
tcga_laml_top = np.log2(tcga_laml_top + 1)

tcga_laml_embedding = umap.UMAP(
    n_neighbors=15,
    min_dist=0.2,
    n_components=2,
    random_state=42,
    n_epochs = 2000
).fit_transform(tcga_laml_top.T)

tcga_laml_labels["x"] = tcga_laml_embedding[:, 0]
tcga_laml_labels["y"] = tcga_laml_embedding[:, 1]
tcga_laml_labels["type"] = "tcga_laml"

sns.set(rc={'figure.figsize':(11.7,8.27)})
filled_markers = ('o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X', 'o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X')

sns_plot = sns.scatterplot(
    x="x",
    y="y",
    hue="subgroup_str",
    palette=color_palette,
    data=tcga_laml_labels,
    style="subgroup_str",
    markers=filled_markers
)
box = sns_plot.get_position() # get position of figure
sns_plot.set_position([box.x0, box.y0, box.width * 0.75, box.height]) # resize position
sns_plot.legend(loc='center right', bbox_to_anchor=(1.45, 0.5), ncol=2)


plt.show()

file_out = os.path.join(dir_out, "train.pdf")
sns_plot = sns_plot.get_figure()
sns_plot.savefig(file_out)
plt.close()
