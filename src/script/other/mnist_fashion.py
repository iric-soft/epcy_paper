
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

mndata = MNIST('fashion-mnist/data/fashion')
train, train_labels = mndata.load_training()
test, test_labels = mndata.load_testing()

test = np.array([np.array(xi) for xi in test])
train = np.array([np.array(xi) for xi in train])

test = test.T
train = train.T

head_test=list()
head_test.extend(["s" + str(x) for x in range(0, test.shape[1])])

head_train=list()
head_train.extend(["s" + str(x) for x in range(0, train.shape[1])])

test = pd.DataFrame(test,
    index=["f" + str(x) for x in range(0, test.shape[0])],
    columns=head_test
)

train = pd.DataFrame(train,
    index=["f" + str(x) for x in range(0, train.shape[0])],
    columns=head_train
)

train.reset_index(inplace=True)
train_header = train.columns.values
train_header[0] = "ID"
train.columns = train_header

test.reset_index(inplace=True)
test_header = test.columns.values
test_header[0] = "ID"
test.columns = test_header

dir_out = os.path.join("data", "mnist_fashion", "train")
if not os.path.exists(dir_out):
    os.makedirs(dir_out)

train.to_csv(os.path.join(dir_out, "tpm.xls"), header=True, index=False, sep="\t")

dir_out = os.path.join("data", "mnist_fashion", "test")
if not os.path.exists(dir_out):
    os.makedirs(dir_out)

test.to_csv(os.path.join(dir_out, "tpm.xls"), header=True, index=False, sep="\t")






test_labels = pd.DataFrame(test_labels,
    columns=["subgroup"]
)
test_labels["subgroup_i"] = test_labels["subgroup"]

train_labels = pd.DataFrame(train_labels,
    columns=["subgroup"]
)
train_labels["subgroup_i"] = train_labels["subgroup"]

test_labels["subgroup"].replace(
    [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
    ["T-shirt", "Trouser", "Pullover", "Dress", "Coat", "Sandal", "Shirt", "Sneaker", "Bag", "Ankle_boot"],
    inplace=True
)

train_labels["subgroup"].replace(
    [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
    ["T-shirt", "Trouser", "Pullover", "Dress", "Coat", "Sandal", "Shirt", "Sneaker", "Bag", "Ankle_boot"],
    inplace=True
)

test_labels.reset_index(inplace=True)
test_header = test_labels.columns.values
test_header[0] = "sample"
test_labels.columns = test_header
test_labels["sample"] = ["s" + str(x) for x in test_labels["sample"].values]

train_labels.reset_index(inplace=True)
train_header = train_labels.columns.values
train_header[0] = "sample"
train_labels.columns = train_header
train_labels["sample"] = ["s" + str(x) for x in train_labels["sample"].values]


dir_out = os.path.join("data", "design", "mnist_fashion_train", "all")
if not os.path.exists(dir_out):
    os.makedirs(dir_out)

train_labels.to_csv(os.path.join(dir_out, "design.tsv"), header=True, index=False, sep="\t")

dir_out = os.path.join("data", "design", "mnist_fashion_test", "all")
if not os.path.exists(dir_out):
    os.makedirs(dir_out)

test_labels.to_csv(os.path.join(dir_out, "design.tsv"), header=True, index=False, sep="\t")













dir_out = os.path.join("data", "mnist_fashion", "train")
train = pd.read_csv(os.path.join(dir_out, "tpm.xls"), sep="\t", index_col=0 )
dir_out = os.path.join("data", "mnist_fashion", "test")
test = pd.read_csv(os.path.join(dir_out, "tpm.xls"), sep="\t", index_col=0)

dir_out = os.path.join("data", "design", "mnist_fashion_test", "all")
test_labels =  pd.read_csv(os.path.join(dir_out, "design.tsv"), sep="\t", index_col=0)
dir_out = os.path.join("data", "design", "mnist_fashion_train", "all")
train_labels =  pd.read_csv(os.path.join(dir_out, "design.tsv"), sep="\t", index_col=0)

def run_train_label(train, train_labels, test, test_labels, features_ids_keep=None, top=None):
    if features_ids_keep is None:
        dir_out = os.path.join("data", "res", "mnist_fashion", "all")
    else:
        dir_out = os.path.join("data", "res", "mnist_fashion", "top"+str(top))
        train = train.loc[features_ids_keep]
        test = test.loc[features_ids_keep]

    if not os.path.exists(dir_out):
        os.makedirs(dir_out)

    train_labels.sort_values(by=["subgroup_i"], inplace=True)
    train = train.reindex(train_labels.index.tolist(), axis=1)
    test_labels.sort_values(by=["subgroup_i"], inplace=True)
    test = test.reindex(test_labels.index.tolist(), axis=1)
    nneighbors = 10
    train_embedding = umap.UMAP(
        n_neighbors=nneighbors,
        min_dist=0.0,
        n_components=2,
        random_state=42,
        n_epochs = 2000
    ).fit_transform(train.T, y=train_labels["subgroup_i"])

    train_labels["x"] = train_embedding[:, 0]
    train_labels["y"] = train_embedding[:, 1]
    train_labels["type"] = "train"

    labels = hdbscan.HDBSCAN(
        min_samples=int(math.floor(1000/4)),
        min_cluster_size=int(math.floor(1000/2))
    ).fit_predict(train_embedding)
    clustered = (labels >= 0)
    p_clustered = round(np.sum(clustered) / train_labels.shape[0], 2)
    full_ARS = round(adjusted_rand_score(train_labels["subgroup_i"], labels), 2)
    full_MIS = round(adjusted_mutual_info_score(train_labels["subgroup_i"], labels, average_method='arithmetic'), 2)
    ARS = round(adjusted_rand_score(train_labels["subgroup_i"][clustered], labels[clustered]), 2)
    MIS = round(adjusted_mutual_info_score(train_labels["subgroup_i"][clustered], labels[clustered], average_method='arithmetic'), 2)

    sns.set(rc={'figure.figsize':(11.7,8.27)})
    filled_markers = ('o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X')

    sns_plot = sns.scatterplot(
        x="x",
        y="y",
        style="subgroup",
        hue="subgroup",
        palette=sns.color_palette("Paired", train_labels["subgroup"].nunique()),
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
        min_dist=0.0,
        n_components=2,
        random_state=42,
        n_epochs = 2000
    ).fit(train.T, np.array(train_labels["subgroup_i"]))

    test_embedding = mapper.transform(test.T)

    test_labels["x"] = test_embedding[:, 0]
    test_labels["y"] = test_embedding[:, 1]
    test_labels["pred"] = test_embedding[:, 0]
    test_labels["type"] = "test"

    #both_labels = pd.concat([train_labels, test_labels])

    labels = hdbscan.HDBSCAN(
        min_samples=int(math.floor(1000/4)),
        min_cluster_size=int(math.floor(1000/2))
    ).fit_predict(test_embedding)
    clustered = (labels >= 0)
    p_clustered = round(np.sum(clustered) / test_labels.shape[0], 2)
    full_ARS = round(adjusted_rand_score(test_labels["subgroup_i"], labels), 2)
    full_MIS = round(adjusted_mutual_info_score(test_labels["subgroup_i"], labels, average_method='arithmetic'), 2)
    ARS = round(adjusted_rand_score(test_labels["subgroup_i"][clustered], labels[clustered]), 2)
    MIS = round(adjusted_mutual_info_score(test_labels["subgroup_i"][clustered], labels[clustered], average_method='arithmetic'), 2)



    from sklearn.metrics.pairwise import euclidean_distances

    all_dist= euclidean_distances(test_embedding, train_embedding)
    print(all_dist)
    for i in range(0,len(test_labels["pred"])):
        selected_subg = train_labels["subgroup"].iloc[np.argpartition(all_dist[i,:], nneighbors)[:nneighbors]]
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
    for subg in set(test_labels["subgroup"]):
        tp = sum(test_labels["pred"][test_labels["subgroup"] == subg] == subg)
        fp = sum(test_labels["pred"][test_labels["subgroup"] != subg] == subg)
        tn = sum(test_labels["pred"][test_labels["subgroup"] != subg] != subg)
        fn = sum(test_labels["pred"][test_labels["subgroup"] == subg] != subg)

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
        style="subgroup",
        hue="subgroup",
        palette=sns.color_palette("Paired", test_labels["subgroup"].nunique()),
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
    path_file = os.path.join(path_file, "train", "tpm", file_name)
    df_diff = pd.read_csv(path_file, sep="\t")
    df_diff["abs_L2FC"] = df_diff.l2fc.abs()
    df_diff = df_diff.sort_values(["kernel_mcc", "abs_L2FC"], ascending=False)

    return(df_diff)



run_train_label(train, train_labels, test, test_labels)


for top in [10,20,30,40,50,60,70]:
    features_ids_keep = []
    path_design="/u/eaudemard/project/epcy_paper/data/design/mnist_fashion_train"
    for design in ["6000_Ankle_boot", "6000_Bag", "6000_Coat", "6000_Dress", "6000_Pullover", "6000_Sandal", "6000_Shirt", "6000_Sneaker", "6000_T-shirt", "6000_Trouser"]:
        print(design)
        dir_design = os.path.join(path_design, design)
        df_diff = read_diff_table("prediction_capability.xls", dir_design)
        features_ids_keep += df_diff["id"][:top].tolist()

    run_train_label(train, train_labels, test, test_labels, features_ids_keep, top)
