import os
import sys
import re
import math
import hdbscan

from collections import defaultdict

import numpy as np
import pandas as pd
import seaborn as sns
import umap

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm

# from scipy.cluster.hierarchy import fclusterdata
from scipy.spatial.distance import pdist, squareform
from scipy.stats import mannwhitneyu

from sklearn.metrics import silhouette_score, silhouette_samples
from sklearn.metrics import f1_score, precision_score, make_scorer
from sklearn.cluster import AffinityPropagation, KMeans
from sklearn import linear_model
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.model_selection import learning_curve
from sklearn.model_selection import ShuffleSplit
from sklearn.naive_bayes import GaussianNB
from sklearn.svm import SVC

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'DejaVu Sans'
mpl.rcParams['pdf.fonttype'] = 42

np.random.seed(42)


def create_LOG(args):
    # creates a log file
    LOG = str(sys.date)
    return LOG


def select_top(args, df_diff, top):
    if (len(df_diff.index)) > top:
        df_diff = df_diff[:top]

    return(df_diff)


def select_epcy(args, df_diff, exp_genes, method, df_ext):
    # print("Read epcy")

    if exp_genes is not None:
        df_diff = df_diff.loc[df_diff["ID"].isin(exp_genes)]

    df_diff = df_diff.loc[df_diff["KERNEL_MCC"] >= args.MCC]
    df_diff = df_diff.loc[abs(df_diff["L2FC"]) >= args.LOG_FC]

    if df_ext is not None:
        df_diff = df_diff.merge(df_ext, on="ID")
        df_diff.L2FC = df_diff.L2FC_ext

    df_diff["abs_L2FC"] = df_diff.L2FC.abs()

    if method == "epcy" or method == "epcy_bs" or method == "epcy_kernel" or method == "epcy_bagging":
        df_diff = df_diff.sort_values(
            ["KERNEL_MCC", "abs_L2FC"], ascending=False
        )
    elif method == "epcy_auc":
        df_diff = df_diff.sort_values("AUC", ascending=False)
    elif method == "epcy_normal":
        df_diff = df_diff.sort_values(
            ["NORMAL_MCC", "abs_L2FC"], ascending=False
        )
    elif method == "epcy_log2FC":
        df_diff = df_diff.sort_values("abs_L2FC", ascending=False)
    else:
        sys.stderr.write('ERROR: in utils/other.py function select_epcy.\n')
        sys.stderr.write(method)
        sys.exit()

    return(df_diff)


def select_mast(args, df_diff, exp_genes, padj):
    # print("Read deseq")
    if exp_genes is not None:
        df_diff = df_diff.loc[df_diff["ID"].isin(exp_genes)]
    df_diff["pval"] = df_diff["pval"].fillna(1)
    df_diff = df_diff.loc[df_diff["pval"] <= padj]

    df_diff = df_diff.sort_values(["pval"], ascending=True)

    # df_diff = select_top(args, df_diff, top)

    return(df_diff)


def select_deseq(args, df_diff, exp_genes, padj):
    # print("Read deseq")
    if exp_genes is not None:
        df_diff = df_diff.loc[df_diff["ID"].isin(exp_genes)]
    df_diff["padj"] = df_diff["padj"].fillna(1)
    df_diff = df_diff.loc[df_diff["padj"] <= padj]
    df_diff = df_diff.loc[abs(df_diff["log2FoldChange"]) >= args.LOG_FC]
    df_diff["abs_L2FC"] = df_diff.log2FoldChange.abs()
    df_diff = df_diff.sort_values(["padj", "abs_L2FC"], ascending=[True, False])

    return(df_diff)


def select_edger(args, df_diff, exp_genes, fdr):
    # print("Read edger")
    if exp_genes is not None:
        df_diff = df_diff.loc[df_diff["ID"].isin(exp_genes)]
    df_diff = df_diff.loc[df_diff["FDR"] <= fdr]
    df_diff = df_diff.loc[abs(df_diff["logFC"]) >= args.LOG_FC]
    df_diff["abs_L2FC"] = df_diff.logFC.abs()
    df_diff["logFC"] = - df_diff["logFC"]
    df_diff = df_diff.sort_values(["FDR", "abs_L2FC"], ascending=[True, False])

    return(df_diff)


def select_limma(args, df_diff, exp_genes, padj):
    # print("Read limma")
    if exp_genes is not None:
        df_diff = df_diff.loc[df_diff["ID"].isin(exp_genes)]
    df_diff = df_diff.loc[df_diff["adj.P.Val"] <= padj]
    df_diff = df_diff.loc[abs(df_diff["logFC"]) >= args.LOG_FC]
    df_diff["abs_L2FC"] = df_diff.logFC.abs()
    df_diff["logFC"] = - df_diff["logFC"]
    df_diff = df_diff.sort_values(["adj.P.Val", "abs_L2FC"], ascending=[True, False])

    # df_diff = select_top(args, df_diff, top)
    return(df_diff)


def read_diff_table(args, file_name, method, path_dir, padj,
                    exp_genes=None, df_ext=None):

    path_file = path_dir

    # TODO: add parameters for STAR and readcounts
    if "bagging" in method:
        path_file = os.path.join(path_file, args.QUANT, args.TYPE_QUANT + "_bagging", file_name)
    else:
        path_file = os.path.join(path_file, args.QUANT, args.TYPE_QUANT, file_name)

    df_diff = pd.read_csv(path_file, sep="\t")

    if "epcy" in method:
        df_diff.rename(str.upper, axis='columns', inplace=True)
        df_diff = select_epcy(args, df_diff, exp_genes, method, df_ext)
    if method == "deseq2":
        df_diff = select_deseq(args, df_diff, exp_genes, padj)
    if method == "deseq2_pvalue":
        df_diff = select_deseq(args, df_diff, exp_genes, padj, by_pvalue=True)
    if method == "edger":
        df_diff = select_edger(args, df_diff, exp_genes, padj)
    if method == "edger_pvalue":
        df_diff = select_edger(args, df_diff, exp_genes, padj, by_pvalue=True)
    if method == "limma":
        df_diff = select_limma(args, df_diff, exp_genes, padj)
    if method == "limma_pvalue":
        df_diff = select_limma(args, df_diff, exp_genes, padj, by_pvalue=True)
    if method == "mast":
        df_diff = select_mast(args, df_diff, exp_genes, padj)

    return(df_diff)


def read_sample_pred_table(args, file_name, df_diff, dataset=None):
    path_file = args
    if type(args) is not str:
        path_file = os.path.join(args.PATH)
    if dataset is not None:
        path_file = os.path.join(args.PATH, dataset)

    path_file = os.path.join(path_file, file_name)

    df_pred = pd.read_table(path_file)
    df_pred = df_pred.merge(df_diff[["gene_id", "ID"]],
                            left_on="IDS", right_on="gene_id")
    df_pred = df_pred.set_index("ID")
    df_pred = df_pred.drop("IDS", 1)
    df_pred = df_pred.drop("gene_id", 1)

    df_pred = df_pred.replace(2, 1)
    df_pred = df_pred.replace(3, -1)
    df_pred = df_pred.replace(4, -1)

    for i, row in df_diff.iterrows():
        df_pred.loc[row["ID"]] = df_pred.loc[row["ID"]] * 10**row["KERNEL_MCC"]

    return(df_pred)


def get_design(args, subg_query,  path=None, dataset=None,
               file_name="design.tsv", replace=True):
    design = path
    if path is None:
        design = os.path.join(args.PATH)
    if dataset is not None:
        design = os.path.join(args.PATH, dataset)

    design = os.path.join(design, file_name)
    design = pd.read_csv(design, sep="\t")
    design["subgroup_str"] = design[args.SUBGROUP]
    if replace:
        design[args.SUBGROUP] = [1 if condition == subg_query else 0 for condition in design[args.SUBGROUP]]
    design = design.sort_values(by=[args.SUBGROUP, 'sample'],
                                ascending=[False, True])

    return(design)


def get_exp(args, file_name, df_design=None):
    if args.BIOTYPE is not None:
        df_biotype = pd.read_csv(args.BF, sep="\t")
        selected_biotype = args.BIOTYPE.split(",")
        df_biotype = df_biotype.loc[df_biotype["gene_biotype"].isin(selected_biotype)]

    df_exp = pd.read_csv(file_name, sep="\t")

    if args.BIOTYPE is not None:
        df_exp = df_exp.loc[df_exp["ID"].isin(df_biotype["ensembl_gene_id"])]

    df_exp = df_exp.set_index("ID")
    if df_design is not None:
        df_exp = df_exp[df_design["sample"]]

    if args.CPM:
        f_norm = 1e6 / df_exp.sum()
        df_exp = df_exp * f_norm

    df_exp = np.log2(df_exp + 1)
    return(df_exp)


def get_silhouette_kmeans(df, metric, num_cluster):
    X = df.T.values
    kmeans = KMeans(n_clusters=num_cluster)
    labels = kmeans.fit_predict(X)

    silh_score = silhouette_score(X, labels, metric=metric)

    return(silh_score, labels)


def get_silhouette_af(df, metric):
    X = df.T.values
    af = AffinityPropagation().fit(X)
    cluster_centers_indices = af.cluster_centers_indices_
    num_cluster_ = len(cluster_centers_indices)
    labels = af.labels_

    silh_score = silhouette_score(X, labels, metric=metric)

    return(silh_score, num_cluster_, labels)


def get_silhouette_graph(df, labels, metric, path_out, method, top,
                         silhouette_avg, num_cluster):
    X = df.T.values
    cluster_labels = labels

    # Create a subplot with 1 row and 2 columns
    fig, (ax1, ax2) = plt.subplots(1, 2)
    fig.set_size_inches(18, 7)

    # The 1st subplot is the silhouette plot
    # The silhouette coefficient can range from -1, 1 but in this example all
    # lie within [-0.1, 1]
    ax1.set_xlim([-0.1, 1])
    # The (num_cluster+1)*10 is for inserting blank space between silhouette
    # plots of individual clusters, to demarcate them clearly.
    ax1.set_ylim([0, len(X) + (num_cluster + 1) * 10])

    # Compute the silhouette scores for each sample
    sample_silhouette_values = silhouette_samples(X, cluster_labels,
                                                  metric=metric)
    y_lower = 10
    for i in range(num_cluster):
        # Aggregate the silhouette scores for samples belonging to
        # cluster i, and sort them
        ith_cluster_silhouette_values = \
            sample_silhouette_values[cluster_labels == i]

        ith_cluster_silhouette_values.sort()

        size_cluster_i = ith_cluster_silhouette_values.shape[0]
        y_upper = y_lower + size_cluster_i

        color = cm.nipy_spectral(float(i) / num_cluster)
        ax1.fill_betweenx(np.arange(y_lower, y_upper),
                          0, ith_cluster_silhouette_values,
                          facecolor=color, edgecolor=color, alpha=0.7)

        # Label the silhouette plots with their cluster numbers at the middle
        ax1.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))

        # Compute the new y_lower for next plot
        y_lower = y_upper + 10  # 10 for the 0 samples

    ax1.set_title("The silhouette plot for the various clusters.")
    ax1.set_xlabel("The silhouette coefficient values")
    ax1.set_ylabel("Cluster label")

    # The vertical line for average silhouette score of all the values
    ax1.axvline(x=silhouette_avg, color="red", linestyle="--")

    ax1.set_yticks([])  # Clear the yaxis labels / ticks
    ax1.set_xticks([-0.1, 0, 0.2, 0.4, 0.6, 0.8, 1])

    # 2nd Plot showing the actual clusters formed
    colors = cm.nipy_spectral(cluster_labels.astype(float) / num_cluster)
    ax2.scatter(X[:, 0], X[:, 1], marker='.', s=30, lw=0, alpha=0.7,
                c=colors, edgecolor='k')

    # Labeling the clusters
    # centers = clusterer.cluster_centers_
    # Draw white circles at cluster centers
    # ax2.scatter(centers[:, 0], centers[:, 1], marker='o',
    #            c="white", alpha=1, s=200, edgecolor='k')

    # for i, c in enumerate(centers):
    #    ax2.scatter(c[0], c[1], marker='$%d$' % i, alpha=1,
    #                s=50, edgecolor='k')

    ax2.set_title("The visualization of the clustered data.")
    ax2.set_xlabel("Feature space for the 1st feature")
    ax2.set_ylabel("Feature space for the 2nd feature")

    plt.suptitle(("Silhouette analysis for KMeans clustering on sample data "
                  "with num_cluster = %d" % num_cluster),
                 fontsize=14, fontweight='bold')

    fig_dir = os.path.join(path_out, method, str(num_cluster))
    if not os.path.exists(fig_dir):
        os.makedirs(fig_dir)
    file_out = os.path.join(fig_dir, str(top) + "_silhouette.pdf")
    plt.savefig(file_out)
    plt.close()


def get_silhouette(df, df_design, metric):
    silh_score = silhouette_score(df.T, df_design[args.SUBGROUP], metric=metric)

    return(silh_score)


def get_ct_closest(query_xy, query_label, ref_xy, ref_label, n_neighbors):
    all_dist = euclidean_distances(xy1, xy2)

    # ct = [tp, fp, fn, tn]
    ct = [0, 0, 0, 0]
    for i in range(0, len(known_label)):
        closest_xy2 = ref_label[np.argpartition(all_dist[i,:], n_neighbors)[:n_neighbors]]
        closest_xy2 = closest_xy2.tolist()

        label_pred = max(closest_xy2, key=closest_xy2.count)
        if query_label[i] == 1 and label_pred == 1:
            ct[0] += 1
        if query_label[i] == 0 and label_pred == 1:
            ct[1] += 1
        if query_label[i] == 1 and label_pred == 0:
            ct[2] += 1
        if query_label[i] == 0 and label_pred == 0:
            ct[3] += 1

    return(ct)


def get_ct_probas(probas, labels):
    # ct = [tp, fp, fn, tn]
    ct = [0, 0, 0, 0]
    ids_query = np.where(labels == 1)
    ids_ref = np.where(labels == 0)

    # compute mcc
    ct[0] = np.sum(probas[ids_query] > 0.5)
    ct[1] = np.sum(probas[ids_query] <= 0.5)
    ct[3] = np.sum(probas[ids_ref] <= 0.5)
    ct[2] = np.sum(probas[ids_ref] > 0.5)

    return(ct)


def get_mcc(ct):
    # ct = [tp, fp, fn, tn]
    d1 = ct[0] + ct[1]
    d2 = ct[0] + ct[2]
    d3 = ct[3] + ct[1]
    d4 = ct[3] + ct[2]

    n1 = ct[0] * ct[3]
    n2 = ct[1] * ct[2]

    mcc = (n1 - n2) / math.sqrt(d1 * d2 * d3 * d4) if d1!=0 and d2!=0 and d3!=0 and d4!=0 else 0
    if mcc < 0:
        mcc = -mcc
    return(mcc)


def auc_u_test(x, ids_query, ids_ref):
    # Don't need to sort when we use mannwhitneyu in auc_u_test

    (u_value, p_value) = mannwhitneyu(x[ids_query], x[ids_ref],
                                      alternative="two-sided")
    auc = u_value / (len(x[ids_query]) * len(x[ids_ref]))

    if auc < 0.5:
        auc = 1 - auc
    return(auc, p_value)


def get_sample_dist(df, df_design, metric):
    samples_dist = pdist(df.T, metric=metric)
    samples_dist = squareform(samples_dist)

    num_sample = len(df_design[args.SUBGROUP])
    sample_query = df_design["sample"].loc[df_design[args.SUBGROUP] == 1]

    close_perc = [0.0, 0.0]

    for group in [0, 1]:
        sample_selected = df_design["sample"].loc[df_design[args.SUBGROUP] == group]
        num_sample_selected = len(sample_selected)
        ind_exp_selected = np.where(df.columns.isin(sample_selected))[0]

        sum_same_group = 0
        for index in ind_exp_selected:
            sample_dist = [samples_dist[index, j] for j in range(0, num_sample)]
            cutoff_dist = sorted(sample_dist)[num_sample_selected]
            closest_samples = np.where(sample_dist <= cutoff_dist)
            closest_samples = df.columns.values[closest_samples]

            num_same_group = np.sum(df_design[args.SUBGROUP].loc[df_design["sample"].isin(closest_samples)] == group)
            num_same_group = np.min([num_same_group, num_sample_selected])
            sum_same_group += num_same_group - 1

        close_perc[group] = sum_same_group / (num_sample_selected * (num_sample_selected - 1))

    return(close_perc, samples_dist)


def display_cluster(args, df, df_design, path_out, method, metric, top,
                    num_cluster, max_exp):
    col_pal = [
        mpl.colors.hex2color('#E62528'),
        mpl.colors.hex2color('#0D6BAC')
    ]
    max_exp = max_exp * 3/4
    samples_palette = dict(zip(df_design[args.SUBGROUP].unique(),
                               [col_pal[0], col_pal[1]]))
    samples_colors = df_design[args.SUBGROUP].map(samples_palette)
    samples_colors = dict(zip(df_design["sample"], samples_colors.get_values()))
    samples_colors = pd.Series(df.columns, index=df.columns).map(samples_colors)

    # default seaborn size is (10, 10)
    figsize = (10, np.max([10, len(df.index)/4]))
    figsize = (10, 10)
    display_label = False
    if args.SCALED:
        figsize = (len(df.columns)/4, len(df.index)/4)
        display_label = True

    sns_plot = sns.clustermap(
        df, figsize=figsize, linewidths=0, metric=metric,
        col_colors=samples_colors, # row_colors=mcc_colors, #cmap="viridis_r",
        xticklabels=display_label, yticklabels=False,
        vmin=0, vmax=max_exp, cmap="Greys"
    )
    sns_plot.fig.suptitle(method)
    # plt.plot(sns_plot)
    fig_dir = os.path.join(path_out, method)
    if not os.path.exists(fig_dir):
        os.makedirs(fig_dir)
    file_out = os.path.join(fig_dir, str(top) + "_heatmap.pdf")
    sns_plot.savefig(file_out)
    plt.close()


def get_dist_by_group(df, df_design, samples_dist, method, top, metric):
    num_sample = len(df_design[args.SUBGROUP])

    dict_dist = defaultdict(list)
    dist_diff = []
    dist_diff_sample = []
    for i in range(0, num_sample):
        if i != num_sample - 1:
            group_i = df_design[args.SUBGROUP].loc[df_design["sample"] == df.columns.values[i]].values[0]
            for j in range(i + 1, num_sample):
                group_j = df_design[args.SUBGROUP].loc[df_design["sample"] == df.columns.values[j]].values[0]
                if group_i == group_j:
                    if group_i == 1:
                        dict_dist[args.SUBGROUP].append("Query")
                    else:
                        dict_dist[args.SUBGROUP].append("Ref")
                else:
                    dict_dist[args.SUBGROUP].append("diff")
                    dist_diff.append(samples_dist[i, j])
                    dist_diff_sample.append([df.columns.values[i], df.columns.values[j]])
                dict_dist['dist'].append(samples_dist[i, j])
                dict_dist['metric'].append(metric)
                dict_dist['method'].append(method)
                dict_dist['top'].append(top)
    return(dict_dist, dist_diff, dist_diff_sample)


def run_kmeans(XTEST, YTEST, design, method, dataset, top):
    silh_score, labels = get_silhouette_kmeans(XTEST, "euclidean", 2)
    YTEST = np.array(YTEST)
    labels = np.array(labels)
    pred = [labels[i] == YTEST[i] for i in range(len(YTEST))]
    pred = np.sum(pred) / len(YTEST)
    if pred < 0.5:
        pred = 1 - pred

    id_dataset = re.findall(r'\d+', dataset)[0]

    dict_kmean = defaultdict(list)
    dict_kmean['method'].append(method)
    dict_kmean['top'].append(top)
    dict_kmean['pred'].append(pred)
    dict_kmean['dataset'].append(id_dataset)
    dict_kmean['design'].append(design)
    dict_kmean['learn'].append("k_mean")

    return(dict_kmean)


def run_model_by_feature(model, XTRAIN, YTRAIN, XTEST, YTEST, id_dataset,
                         param, learn_meth=""):
    for ids, (feature_name, feature_data) in enumerate(XTRAIN.T.iteritems()):
        model.fit(XTEST.T[:, [ids]], YTRAIN)

        proba_test = model.predict_proba(XTEST.T[:, [ids]])
        proba_train = model.predict_proba(XTRAIN.T[:, [ids]])

        ct_test = get_ct_probas(proba_test)
        mcc_test = get_mcc(ct_test)

        ct_train = get_ct_probas(proba_train)
        mcc_train = get_mcc(ct_train)

    dict_res = defaultdict(list)
    dict_res['ct_train'].append(ct_train)
    dict_res['mcc_train'].append(mcc_train)
    dict_res['ct_test'].append(ct_test)
    dict_res['mcc_test'].append(mcc_test)
    dict_res['id_dataset'].append(id_dataset)
    dict_res['learn_meth'].append(learn_meth)
    for key in param.keys():
        dict_res[key] = dict_res[key] + param[key]

    return(dict_res)


def run_rand_forest(XTRAIN, YTRAIN, XTEST, YTEST, design, method, dataset,
                    top, pvalue, n_estimators=100,
                    class_weight="balanced_subsample"):
    rf = RandomForestClassifier(n_estimators=n_estimators,
                                class_weight=class_weight)
    rf.fit(XTRAIN.T, YTRAIN)

    proba_test = rf.predict_proba(XTEST.T)

    id_dataset = re.findall(r'\d+', dataset)[0]
    id_dataset = int(id_dataset)

    dict_rf = defaultdict(list)
    dict_rf['method'].append(method)
    dict_rf['top'].append(top)
    dict_rf['label'].append(YTEST.values)
    dict_rf['proba'].append(proba_test)
    dict_rf['dataset'].append(id_dataset)
    dict_rf['design'].append(design)
    dict_rf['learn'].append("R-Forest")
    dict_LR['pvalue'].append(pvalue)

    return(dict_rf)


def run_LR(XTRAIN, YTRAIN, XTEST, YTEST, design, method, dataset, top, pvalue,
           solver="liblinear", C=1e10, penalty='l2',
           max_iter=100, class_weight='balanced'):

    # large C ==> no regularization
    LR = linear_model.LogisticRegression(
        solver="liblinear", C=1e10, penalty='l2', max_iter=1000,
        class_weight='balanced')
    LR.fit(XTRAIN.T, YTRAIN)

    # tpm_test = np.reshape(df_exp_test, (-1, 1))
    proba_test = LR.predict_proba(XTEST.T)
    proba_train = LR.predict_proba(XTRAIN.T)

    id_dataset = re.findall(r'\d+', dataset)[0]
    id_dataset = int(id_dataset)

    dict_LR = defaultdict(list)
    dict_LR['method'].append(method)
    dict_LR['top'].append(top)
    dict_LR['label'].append(YTEST.values)
    dict_LR['proba'].append(proba_test)
    dict_LR['label_train'].append(YTRAIN.values)
    dict_LR['proba_train'].append(proba_train)
    dict_LR['dataset'].append(id_dataset)
    dict_LR['design'].append(design)
    dict_LR['learn'].append("L-Reg")
    dict_LR['pvalue'].append(pvalue)

    return(dict_LR)


def run_umap(XTRAIN, YTRAIN, XTEST, YTEST, id_dataset, param,
             n_neighbors=15, n_epochs=1000, random_state=42):

    mapper = umap.UMAP(
        n_neighbors=n_neighbors,
        min_dist=0.2,
        n_components=2,
        random_state=random_state,
        n_epochs=n_epochs
    ).fit(XTRAIN.T, YTRAIN)

    train_embedding = mapper.transform(XTRAIN.T)
    test_embedding = mapper.transform(XTEST.T)

    test_labels["x"] = test_embedding[:, 0]
    test_labels["y"] = test_embedding[:, 1]

    # ct = [tp, fp, fn, tn]
    ct_test = get_ct_closest(test_embedding, YTEST, train_embedding, YTRAIN,
                             n_neighbors)
    mcc_test = get_mcc(ct_test)

    ct_train = get_ct_closest(train_embedding, YTRAIN, train_embedding, YTRAIN,
                              n_neighbors)
    mcc_train = get_mcc(ct_train)

    dict_res = defaultdict(list)
    dict_res['ct_train'].append(ct_train)
    dict_res['mcc_train'].append(mcc_train)
    dict_res['ct_test'].append(ct_test)
    dict_res['mcc_test'].append(mcc_test)
    dict_res['id_dataset'].append(id_dataset)
    dict_res['learn_meth'].append("UMAP")
    for key in param.keys():
        dict_res[key] = dict_res[key] + param[key]

    return(dict_res)


def get_exp_cluster(args, df_exp, df_design, df_diff, path_out, method,
                    design, top, num_cluster=2, max_exp=20, metric='euclidean',
                    type='heatmap'):
    df_exp = df_exp.loc[df_exp.index.isin(df_diff["ID"][:top])]
    return(get_cluster(args, df_exp, df_design, path_out, method, design, top,
                       num_cluster, max_exp, metric=metric, type=type))


def get_cluster(args, df_exp, df_design, path_out, method, design, top,
                num_cluster, max_exp, metric='euclidean', type='heatmap'):

    if type == "heatmap":
        display_cluster(args, df_exp, df_design, path_out, method, metric,
                        top, num_cluster, max_exp)
        (silh_score, labels) = get_silhouette_kmeans(df, metric, num_cluster)
        # get_silhouette_graph(df, labels, metric, path_out, method, top, silh_score, num_cluster)

        # pred = [labels[i] == df_design["group"][i] for i in range(len(df_design["group"]))]
        # pred = np.sum(pred) / len(df_design["group"])
        # if pred < 0.5:
        #    pred = 1 - pred

        return(silh_score)

    if type == "umap":
        return(cluster_umap(args, df_exp, df_design, path_out, method, design,
                            top, metric='euclidean'))


def cluster_umap(args, df_exp, df_design, path_out, method, design, top,
                 metric='euclidean', r_state=42, num_epochs=2000,
                 num_neighbors=15):

    filled_markers = ('o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H',
                      'D', 'd', 'P', 'X')

    # df_exp = df_exp.loc[df_exp.index.isin(df_diff["ID"][:top])]
    num_query = len(np.where(df_design[args.SUBGROUP] == 1)[0])
    num_ref = len(np.where(df_design[args.SUBGROUP] == 0)[0])
    num_small = min(num_query, num_ref)
    if num_small == 0:
        num_small = 20

    # STEP1: display UMAP with subgroup
    reducer = umap.UMAP(
      n_neighbors=num_neighbors,
      min_dist=0.2,
      metric=metric,
      n_components=2,
      random_state=r_state,
      n_epochs=num_epochs)
    sample_embedding = reducer.fit_transform(df_exp.T)
    df_design["x"] = sample_embedding[:, 0]
    df_design["y"] = sample_embedding[:, 1]

    # STEP2: Find cluster
    clusterable_embedding = umap.UMAP(
        n_neighbors=num_neighbors,
        min_dist=0.0,
        metric=metric,
        n_components=2,
        random_state=r_state,
        n_epochs=num_epochs
    ).fit_transform(df_exp.T)

    labels = hdbscan.HDBSCAN(
        min_samples=int(math.floor(num_small/4)),
        min_cluster_size=int(math.floor(num_small/2))
    ).fit_predict(clusterable_embedding)
    clustered = (labels >= 0)

    # STEP3: compute metric
    p_clustered = round(np.sum(clustered) / df_exp.shape[1], 2)
    if p_clustered == 0:
        full_ARS = 0
        full_MIS = 0
        ARS = 0
        MIS = 0
    else:
        full_ARS = round(adjusted_rand_score(df_design[args.SUBGROUP], labels), 2)
        full_MIS = round(adjusted_mutual_info_score(df_design[args.SUBGROUP], labels, average_method='arithmetic'), 2)
        ARS = round(adjusted_rand_score(df_design[args.SUBGROUP][clustered], labels[clustered]), 2)
        MIS = round(adjusted_mutual_info_score(df_design[args.SUBGROUP][clustered], labels[clustered], average_method='arithmetic'), 2)

    dict_umap = defaultdict(list)
    dict_umap['full_ARS'].append(full_ARS)
    dict_umap['full_MIS'].append(full_MIS)
    dict_umap['ARS'].append(ARS)
    dict_umap['MIS'].append(MIS)
    dict_umap['p_clustered'].append(p_clustered)
    dict_umap['top'].append(top)
    dict_umap['method'].append(method)
    dict_umap['design'].append(design)

    # STEP3: display cluster

    sns.set(rc={'figure.figsize': (11.7, 8.27)})
    df_design["cluster"] = [str(x) for x in labels]

    for hue in ["cluster", "subgroup_str", "stranded", "project", "tissue", "blasts"]:

        if hue == "blasts":
            sns_plot = sns.scatterplot(
                x="x",
                y="y",
                style="subgroup_str",
                hue=hue,
                markers=filled_markers,
                data=df_design
            )
        else:
            sns_plot = sns.scatterplot(
                x="x",
                y="y",
                style="subgroup_str",
                hue=hue,
                palette=sns.color_palette("Paired", df_design[hue].nunique()),
                markers=filled_markers,
                data=df_design
            )
        # resize figure box to -> put the legend out of the figure
        box = sns_plot.get_position()
        sns_plot.set_position([box.x0, box.y0, box.width * 0.75, box.height])

        # Put a legend to the right side
        sns_plot.legend(loc='center right', bbox_to_anchor=(1.45, 0.5), ncol=2)

        sns_plot.set_title(method + ": " + str(full_ARS) + ", " + str(full_MIS) + ", " + str(p_clustered))
        fig_dir = os.path.join(path_out, method)
        if not os.path.exists(fig_dir):
            os.makedirs(fig_dir)
        file_out = os.path.join(fig_dir, str(top) + "_" + hue + "_umap." + args.EXT)
        sns_plot = sns_plot.get_figure()
        sns_plot.savefig(file_out)
        plt.close()

    file_out = os.path.join(fig_dir, str(top) + "_count_by_subg." + args.EXT)
    sns_plot = sns.countplot(
        x="subgroup_str",
        hue="cluster",
        data=df_design,
        palette=sns.color_palette("Paired", df_design["subgroup_str"].nunique())
    )
    sns_plot.legend(bbox_to_anchor=(1.01, 1), loc='upper left',
                    borderaxespad=0.)
    sns_plot = sns_plot.get_figure()
    sns_plot.savefig(file_out)
    plt.close()


    file_out = os.path.join(fig_dir, str(top) + "_count_by_cluster." + args.EXT)
    sns_plot = sns.countplot(
        x="cluster",
        hue="subgroup_str",
        data=df_design,
        palette=sns.color_palette("Paired", df_design["subgroup_str"].nunique())
    )
    sns_plot.legend(bbox_to_anchor=(1.01, 1), loc='upper left',
                    borderaxespad=0.)
    sns_plot = sns_plot.get_figure()
    sns_plot.savefig(file_out)
    plt.close()

    csv_out = os.path.join(fig_dir, str(top) + "_design.csv")
    df_design.to_csv(csv_out, index=False, sep="\t")
    df_design

    csv_out = os.path.join(fig_dir, str(top) + "_exp.csv")
    df_exp.to_csv(csv_out, index=True, sep="\t")

    return(dict_umap)


def run_train_label(args, method, train, train_labels, test, test_labels, top,
                    pvalue, p_train,
                    features_ids_keep=None, nneighbors=15, num_epochs=500):

    if features_ids_keep is None:
        dir_out = os.path.join(args.OUTDIR, "all")
    else:
        if top is not None:
            dir_out = os.path.join(args.OUTDIR, method, "top" + str(top), "ptrain" + str(p_train))
            train = train.loc[features_ids_keep]
            test = test.loc[features_ids_keep]

    if not os.path.exists(dir_out):
        os.makedirs(dir_out)

    train_embedding = umap.UMAP(
        n_neighbors=nneighbors,
        min_dist=0.2,
        n_components=2,
        random_state=42,
        n_epochs=num_epochs
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

    sns.set(rc={'figure.figsize': (11.7, 8.27)})
    filled_markers = ('o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H',
                      'D', 'd', 'P', 'X')
    color_palette = sns.color_palette("Paired",
                                      train_labels["subgroup_str"].nunique())

    sns_plot = sns.scatterplot(
        x="x",
        y="y",
        style="subgroup_str",
        hue="subgroup_str",
        palette=color_palette,
        markers=filled_markers,
        data=train_labels
    )
    box = sns_plot.get_position()
    sns_plot.set_position([box.x0, box.y0, box.width * 0.75, box.height])
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
        n_epochs=num_epochs
    ).fit(train.T, np.array(train_labels["subgroup"]))

    test_embedding = mapper.transform(test.T)

    test_labels["x"] = test_embedding[:, 0]
    test_labels["y"] = test_embedding[:, 1]
    test_labels["pred"] = test_embedding[:, 0]
    test_labels["type"] = "test"

    # both_labels = pd.concat([train_labels, test_labels])

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

    all_dist = euclidean_distances(test_embedding, train_embedding)
    for i in range(0, len(test_labels["pred"])):
        selected_subg = train_labels["subgroup_str"].iloc[np.argpartition(all_dist[i,:], nneighbors)[:nneighbors]]
        selected_subg = selected_subg.tolist()
        test_labels.loc[test_labels.iloc[i].name, "pred"] = max(selected_subg, key=selected_subg.count)

    dict_res_subg = defaultdict(list)
    dict_res_all = defaultdict(list)
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

        dict_res_subg['tp'].append(tp)
        dict_res_subg['tn'].append(tn)
        dict_res_subg['fn'].append(fn)
        dict_res_subg['fp'].append(fp)
        dict_res_subg['acc'].append(acc)
        dict_res_subg['mcc'].append(mcc)
        dict_res_subg['subg'].append(subg)
        dict_res_subg['top'].append(top)
        dict_res_subg['method'].append(method)
        dict_res_subg['pvalue'].append(pvalue)
        dict_res_subg['p_train'].append(p_train)

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

    dict_res_all['tp'].append(all_tp)
    dict_res_all['tn'].append(all_tn)
    dict_res_all['fp'].append(all_fp)
    dict_res_all['fn'].append(all_fn)
    dict_res_all['acc'].append(all_acc)
    dict_res_all['mcc'].append(all_mcc)
    dict_res_all['full_ARS'].append(full_ARS)
    dict_res_all['full_MIS'].append(full_MIS)
    dict_res_all['p_clustered'].append(p_clustered)
    dict_res_all['num_feature'].append(test.shape[0])
    dict_res_all['top'].append(top)
    dict_res_all['method'].append(method)
    dict_res_all['pvalue'].append(pvalue)
    dict_res_all['p_train'].append(p_train)

    color_palette = sns.color_palette("Paired", test_labels["subgroup_str"].nunique())
    sns_plot = sns.scatterplot(
        x="x",
        y="y",
        style="subgroup_str",
        hue="subgroup_str",
        palette=color_palette,
        markers=filled_markers,
        data=test_labels
    )
    box = sns_plot.get_position()
    sns_plot.set_position([box.x0, box.y0, box.width * 0.75, box.height])
    sns_plot.legend(loc='center right', bbox_to_anchor=(1.45, 0.5), ncol=2)
    sns_plot.set_title(
        method +
        " ACC=" + str(all_acc) +
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

    return(dict_res_subg, dict_res_all)
