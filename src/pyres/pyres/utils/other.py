import os
import sys
import re
import math

from collections import defaultdict

import numpy as np
import pandas as pd
import seaborn as sns

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm

#from scipy.cluster.hierarchy import fclusterdata
from scipy.spatial.distance import pdist, squareform
from scipy.stats import mannwhitneyu

from sklearn.metrics import silhouette_score, silhouette_samples
from sklearn.cluster import AffinityPropagation, KMeans
from sklearn import linear_model
from sklearn.ensemble import RandomForestClassifier

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'DejaVu Sans'
mpl.rcParams['pdf.fonttype'] = 42

def create_LOG(args):
    # creates a log file
    LOG = str(sys.date)
    return LOG

def select_top(args, df_diff, top):
    if (len(df_diff.index)) > top:
        df_diff = df_diff[:top]

    return(df_diff)

def select_epcy(args, df_diff, exp_genes, method, df_ext):
    #print("Read epcy")

    df_diff = df_diff.loc[df_diff["ID"].isin(exp_genes)]
    df_diff = df_diff.loc[df_diff["KERNEL_MCC"] >= args.MCC]
    df_diff = df_diff.loc[abs(df_diff["L2FC"]) >= args.LOG_FC]
    if df_ext is not None:
        df_diff = df_diff.merge(df_ext, on="ID")
        df_diff.L2FC = df_diff.L2FC_ext

    df_diff["abs_L2FC"] = df_diff.L2FC.abs()


    if method == "epcy" or method == "epcy_bs" or method == "epcy_kernel":
        df_diff = df_diff.sort_values(["KERNEL_MCC", "abs_L2FC"], ascending=False)
    elif method == "epcy_auc":
        df_diff = df_diff.sort_values("AUC", ascending=False)
    elif method == "epcy_normal":
        df_diff = df_diff.sort_values(["NORMAL_MCC", "abs_L2FC"], ascending=False)
    elif method == "epcy_log2FC":
        df_diff = df_diff.sort_values("abs_L2FC", ascending=False)
    else:
        sys.stderr.write('ERROR: in utils/other.py function select_epcy.\n')
        sys.stderr.write(method)
        sys.exit()

    return(df_diff)

def select_deseq(args, df_diff, exp_genes):
    #print("Read deseq")
    df_diff = df_diff.loc[df_diff["ID"].isin(exp_genes)]
    df_diff = df_diff.loc[df_diff["pvalue"] <= args.PVALUE]
    df_diff = df_diff.loc[abs(df_diff["log2FoldChange"]) >= args.LOG_FC]
    df_diff = df_diff.reindex(df_diff.log2FoldChange.abs().sort_values(ascending=False).index)

    # df_diff = select_top(args, df_diff, top)

    return(df_diff)

def select_edger(args, df_diff, exp_genes):
    #print("Read edger")
    df_diff = df_diff.loc[df_diff["ID"].isin(exp_genes)]
    df_diff = df_diff.loc[df_diff["PValue"] <= args.PVALUE]
    df_diff = df_diff.loc[abs(df_diff["logFC"]) >= args.LOG_FC]
    df_diff = df_diff.reindex(df_diff.logFC.abs().sort_values(ascending=False).index)

    # df_diff = select_top(args, df_diff, top)

    return(df_diff)

def select_limma(args, df_diff, exp_genes):
    #print("Read limma")
    df_diff = df_diff.loc[df_diff["ID"].isin(exp_genes)]
    df_diff = df_diff.loc[df_diff["P.Value"] <= args.PVALUE]
    df_diff = df_diff.loc[abs(df_diff["logFC"]) >= args.LOG_FC]
    df_diff = df_diff.reindex(df_diff.logFC.abs().sort_values(ascending=False).index)

    # df_diff = select_top(args, df_diff, top)

    return(df_diff)

def read_diff_table(args, file_name, method, exp_genes, path_dir=None, design=None, df_ext=None):
    path_file = path_dir
    if path_dir is None:
        path_file = os.path.join(args.PATH)
    if design is not None:
        path_file = os.path.join(args.PATH, design)

    #TODO: add parameters for STAR and readcounts
    path_file = os.path.join(path_file, "STAR", "readcounts", file_name)

    df_diff = pd.read_csv(path_file, sep="\t")

    if "epcy" in method:
        df_diff = select_epcy(args, df_diff, exp_genes, method, df_ext)
    if method == "deseq2":
        df_diff = select_deseq(args, df_diff, exp_genes)
    if method == "edger":
        df_diff = select_edger(args, df_diff, exp_genes)
    if method == "limma":
        df_diff = select_limma(args, df_diff, exp_genes)

    return(df_diff)

def read_sample_pred_table(args, file_name, df_diff, dataset=None):
    path_file = args
    if type(args) is not str:
        path_file = os.path.join(args.PATH)
    if dataset is not None:
        path_file = os.path.join(args.PATH, dataset)

    path_file = os.path.join(path_file, file_name)

    df_pred = pd.read_table(path_file)
    df_pred = df_pred.merge(df_diff[["gene_id", "ID"]], left_on="IDS", right_on="gene_id")
    df_pred = df_pred.set_index("ID")
    df_pred = df_pred.drop("IDS", 1)
    df_pred = df_pred.drop("gene_id", 1)

    df_pred = df_pred.replace(2, 1)
    df_pred = df_pred.replace(3, -1)
    df_pred = df_pred.replace(4, -1)

    for i, row in df_diff.iterrows():
        df_pred.loc[row["ID"]] = df_pred.loc[row["ID"]] * 10**row["KERNEL_MCC"]

    return(df_pred)

def get_design(args, path=None, dataset=None, file_name="design.tsv"):
    design = path
    if path is None:
        design = os.path.join(args.PATH)
    if dataset is not None:
        design = os.path.join(args.PATH, dataset)

    design = os.path.join(design, file_name)
    design = pd.read_csv(design, sep="\t")
    design[args.SUBGROUP] = [1 if condition == args.QUERY else 0 for condition in design[args.SUBGROUP]]
    design = design.sort_values(by=[args.SUBGROUP, 'sample'], ascending=[False, True])

    return(design)

def get_exp(args, file_name, df_design=None):
    df_biotype = pd.read_csv(args.BF, sep="\t")

    if (not args.BIOTYPE is None):
        selected_biotype = args.BIOTYPE.split(",")
        df_biotype = df_biotype.loc[df_biotype["gene_biotype"].isin(selected_biotype)]

    df_exp = pd.read_csv(file_name, sep="\t")
    df_exp = df_exp.loc[df_exp["ID"].isin(df_biotype["ensembl_gene_id"])]

    df_exp = df_exp.set_index("ID")
    if df_design is not None:
        df_exp = df_exp[df_design["sample"]]

    if args.CPM:
        f_norm = 1e6 /  df_exp.sum()
        df_exp = df_exp * f_norm

    df_exp = np.log2(df_exp + 1)
    return(df_exp)

def get_silhouette_kmeans(df, metric):
    X = df.T.values
    n_clusters_ = 2
    kmeans = KMeans(n_clusters=n_clusters_)
    labels = kmeans.fit_predict(X)

    silh_score = silhouette_score(X, labels, metric=metric)

    return(silh_score, n_clusters_, labels)


def get_silhouette_af(df, df_design, metric, prefix_file_out, method, top):
    X = df.T.values
    af = AffinityPropagation().fit(X)
    cluster_centers_indices = af.cluster_centers_indices_
    n_clusters_ = len(cluster_centers_indices)
    labels = af.labels_

    silh_score = silhouette_score(X, labels, metric=metric)

    return(silh_score, n_clusters_, labels)

def get_silhouette_graph(df, labels, metric, prefix_file_out, method, top, silhouette_avg, n_clusters):
    X = df.T.values
    cluster_labels = labels

    # Create a subplot with 1 row and 2 columns
    fig, (ax1, ax2) = plt.subplots(1, 2)
    fig.set_size_inches(18, 7)

    # The 1st subplot is the silhouette plot
    # The silhouette coefficient can range from -1, 1 but in this example all
    # lie within [-0.1, 1]
    ax1.set_xlim([-0.1, 1])
    # The (n_clusters+1)*10 is for inserting blank space between silhouette
    # plots of individual clusters, to demarcate them clearly.
    ax1.set_ylim([0, len(X) + (n_clusters + 1) * 10])

    # Compute the silhouette scores for each sample
    sample_silhouette_values = silhouette_samples(X, cluster_labels, metric=metric)
    y_lower = 10
    for i in range(n_clusters):
        # Aggregate the silhouette scores for samples belonging to
        # cluster i, and sort them
        ith_cluster_silhouette_values = \
            sample_silhouette_values[cluster_labels == i]

        ith_cluster_silhouette_values.sort()

        size_cluster_i = ith_cluster_silhouette_values.shape[0]
        y_upper = y_lower + size_cluster_i

        color = cm.nipy_spectral(float(i) / n_clusters)
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
    colors = cm.nipy_spectral(cluster_labels.astype(float) / n_clusters)
    ax2.scatter(X[:, 0], X[:, 1], marker='.', s=30, lw=0, alpha=0.7,
                c=colors, edgecolor='k')

    # Labeling the clusters
    # centers = clusterer.cluster_centers_
    # Draw white circles at cluster centers
    #ax2.scatter(centers[:, 0], centers[:, 1], marker='o',
    #            c="white", alpha=1, s=200, edgecolor='k')

    #for i, c in enumerate(centers):
    #    ax2.scatter(c[0], c[1], marker='$%d$' % i, alpha=1,
    #                s=50, edgecolor='k')

    ax2.set_title("The visualization of the clustered data.")
    ax2.set_xlabel("Feature space for the 1st feature")
    ax2.set_ylabel("Feature space for the 2nd feature")

    plt.suptitle(("Silhouette analysis for KMeans clustering on sample data "
                  "with n_clusters = %d" % n_clusters),
                 fontsize=14, fontweight='bold')

    file_out = prefix_file_out + "_" + method + "_" + str(top) + "_silhouette.pdf"
    plt.savefig(file_out)
    plt.close()



def get_silhouette(df, df_design, metric):
    silh_score = silhouette_score(df.T, df_design["group"], metric=metric)

    return(silh_score)

def get_mcc(ct):
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

    (u_value, p_value) = mannwhitneyu(x[ids_query], x[ids_ref], alternative="two-sided")
    auc = u_value / (len(x[ids_query]) * len(x[ids_ref]))

    if auc < 0.5:
        auc = 1 - auc
    return(auc, p_value)

def get_sample_dist(df, df_design, metric):
    samples_dist = pdist(df.T, metric=metric)
    samples_dist = squareform(samples_dist)

    num_sample = len(df_design["group"])
    sample_query = df_design["sample"].loc[df_design["group"] == 1]

    close_perc = [0.0, 0.0]

    for group in [0, 1]:
        sample_selected = df_design["sample"].loc[df_design["group"] == group]
        num_sample_selected = len(sample_selected)
        ind_exp_selected = np.where(df.columns.isin(sample_selected))[0]

        sum_same_group = 0
        for index in ind_exp_selected:
            sample_dist = [samples_dist[index, j] for j in range(0,num_sample)]
            cutoff_dist = sorted(sample_dist)[num_sample_selected]
            closest_samples = np.where(sample_dist <= cutoff_dist)
            closest_samples = df.columns.values[closest_samples]

            num_same_group = np.sum(df_design["group"].loc[df_design["sample"].isin(closest_samples)] == group)
            num_same_group = np.min([num_same_group, num_sample_selected])
            sum_same_group += num_same_group - 1

        close_perc[group] = sum_same_group / (num_sample_selected * (num_sample_selected - 1))

    return(close_perc, samples_dist)

def display_cluster(args, df, df_design, prefix_file_out, method, metric, top, max_exp):
    col_pal = sns.color_palette("vlag", 200)

    samples_palette = dict(zip(df_design.group.unique(), [col_pal[199], col_pal[0]]))
    samples_colors = df_design.group.map(samples_palette)
    samples_colors = dict(zip(df_design["sample"], samples_colors.get_values()))
    samples_colors = pd.Series(df.columns, index=df.columns).map(samples_colors)

    #default seaborn size is (10, 10)
    figsize = (10, np.max([10, len(df.index)/4]))
    figsize = (10,10)
    display_label = False
    if args.SCALED:
        figsize=(len(df.columns)/4, len(df.index)/4)
        display_label = True
        prefix_file_out += "_scaled"

    sns_plot = sns.clustermap(
        df, figsize=figsize, linewidths=0, metric=metric,
        col_colors=samples_colors, #row_colors=mcc_colors, #cmap="viridis_r",
        xticklabels=display_label, yticklabels=False,
        vmin=0, vmax=max_exp#, cmap="Blues"
    )
    sns_plot.fig.suptitle(method)
    # plt.plot(sns_plot)
    file_out = prefix_file_out + "_" + method + "_" + str(top) + "_heatmap.pdf"
    sns_plot.savefig(file_out)
    plt.close()

def get_dist_by_group(df, df_design, samples_dist, method, top, metric):
    num_sample = len(df_design["group"])

    dict_dist = defaultdict(list) #pd.DataFrame(columns=['group', 'dist', 'metric'])
    dist_diff = []
    dist_diff_sample = []
    for i in range(0,num_sample):
        if i != num_sample - 1:
            group_i = df_design["group"].loc[df_design["sample"] == df.columns.values[i] ].values[0]
            for j in range(i + 1,num_sample):
                group_j = df_design["group"].loc[df_design["sample"] == df.columns.values[j] ].values[0]
                if group_i == group_j:
                    if group_i == 1:
                        dict_dist['group'].append("Query")
                    else:
                        dict_dist['group'].append("Ref")
                else:
                    dict_dist['group'].append("diff")
                    dist_diff.append(samples_dist[i, j])
                    dist_diff_sample.append([df.columns.values[i], df.columns.values[j]])
                dict_dist['dist'].append(samples_dist[i, j])
                dict_dist['metric'].append(metric)
                dict_dist['method'].append(method)
                dict_dist['top'].append(top)
    return(dict_dist, dist_diff, dist_diff_sample)


def run_kmeans(XTEST, YTEST, design, method, dataset, top):
    silh_score, n_clusters_, labels = get_silhouette_kmeans(XTEST, "euclidean")
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

def run_rand_forest(XTRAIN, YTRAIN, XTEST, YTEST, design, method, dataset, top):
    rf = RandomForestClassifier(n_estimators=100, class_weight="balanced_subsample")
    rf.fit(XTRAIN.T, YTRAIN)

    proba_test = rf.predict_proba(XTEST.T)
    proba_train = rf.predict_proba(XTRAIN.T)

    id_dataset = re.findall(r'\d+', dataset)[0]
    id_dataset = int(id_dataset)

    dict_rf = defaultdict(list)
    dict_rf['method'].append(method)
    dict_rf['top'].append(top)
    dict_rf['label'].append(YTEST.values)
    dict_rf['proba'].append(proba_test)
    dict_rf['label_train'].append(YTRAIN.values)
    dict_rf['proba_train'].append(proba_train)
    dict_rf['dataset'].append(id_dataset)
    dict_rf['design'].append(design)
    dict_rf['learn'].append("R-Forest")

    return(dict_rf)


def run_LR(XTRAIN, YTRAIN, XTEST, YTEST, design, method, dataset, top):
    ### large C ==> no regularization
    LR = linear_model.LogisticRegression(solver="liblinear", C = 1e10, penalty = 'l2', max_iter=1000)
    LR.fit(XTRAIN.T, YTRAIN)


    #tpm_test = np.reshape(df_exp_test, (-1, 1))
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

    return(dict_LR)

def get_exp_cluster(args, df_exp, df_design, df_diff, prefix_file_out, method, design, top, max_exp, metric='euclidean'):
    df_exp = df_exp.loc[df_exp.index.isin(df_diff["ID"][:top])]
    return(get_cluster(args, df_exp, df_design, prefix_file_out, method, design, top, max_exp, metric=metric))

def get_cluster(args, df, df_design, prefix_file_out, method, design, top, max_exp, metric='euclidean'):
    display_cluster(args, df, df_design, prefix_file_out, method, metric, top, max_exp)
    (silh_score, num_cluster, labels) = get_silhouette_kmeans(df, metric)
    get_silhouette_graph(df, labels, metric, prefix_file_out, method, top, silh_score, num_cluster)

    #pred = [labels[i] == df_design["group"][i] for i in range(len(df_design["group"]))]
    #pred = np.sum(pred) / len(df_design["group"])
    #if pred < 0.5:
    #    pred = 1 - pred

    return(silh_score)
