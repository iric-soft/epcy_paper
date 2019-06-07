import os

from collections import defaultdict

import seaborn as sns
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'DejaVu Sans'
mpl.rcParams['pdf.fonttype'] = 42

from .. utils import other as uo

def main_compare_mcc(args, argparser):

    path_file_out = os.path.join(args.PATH, "res", "mcc")
    if not os.path.exists(path_file_out):
        os.makedirs(path_file_out)

    if args.GENE:
        prefix_file_out = "genes"

        file_pred = "genes_prediction.xls"
        file_pred_sample = "samples_genes_prediction.xls"
        file_pred_bs = "bs_100/kal/genes_prediction.xls"

        col_ids_name = "gene_id"
        file_exp = "/u/eaudemard/project/kt_paper/data/all_exp_404.csv"
        #file_exp = "/u/eaudemard/project/kt_paper/data/all_exp_tcga_brca.csv"
    else:
        prefix_file_out = "trans"

        file_pred = "trans_prediction.xls"
        file_pred_sample = "samples_trans_prediction.xls"
        file_pred_bs = "bs_100/kal/trans_prediction.xls"

        col_ids_name = "trans_id"
        file_exp = "/u/eaudemard/project/kt_paper/data/"

    prefix_file_out += "_" + args.PREFIX
    prefix_file_out += "_lfc" + str(args.LOG_FC)
    prefix_file_out += "_MCC" + str(args.MCC)

    datasets = args.DATASET.split(",")

    dict_perf = None
    dict_silh = None
    dict_diff = defaultdict(list)
    for dataset in datasets:
        print(dataset)
        path_dataset_file_out = os.path.join(args.PATH, dataset, "res", "mcc")

        if not os.path.exists(path_dataset_file_out):
            os.makedirs(path_dataset_file_out)

        path_out = os.path.join(path_dataset_file_out, prefix_file_out)

        df_design = uo.get_design(args, dataset)
        df_exp = uo.get_exp(args, file_exp, df_design)
        list_genes = df_exp.index.tolist()

        dict_dist = None
        for top in [2, 5, 10, 20, 50, 100, 200]:
            print(top)

            for method in ["kt", "kt_mcc", "kt_binom", "kt_auc"]:
                dict_diff[method] = uo.read_diff_table(args, file_pred, method, list_genes, dataset)
            df_pred_sample = uo.read_sample_pred_table(args, file_pred_sample, dict_diff["kt"], dataset)
            dict_diff["kt_bs"] = uo.read_diff_table(args, file_pred_bs, "kt_bs", list_genes, dataset)

            for method in ["kt", "kt_mcc", "kt_binom", "kt_auc", "kt_bs"]:
                if method == "kt_mcc":
                    (dict_dist_tmp, dict_perf_tmp, dict_silh_tmp) = uo.get_cluster(args, df_pred_sample, df_design, path_out, method, dataset, top)
                else:
                    (dict_dist_tmp, dict_perf_tmp, dict_silh_tmp) = uo.get_exp_cluster(args, df_exp, df_design, dict_diff[method], path_out, method, dataset, top)

                if dict_dist is None:
                    dict_dist = dict_dist_tmp
                else:
                    for key in dict_dist.keys():
                        dict_dist[key] = dict_dist[key] + dict_dist_tmp[key]

                if dict_perf is None:
                    dict_perf = dict_perf_tmp
                else:
                    for key in dict_perf.keys():
                        dict_perf[key] = dict_perf[key] + dict_perf_tmp[key]

                if dict_silh is None:
                    dict_silh = dict_silh_tmp
                else:
                    for key in dict_silh.keys():
                        dict_silh[key] = dict_silh[key] + dict_silh_tmp[key]

        df_dist = pd.DataFrame(dict_dist)
        #print(df_dist)
        bw = 2
        kind = "violin"
        if dataset == "28_inv16_vs28":
            kind = "swarm"

        dist_density = sns.catplot(data=df_dist, x="group", y="dist", col="method", row="top", kind=kind)#, facet_kws=dict(subplot_kws=dict(alpha=.1)))
        file_out = path_out + "_dist_density.pdf"
        dist_density.savefig(file_out)
        plt.close()

    path_out = os.path.join(path_file_out, prefix_file_out)

    df_perf = pd.DataFrame(dict_perf)
    #print(df_dist)
    bw = 2
    dist_density = sns.catplot(data=df_perf, x="top", y="perc", col="dataset", hue="method", row="group", kind="point")#, facet_kws=dict(subplot_kws=dict(alpha=.1)))
    file_out = path_out + "_perc_query.pdf"
    dist_density.savefig(file_out)
    plt.close()

    df_silh = pd.DataFrame(dict_silh)

#    bw = 2
#    dist_density = sns.catplot(data=df_silh, x="top", y="min_dist", col="dataset", hue="method", row="metric", kind="point")#, facet_kws=dict(subplot_kws=dict(alpha=.1)))
#    file_out = path_out + "_min_dist.pdf"
#    dist_density.savefig(file_out)
#    plt.close()

    bw = 2
    dist_density = sns.catplot(data=df_silh, x="top", y="silhouette_score", col="dataset", hue="method", row="metric", kind="point")#, facet_kws=dict(subplot_kws=dict(alpha=.1)))
    file_out = path_out + "_silhouette.pdf"
    dist_density.savefig(file_out)
    plt.close()

    bw = 2
    dist_density = sns.catplot(data=df_silh, x="top", y="pred", col="dataset", hue="method", row="metric", kind="point")#, facet_kws=dict(subplot_kws=dict(alpha=.1)))
    file_out = path_out + "_pred.pdf"
    dist_density.savefig(file_out)
    plt.close()

    for dataset in datasets:
        bw = 2
        df_silh_tmp = df_silh.loc[df_silh["dataset"] == dataset]
        dist_density = sns.pointplot(data=df_silh_tmp, x="top", y="silhouette_score", hue="method")#, facet_kws=dict(subplot_kws=dict(alpha=.1)))
        file_out = path_out + "_" + dataset + "_silhouette.pdf"
        plt.savefig(file_out)
        plt.close()
