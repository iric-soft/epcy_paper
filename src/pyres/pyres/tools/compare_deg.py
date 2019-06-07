import os

from collections import defaultdict

import seaborn as sns
import pandas as pd
import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'DejaVu Sans'
mpl.rcParams['pdf.fonttype'] = 42

from .. utils import other as uo

def main_compare_deg(args, argparser):
    path_file_out = os.path.join(args.PATH, "res", "deg")
    if not os.path.exists(path_file_out):
        os.makedirs(path_file_out)

    if args.GENE:
        prefix_file_out = "genes"

        file_pred = "genes_prediction.xls"
        file_deseq2 = "deseq2_genes_h5.xls"
        file_deseq2_lrt = "deseq2_genes_h5.xls"
        file_edger_lrt = "edger_lrt_genes_h5.xls"
        file_edger_qlf = "edger_qlf_genes_h5.xls"
        file_limma = "limma_voom_genes_h5.xls"

        col_ids_name = "gene_id"
        file_exp = "/u/eaudemard/project/kt_paper/data/all_exp_404.csv"
        #file_exp = "/u/eaudemard/project/kt_paper/data/all_exp_tcga_brca.csv"
    else:
        prefix_file_out = "trans"

        file_pred = "trans_prediction.xls"
        file_deseq2 = "deseq2_trans_h5.xls"
        file_deseq2_lrt = "deseq2_trans_h5.xls"
        file_edger_lrt = "edger_lrt_trans_h5.xls"
        file_edger_qlf = "edger_qlf_trans_h5.xls"
        file_limma = "limma_voom_trans_h5.xls"

        col_ids_name = "trans_id"
        file_exp = "/u/eaudemard/project/kt_paper/data/"

    prefix_file_out += "_" + args.PREFIX
    prefix_file_out += "_lfc" + str(args.LOG_FC)
    prefix_file_out += "_MCC" + str(args.MCC)

    subgroups = args.SUBGROUP.split(",")

    dict_perf = None
    dict_silh = None
    dict_diff = defaultdict(list)
    for subgroup in subgroups:
        print(subgroup)
        path_subgroup_file_out = os.path.join(args.PATH, subgroup, "res", "deg")

        if not os.path.exists(path_subgroup_file_out):
            os.makedirs(path_subgroup_file_out)

        path_out = os.path.join(path_subgroup_file_out, prefix_file_out)

        dir_subgroup = os.path.join(args.PATH, subgroup)
        df_design = uo.get_design(args, path=dir_subgroup)
        df_exp = uo.get_exp(args, file_exp)
        max_exp = np.amax(df_exp.values)

        #print(df_exp)
        #print(df_design["sample"])
        df_exp = df_exp[df_design["sample"]]
        list_genes = df_exp.index.tolist()

        dict_dist = None
        for top in [50]:
            print(top)
            dict_diff["kt"] = uo.read_diff_table(args, file_pred, "kt", list_genes, path_dir=dir_subgroup)
            dict_diff["deseq2"] = uo.read_diff_table(args, file_deseq2, "deseq2", list_genes, path_dir=dir_subgroup)
            #dict_diff["deseq2_lrf"] = uo.read_diff_table(args, file_deseq2_lrt, "deseq2", list_genes, path_dir=dir_subgroup)
            dict_diff["edger"] = uo.read_diff_table(args, file_edger_qlf, "edger", list_genes, path_dir=dir_subgroup)
            #dict_diff["edger_qlf"] = uo.read_diff_table(args, file_edger_qlf, "edger", list_genes, path_dir=dir_subgroup)
            dict_diff["limma"] = uo.read_diff_table(args, file_limma, "limma", list_genes, path_dir=dir_subgroup)

            for method in ["edger", "kt", "deseq2", "limma"]: #deseq2_lrf, edger_qlf
                print(method)
                (dict_dist_tmp, dict_perf_tmp, dict_silh_tmp) = uo.get_exp_cluster(args, df_exp, df_design, dict_diff[method], path_out, method, subgroup, top, max_exp)

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
        dist_density = sns.FacetGrid(df_dist, col="method", hue="group", row="metric") #, col="overlap", hue="group"

        kind = "violin"
        if subgroup == "28_inv16_vs28":
            kind = "swarm"

        dist_density = sns.catplot(data=df_dist, x="group", y="dist", col="method", row="top", kind=kind)#, facet_kws=dict(subplot_kws=dict(alpha=.1)))
        plt.yscale('log')
        plt.legend();

        file_out = path_out + "_dist_density.pdf"
        dist_density.savefig(file_out)
        plt.close()

    path_out = os.path.join(path_file_out, prefix_file_out)

    df_perf = pd.DataFrame(dict_perf)
    #print(df_dist)
    dist_density = sns.catplot(data=df_perf, x="top", y="perc", col="subgroup", hue="method", row="group", kind="point")#, facet_kws=dict(subplot_kws=dict(alpha=.1)))
    file_out = path_out + "_perc_query.pdf"
    dist_density.savefig(file_out)
    plt.close()

    df_silh = pd.DataFrame(dict_silh)

    dist_density = sns.catplot(data=df_silh, x="top", y="min_dist", col="subgroup", hue="method", row="metric", kind="point")#, facet_kws=dict(subplot_kws=dict(alpha=.1)))
    file_out = path_out + "_min_dist.pdf"
    dist_density.savefig(file_out)
    plt.close()

    dist_density = sns.catplot(data=df_silh, x="top", y="silhouette_score", col="subgroup", hue="method", row="metric", kind="point")#, facet_kws=dict(subplot_kws=dict(alpha=.1)))
    file_out = path_out + "_silhouette.pdf"
    dist_density.savefig(file_out)
    plt.close()

    dist_density = sns.catplot(data=df_silh, x="top", y="pred", col="subgroup", hue="method", row="metric", kind="point")#, facet_kws=dict(subplot_kws=dict(alpha=.1)))
    file_out = path_out + "_pred.pdf"
    dist_density.savefig(file_out)
    plt.close()

    for subgroup in subgroups:
        df_silh_tmp = df_silh.loc[df_silh["subgroup"] == subgroup]
        dist_density = sns.pointplot(data=df_silh_tmp, x="top", y="silhouette_score", hue="method")#, facet_kws=dict(subplot_kws=dict(alpha=.1)))
        file_out = path_out + "_" + subgroup + "_silhouette.pdf"
        plt.savefig(file_out)
        plt.close()
