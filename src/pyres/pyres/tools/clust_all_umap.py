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

def main_clust_all_umap(args, argparser):
    path_file_out = os.path.join(args.PATH, "res", "deg")
    if not os.path.exists(path_file_out):
        os.makedirs(path_file_out)

    file_dict = {
        'epcy' : "prediction_capability.xls",
        'deseq2' : "deseq2_genes.xls",
        'edger' : "edger_genes.xls",
        'limma' : "limma_voom_genes.xls",
    }

    # create dict with search params
    search_params_dict = {
        'methods' : args.METHODS,
        'designs' : args.DESIGN
    }

    dir_design = os.path.join(args.PATH, "all")
    df_design = uo.get_design(args, path=dir_design, replace=False)

    df_exp = uo.get_exp(args, args.MATRIX)
    df_exp = df_exp[df_design["sample"]]
    list_genes = df_exp.index.tolist()

    path_out = os.path.join(args.OUTDIR, "clust_all", "umap")
    if not os.path.exists(path_out):
        os.makedirs(path_out)

    dict_umap = uo.cluster_umap(args, df_exp, df_design, path_out, "no_method", "all_design", "all_exp", metric='euclidean')


    dict_mcc = defaultdict(list)
    for method in search_params_dict['methods']:
        gene_ids_keep = []
        for design in search_params_dict['designs']:
            print("##{}, {}".format(method.upper(), design))
            dir_design = os.path.join(args.PATH, design)
            df_diff = uo.read_diff_table(args, file_dict[method], method, dir_design, args.PVALUE, list_genes)
            gene_ids_keep += df_diff["ID"][:10].tolist()

            if method == "epcy":
                for mcc in df_diff['KERNEL_MCC'][:10].tolist():
                    dict_mcc['design'].append(design)
                    dict_mcc['MCC'].append(mcc)

        gene_ids_keep = list(set(gene_ids_keep))
        print(len(gene_ids_keep))
        df_exp_method = df_exp.loc[gene_ids_keep]
        df_exp_method.sort_index(inplace=True)

        print(str(df_exp_method.shape[0]) + " " + str(df_exp_method.shape[1]))
        dict_umap_tmp = uo.cluster_umap(args, df_exp_method, df_design, path_out, method, "all_design", "selected_exp", metric='euclidean')
        for key in dict_umap_tmp.keys():
            dict_umap[key] = dict_umap[key] + dict_umap_tmp[key]

    # create subfolder for log and res
    fig_dir = os.path.join(args.OUTDIR, "clust_all", "umap")
    if not os.path.exists(fig_dir):
        os.makedirs(fig_dir)

    # convert dict to df for analysis
    df_res = pd.DataFrame(dict_umap)
    csv_out = os.path.join(fig_dir, "umap_table.csv")
    df_res.to_csv(csv_out, index=False, sep="\t")

    # plot results
    plt_fig = sns.catplot(data=df_res, x="top", y="full_ARS", hue="method", col="design", kind="point", facet_kws=dict(subplot_kws=dict(ylim=[-0.05,1.05])))
    fig_out =  os.path.join(fig_dir, "full_ARS." + args.EXT)
    plt_fig.savefig(fig_out)
    plt.close()

    plt_fig = sns.catplot(data=df_res, x="top", y="full_MIS", hue="method", col="design", kind="point", facet_kws=dict(subplot_kws=dict(ylim=[-0.05,1.05])))
    fig_out =  os.path.join(fig_dir, "full_MIS." + args.EXT)
    plt_fig.savefig(fig_out)
    plt.close()

    plt_fig = sns.catplot(data=df_res, x="top", y="ARS", hue="method", col="design", kind="point", facet_kws=dict(subplot_kws=dict(ylim=[-0.05,1.05])))
    fig_out =  os.path.join(fig_dir, "ARS." + args.EXT)
    plt_fig.savefig(fig_out)
    plt.close()

    plt_fig = sns.catplot(data=df_res, x="top", y="MIS", hue="method", col="design", kind="point", facet_kws=dict(subplot_kws=dict(ylim=[-0.05,1.05])))
    fig_out =  os.path.join(fig_dir, "MIS." + args.EXT)
    plt_fig.savefig(fig_out)
    plt.close()

    plt_fig = sns.catplot(data=df_res, x="top", y="p_clustered", hue="method", col="design", kind="point", facet_kws=dict(subplot_kws=dict(ylim=[-0.05,1.05])))
    fig_out =  os.path.join(fig_dir, "p_clustered." + args.EXT)
    plt_fig.savefig(fig_out)
    plt.close()

    df_res = pd.DataFrame(dict_mcc)
    plt_fig = sns.boxplot(data=df_res, x="design", y="MCC")
    plt_fig.set_xticklabels(plt_fig.get_xticklabels(),rotation=30)
    plt_fig.set(ylim=(-0.05, 1.05))
    fig_out =  os.path.join(fig_dir, "mcc_boxplot." + args.EXT)
    plt_fig = plt_fig.get_figure()
    plt_fig.savefig(fig_out)
    plt.close()
