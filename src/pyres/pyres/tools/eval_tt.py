import datetime
# get time and date
start_time = datetime.datetime.now()
# make rest of imports
import os

from collections import defaultdict
from collections import Counter

import pdb

import pandas as pd
import numpy as np

import matplotlib as mpl
# set no display
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

from .. utils import other as uo
from sklearn.utils import shuffle as skl_shuffle

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'DejaVu Sans'
mpl.rcParams['pdf.fonttype'] = 42

def main_eval_tt(args, argparser):

    file_dict = {
        'epcy' : "prediction_capability.xls",
        'deseq2' : "deseq2_genes.xls",
        'edger' : "edger_genes.xls",
        'limma' : "limma_voom_genes.xls",
    }

    # create dict with search params
    search_params_dict = {
        'tops' : args.TOP_VALUES,
        'methods' : args.METHODS,
        'pvalues' : args.PVALUES
    }


    df_exp = uo.get_exp(args, args.MATRIX)
    list_genes = df_exp.index.tolist()

    design = os.path.join(args.PATH, "leucegene_train", "all", "design.tsv")
    all_train_labels = pd.read_csv(design, sep="\t")
    all_train_labels["subgroup_str"] = all_train_labels["subgroup"]
    all_train_labels["subgroup"].replace(
        ["inv16", "MLL", "CK", "EVI1", "Inter", "Mono5", "NK", "Tri8", "NUP98NSD1", "t8_21", "t15_17"],
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
        inplace=True
    )
    all_train_labels = all_train_labels.sort_values(["subgroup"], ascending=True)

    design = os.path.join(args.PATH, "leucegene_test", "all", "design.tsv")
    test_labels = pd.read_csv(design, sep="\t")
    test_labels["subgroup_str"] = test_labels["subgroup"]
    test_labels["subgroup"].replace(
        ["inv16", "MLL", "CK", "EVI1", "Inter", "Mono5", "NK", "Tri8", "NUP98NSD1", "t8_21", "t15_17"],
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
        inplace=True
    )
    test_labels = test_labels.sort_values(["subgroup"], ascending=True)

    dict_all = defaultdict(list)
    dict_subg = defaultdict(list)
    for p_train in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]:
        train_labels = pd.DataFrame()
        count_subg = Counter(all_train_labels["subgroup_str"])
        for i in count_subg:
            num_keep = int(count_subg[i]*p_train)
            tmp = all_train_labels[all_train_labels["subgroup_str"] == i][:num_keep]
            train_labels = pd.concat([train_labels, tmp])

        test = uo.get_exp(args, args.MATRIX, test_labels)
        train = uo.get_exp(args, args.MATRIX, train_labels)

        dir_design = os.path.join(args.PATH, "leucegene_train")
        designs = os.listdir(dir_design)
        designs.remove("all")

        for pvalue in search_params_dict['pvalues']:
            for top in search_params_dict['tops']:
                for method in search_params_dict['methods']:
                    print('TRAINING/TESTING p_train: {}, pvalue: {}, top: {}, method: {}'.format(p_train, pvalue, top, method))
                    features_ids_keep = []
                    for design in designs:
                        path_design = os.path.join(dir_design, design)
                        df_design_train = uo.get_design(args, path_design)
                        df_diff = uo.read_diff_table(args, file_dict[method], method, path_design, pvalue, list_genes)
                        features_ids_keep += df_diff["ID"][:top].tolist()

                    dict_subg_tmp, dict_all_tmp = uo.run_train_label(args, method, train, train_labels, test, test_labels, top, pvalue, p_train, features_ids_keep=features_ids_keep)

                    for key in dict_subg_tmp.keys():
                       dict_subg[key] = dict_subg[key] + dict_subg_tmp[key]

                    for key in dict_all_tmp.keys():
                       dict_all[key] = dict_all[key] + dict_all_tmp[key]


    df_all = pd.DataFrame(dict_all)
    df_subg = pd.DataFrame(dict_subg)
    fig_dir = os.path.join(args.OUTDIR)

    plt_fig = sns.catplot(data=df_subg, x="p_train", y="mcc", hue="method", row="subg", col = "top", kind="point", facet_kws=dict(subplot_kws=dict(ylim=[-0.05,1.05], xlim=[-0.05,1.05])))
    fig_out = os.path.join(fig_dir, "pred_by_subg_mcc.pdf")
    plt_fig.savefig(fig_out)
    plt.close()

    plt_fig = sns.catplot(data=df_all, x="p_train", y="mcc", hue="method", row = "pvalue", col = "top", kind="point", facet_kws=dict(subplot_kws=dict(ylim=[-0.05,1.05], xlim=[-0.05,1.05])))
    fig_out = os.path.join(fig_dir, "pred_mcc.pdf")
    plt_fig.savefig(fig_out)
    plt.close()

    plt_fig = sns.catplot(data=df_all, x="p_train", y="full_ARS", hue="method", row = "pvalue", col = "top", kind="point", facet_kws=dict(subplot_kws=dict(ylim=[-0.05,1.05], xlim=[-0.05,1.05])))
    fig_out = os.path.join(fig_dir, "pred_ARS.pdf")
    plt_fig.savefig(fig_out)
    plt.close()

    plt_fig = sns.catplot(data=df_all, x="p_train", y="full_MIS", hue="method", row = "pvalue", col = "top", kind="point", facet_kws=dict(subplot_kws=dict(ylim=[-0.05,1.05], xlim=[-0.05,1.05])))
    fig_out =  os.path.join(fig_dir, "pred_MIS.pdf")
    plt_fig.savefig(fig_out)
    plt.close()

    plt_fig = sns.catplot(data=df_all, x="p_train", y="num_feature", hue="method", row = "pvalue", col = "top", kind="point")
    fig_out = os.path.join(fig_dir, "pred_num_feature.pdf")
    plt_fig.savefig(fig_out)
    plt.close()
