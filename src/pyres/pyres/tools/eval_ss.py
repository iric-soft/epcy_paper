

# get time and date
# make rest of imports
import os
import sys

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

import datetime

start_time = datetime.datetime.now()

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'DejaVu Sans'
mpl.rcParams['pdf.fonttype'] = 42


def main_eval_ss(args, argparser):

    file_dict = {
        'epcy': "predictive_capability.xls",
        'epcy_bagging': "predictive_capability.xls",
        'deseq2': "deseq2_genes.xls",
        'edger': "edger_genes.xls",
        'voom': "limma_voom_genes.xls",
        'trend': "limma_trend_genes.xls",
        'mast': "mast_genes.xls",
    }

    # create dict with search params
    search_params_dict = {
        'reps': args.REPS,
        'p_ss': args.P_SS,
        'methods': args.METHODS
    }
    design = args.DESIGN

    cohort_order = []
    dict_cohort = defaultdict(list)
    for p_ss in search_params_dict['p_ss']:
        if p_ss == 0.0:
            p_ss = 0
        rep = 1
        path_design = os.path.join(args.PATH, design, str(p_ss), str(rep))
        df_design = uo.get_design(args, "Query", path_design)
        num_ref = df_design[df_design[args.SUBGROUP] == 0].count()[0]
        num_query = df_design[df_design[args.SUBGROUP] == 1].count()[0]
        cohort = str(num_query) + "_vs_" + str(num_ref)
        cohort_order.append(cohort)

        dict_cohort_tmp = defaultdict(list)
        dict_cohort_tmp['Query'].append(num_query)
        dict_cohort_tmp['Ref'].append(num_ref)

        for key in dict_cohort_tmp.keys():
            dict_cohort[key] = dict_cohort[key] + dict_cohort_tmp[key]

    df_cohort = pd.DataFrame(dict_cohort)

    log_c = sys.float_info.min
    dict_diff = defaultdict(list)
    for p_ss in search_params_dict['p_ss']:
        for rep in search_params_dict['reps']:
            for method in search_params_dict['methods']:
                if p_ss == 0.0:
                    p_ss = 0
                print('SUBSAMPLING: design: {}, p_ss: {}, rep: {}, method: {}'.format(design, p_ss, rep, method))
                path_design = os.path.join(args.PATH, design, str(p_ss), str(rep))
                df_design = uo.get_design(args, "Query", path_design)
                num_ref = df_design[df_design[args.SUBGROUP] == 0].count()[0]
                num_query = df_design[df_design[args.SUBGROUP] == 1].count()[0]
                cohort = str(num_query) + "_vs_" + str(num_ref)
                df_diff = uo.read_diff_table(
                    args, file_dict[method], method,
                    path_design, 1
                )

                dict_diff_tmp = defaultdict(list)
                dict_diff_tmp['ID'] = df_diff["ID"].tolist()
                dict_diff_tmp['rep'] = [rep] * df_diff["ID"].shape[0]
                dict_diff_tmp['method'] = [method] * df_diff["ID"].shape[0]
                dict_diff_tmp['cohort'] = [cohort] * df_diff["ID"].shape[0]

                if method == "epcy" or method == "epcy_bagging":
                    dict_diff_tmp['value'] = df_diff["KERNEL_MCC"].tolist()
                elif method == "deseq2":
                    dict_diff_tmp['value'] = (-np.log10(df_diff["padj"] + log_c)).tolist()
                elif method == "edger":
                    dict_diff_tmp['value'] = (-np.log10(df_diff["FDR"] + log_c)).tolist()
                elif method == "voom":
                    dict_diff_tmp['value'] = (-np.log10(df_diff["adj.P.Val"] + log_c)).tolist()
                elif method == "trend":
                    dict_diff_tmp['value'] = (-np.log10(df_diff["adj.P.Val"] + log_c)).tolist()
                elif method == "mast":
                    dict_diff_tmp['value'] = (-np.log10(df_diff["pval"] + log_c)).tolist()

                dict_diff_tmp['top'] = [x for x in range(1, df_diff.shape[0]+1, 1)]

                for key in dict_diff_tmp.keys():
                    dict_diff[key] = dict_diff[key] + dict_diff_tmp[key]

    df_diff = pd.DataFrame(dict_diff)

    fig_dir = os.path.join(args.OUTDIR, "eval_ss", args.DESIGN)
    if not os.path.exists(fig_dir):
        os.makedirs(fig_dir)

    bar_plot = df_cohort.plot(kind="bar", stacked=True)
    fig_out = os.path.join(fig_dir, "sub_sampling.pdf")
    bar_plot.figure.savefig(fig_out)
    plt.close()

    for id in args.IDS:
        df_diff_id = df_diff[df_diff["ID"] == id]
        df_diff_id_epcy = df_diff_id[(df_diff_id["method"] == "epcy") | (df_diff_id["method"] == "epcy_bagging")]
        df_diff_id_deg = df_diff_id[(df_diff_id["method"] != "epcy") & (df_diff_id["method"] != "epcy_bagging")]

        sns_plot = sns.pointplot(
            data=df_diff_id_epcy,
            x="cohort", y="value", hue="method", dodge=0.3,
            linestyles=[':'],
            markers=['x']
        )
        sns_plot.set_title("EPCY " + id)
        sns_plot.set(ylim=(-1.05, 1.05))
        sns_plot.set(ylabel="MCC")
        sns_plot.set_xticklabels(sns_plot.get_xticklabels(), rotation=90)
        fig_out = os.path.join(fig_dir, id + "_epcy.pdf")
        sns_plot.figure.savefig(fig_out)
        plt.close()

        sns_plot = sns.pointplot(
            data=df_diff_id_deg,
            x="cohort", y="value", hue="method", dodge=0.3,
            linestyles=['-', '--', '-.'],
            markers=['o', 'v', 's']
        )

        sns_plot.set_title("DEG " + id)
        sns_plot.axes.axhline(-np.log10(0.05), ls='--', color="r")
        sns_plot.set(ylabel="-log10(padj)")
        sns_plot.set_xticklabels(sns_plot.get_xticklabels(), rotation=90)
        fig_out = os.path.join(fig_dir, id + "_DEG.pdf")
        sns_plot.figure.savefig(fig_out)
        plt.close()

    df_diff_epcy = df_diff[(df_diff["method"] == "epcy") & (df_diff["value"] >= 0.2)]
    df_diff_epcy = df_diff_epcy.groupby(["method", "cohort", "rep"]).size().reset_index(name='counts')

    df_diff_DEG = df_diff[(df_diff["method"] != "epcy") & (df_diff["value"] >= -np.log10(0.05 + log_c))]
    df_diff_DEG = df_diff_DEG.groupby(["method", "cohort", "rep"]).size().reset_index(name='counts')

    df_merge_count = pd.concat([df_diff_epcy, df_diff_DEG], sort=False)

    sns_plot = sns.pointplot(
        x="cohort", y="counts", hue="method",
        hue_order=args.METHODS, order=cohort_order, dodge=0.3,
        linestyles=['-', '--', '-.', ':'], data=df_merge_count,
        markers=['o', 'v', 's', 'x']
    )
    sns_plot.set_title("Number of selected genes")
    sns_plot.set_xticklabels(sns_plot.get_xticklabels(), rotation=90)
    fig_out = os.path.join(fig_dir, "num_gene_selected_cat.pdf")
    sns_plot.figure.savefig(fig_out)
    plt.close()

    df_diff_top = df_diff[df_diff["top"] <= 20]
    df_diff_top_epcy = df_diff_top[(df_diff_top["method"] == "epcy") | (df_diff_top["method"] == "epcy_bagging")]
    df_diff_top_deg = df_diff_top[(df_diff_top["method"] != "epcy") & (df_diff_top["method"] != "epcy_bagging")]

    sns_plot = sns.pointplot(
        data=df_diff_top_epcy,
        x="cohort", y="value", hue="method", dodge=0.3,
        linestyles=[':'], order=cohort_order,
        markers=['x']
    )
    sns_plot.set_title("EPCY top 20")
    sns_plot.set(ylabel="MCC")
    sns_plot.set_xticklabels(sns_plot.get_xticklabels(), rotation=90)
    fig_out = os.path.join(fig_dir, "Top150_values_epcy.pdf")
    sns_plot.figure.savefig(fig_out)
    plt.close()

    sns_plot = sns.pointplot(
        data=df_diff_top_deg,
        x="cohort", y="value", hue="method", dodge=0.3,
        linestyles=['-', '--', '-.'], order=cohort_order,
        markers=['o', 'v', 's']
    )

    sns_plot.set_title("DEG top 20")
    sns_plot.axes.axhline(-np.log10(0.05), ls='--', color="r")
    sns_plot.set(ylabel="-log10(padj)")
    sns_plot.set_xticklabels(sns_plot.get_xticklabels(), rotation=90)
    fig_out = os.path.join(fig_dir, "Top150_values_DEG.pdf")
    sns_plot.figure.savefig(fig_out)
    plt.close()

    plt_fig = sns.swarmplot(
        x="cohort", y="value",
        data=df_diff_top_epcy, hue="method", order=cohort_order
    )
    plt_fig.set(ylabel="MCC")
    plt_fig.set_title("EPCY top 20")
    plt_fig.set_xticklabels(plt_fig.get_xticklabels(), rotation=90)
    fig_out = os.path.join(fig_dir, "Top150_values_epcy_dot.pdf")
    plt_fig.figure.savefig(fig_out)
    plt.close()

    df_diff_top_deseq = df_diff_top_deg[df_diff_top_deg["method"] == "deseq2"]
    plt_fig = sns.swarmplot(
        x="cohort", y="value",
        data=df_diff_top_deseq, hue="method", order=cohort_order
    )
    plt_fig.set_title("DESeq2 top 20")
    plt_fig.set(ylabel="-log10(padj)")
    plt_fig.set_xticklabels(plt_fig.get_xticklabels(), rotation=90)
    fig_out = os.path.join(fig_dir, "Top150_values_deseq2_dot.pdf")
    plt_fig.figure.savefig(fig_out)
    plt.close()

    df_diff_top_edger = df_diff_top_deg[df_diff_top_deg["method"] == "edger"]
    plt_fig = sns.swarmplot(
        x="cohort", y="value",
        data=df_diff_top_edger, hue="method", order=cohort_order
    )
    plt_fig.set_title("EdgeR top 10")
    plt_fig.set(ylabel="-log10(padj)")
    plt_fig.set_xticklabels(plt_fig.get_xticklabels(), rotation=90)
    fig_out = os.path.join(fig_dir, "Top150_values_edger_dot.pdf")
    plt_fig.figure.savefig(fig_out)
    plt.close()

    df_diff_top_voom = df_diff_top_deg[df_diff_top_deg["method"] == "voom"]
    plt_fig = sns.swarmplot(
        x="cohort", y="value",
        data=df_diff_top_voom, hue="method", order=cohort_order
    )
    plt_fig.set_title("Voom top 20")
    plt_fig.set(ylabel="-log10(padj)")
    plt_fig.set_xticklabels(plt_fig.get_xticklabels(), rotation=90)
    fig_out = os.path.join(fig_dir, "Top150_values_voom_dot.pdf")
    plt_fig.figure.savefig(fig_out)
    plt.close()
