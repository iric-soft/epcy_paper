

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

                if method == "epcy" or method == "epcy_bagging":
                    dict_diff_tmp['l2fc'] = df_diff["L2FC"].tolist()
                elif method == "deseq2":
                    dict_diff_tmp['l2fc'] = df_diff["log2FoldChange"].tolist()
                elif method == "edger":
                    dict_diff_tmp['l2fc'] = df_diff["logFC"].tolist()
                elif method == "voom":
                    dict_diff_tmp['l2fc'] = df_diff["logFC"].tolist()
                elif method == "trend":
                    dict_diff_tmp['l2fc'] = df_diff["logFC"].tolist()
                elif method == "mast":
                    dict_diff_tmp['l2fc'] = [np.nan for x in range(1, df_diff.shape[0]+1, 1)]

                dict_diff_tmp['top'] = [x for x in range(1, df_diff.shape[0]+1, 1)]

                for key in dict_diff_tmp.keys():
                    dict_diff[key] = dict_diff[key] + dict_diff_tmp[key]

    df_diff = pd.DataFrame(dict_diff)

    df_diff_mean = df_diff.groupby(["method", "cohort", "ID"])["value"].mean()
    df_diff_mean = df_diff_mean.to_frame()
    df_diff_mean = df_diff_mean.rename({"value": "mean_value"}, axis='columns')
    df_diff_mean.reset_index(inplace=True)

    df_diff_std = df_diff.groupby(["method", "cohort", "ID"])["value"].std()
    df_diff_std = df_diff_std.to_frame()
    df_diff_std = df_diff_std.rename({"value": "y_std"}, axis='columns')
    df_diff_std.reset_index(inplace=True)

    df_diff_l2fc = df_diff.groupby(["method", "cohort", "ID"])["l2fc"].mean()
    df_diff_l2fc = df_diff_l2fc.to_frame()
    df_diff_l2fc = df_diff_l2fc.rename({"l2fc": "mean_l2fc"}, axis='columns')
    df_diff_l2fc.reset_index(inplace=True)

    df_diff_resume = pd.merge(df_diff_mean, df_diff_std, on=['method', 'cohort', 'ID'])
    df_diff_resume = pd.merge(df_diff_resume, df_diff_l2fc, on=['method', 'cohort', 'ID'])


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

        fig = plt.figure(figsize=(5, 5))
        gs = plt.GridSpec(1, 2)

        ax_deg = fig.add_subplot(gs[0, 0])
        ax_epcy = fig.add_subplot(gs[0, 1], sharex=ax_kde)

        fig.subplots_adjust(hspace=0)
        sns.despine(ax=ax_deg, top=True, right=True, left=True, bottom=True)
        sns.despine(ax=ax_epcy, top=True, right=True, left=True, bottom=False)


        sns_plot = sns.pointplot(
            data=df_diff_id_epcy,
            y="cohort", x="value", hue="method", dodge=0.3,
            linestyles=[':'],
            markers=['x']
        )
        sns_plot.set_title("EPCY " + id)
        sns_plot.set(xlim=(-1.05, 1.05))
        sns_plot.set(xlabel="MCC")
        sns_plot.set_xticklabels(sns_plot.get_xticklabels(), rotation=90)
        fig_out = os.path.join(fig_dir, id + "_epcy.pdf")
        sns_plot.figure.savefig(fig_out)
        plt.close()

        sns_plot = sns.pointplot(
            data=df_diff_id_deg,
            y="cohort", x="value", hue="method", dodge=0.3,
            linestyles=['-', '--', '-.'],
            markers=['o', 'v', 's']
        )

        sns_plot.set_title("DEG " + id)
        sns_plot.axes.axvline(-np.log10(0.05), ls='--', color="b")
        sns_plot.set(xlabel="-log10(padj)")
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

    if len(cohort_order) == 10:
        cohort_order = [cohort_order[i] for i in [0, 2, 9]]
        df_diff_resume = df_diff_resume[df_diff_resume.isin(["3_vs_3", "8_vs_127", "27_vs_563"])]

    for method in search_params_dict['methods']:
        df_diff_method = df_diff_resume[(df_diff_resume["method"] == method)]
        df_diff_method["group"] = "Query"
        df_diff_method["group"].loc[df_diff_method["mean_l2fc"] < 0] = "Ref"


        plt_fig = sns.relplot(
            x="mean_l2fc", y="mean_value", hue="group",
            row="method", col="cohort", size="y_std", sizes=(10, 200),
            col_order=cohort_order,
            data=df_diff_method
        )

        fig_out = os.path.join(fig_dir, method + "_volcano.pdf")
        plt_fig.savefig(fig_out)
        plt.close()
