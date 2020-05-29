

# get time and date
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

import datetime

start_time = datetime.datetime.now()

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'DejaVu Sans'
mpl.rcParams['pdf.fonttype'] = 42


def main_eval_ss(args, argparser):

    file_dict = {
        'epcy': "predictive_capability.xls",
        'deseq2': "deseq2_genes.xls",
        'edger': "edger_genes.xls",
        'voom': "limma_voom_genes.xls",
        'trend': "limma_trend_genes.xls",
    }

    # create dict with search params
    search_params_dict = {
        'reps': args.REPS,
        'designs': args.DESIGNS,
        'p_ss': args.P_SS,
        'methods': args.METHODS
    }

    dict_diff = defaultdict(list)
    for design in search_params_dict['designs']:
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
                    if method == "epcy":
                        dict_diff_tmp['value'] = df_diff["KERNEL_MCC"].tolist()
                    elif method == "deseq2":
                        dict_diff_tmp['value'] = (-np.log10(df_diff["padj"])).tolist()
                    elif method == "edger":
                        dict_diff_tmp['value'] = (-np.log10(df_diff["FDR"])).tolist()
                    elif method == "voom":
                        dict_diff_tmp['value'] = (-np.log10(df_diff["adj.P.Val"])).tolist()
                    elif method == "trend":
                        dict_diff_tmp['value'] = (-np.log10(df_diff["adj.P.Val"])).tolist()

                    for key in dict_diff_tmp.keys():
                        dict_diff[key] = dict_diff[key] + dict_diff_tmp[key]

    df_diff = pd.DataFrame(dict_diff)
    df_diff_std = df_diff[["method", "ID", "value"]]
    df_diff_std = df_diff_std.groupby(["method", "ID"]).std()
    df_diff_std = df_diff_std.reset_index()

    fig_dir = os.path.join(args.OUTDIR, "eval_ss")
    if not os.path.exists(fig_dir):
        os.makedirs(fig_dir)

    df_diff_std_epcy = df_diff_std[df_diff_std["method"] == "epcy"]
    df_diff_std_epcy = df_diff_std_epcy.rename(columns={"value": "epcy"})
    for method in search_params_dict['methods']:
        if method != "epcy":
            df_diff_std_DEG = df_diff_std[df_diff_std["method"] == method]
            df_diff_std_DEG = df_diff_std_DEG.rename(columns={"value": method})
            df_merge = df_diff_std_epcy.merge(
                df_diff_std_DEG, left_on='ID', right_on='ID'
            )
            sns_plot = sns.scatterplot(
                x="epcy", y=method,
                data=df_merge, size=args.SIZE
            )
            df_ids = df_merge.loc[df_merge['ID'].isin(args.IDS)]
            for index, row in df_ids.iterrows():
                sns_plot.text(
                    x=row["epcy"]+0.01,
                    y=row[method],
                    s=row["ID"],
                    horizontalalignment='left',
                    size='medium', color='black', weight='semibold'
                )

            fig_out = os.path.join(fig_dir, "std_epcy_vs_" + method + ".pdf")
            sns_plot.figure.savefig(fig_out)
            plt.close()


    file_out = os.path.join(fig_dir, "all.tsv")
    df_diff.to_csv(file_out, index=False, sep="\t")

    file_out = os.path.join(fig_dir, "std.tsv")
    df_diff_std.to_csv(file_out, index=False, sep="\t")

    for id in args.IDS:
        df_diff_id = df_diff[df_diff["ID"] == id]
        df_diff_id_epcy = df_diff_id[df_diff_id["method"] == "epcy"]
        df_diff_id_deg = df_diff_id[df_diff_id["method"] != "epcy"]

        sns_plot = sns.swarmplot(
            data=df_diff_id_epcy,
            x="cohort", y="value", color=".2", size=args.SIZE
        )
        sns_plot.set_title("EPCY " + id)
        sns_plot.set(ylim=(-1.05, 1.05))
        sns_plot.set(ylabel="MCC")
        sns_plot.set_xticklabels(sns_plot.get_xticklabels(), rotation=90)
        fig_out = os.path.join(fig_dir, id + "_epcy.pdf")
        sns_plot.figure.savefig(fig_out)
        plt.close()

        sns_plot = sns.swarmplot(
            data=df_diff_id_deg,
            x="cohort", y="value", hue="method", dodge=True, size=args.SIZE
        )

        sns_plot.set_title("DEG " + id)
        sns_plot.axes.axhline(-np.log10(0.05), ls='--', color="r")
        sns_plot.set(ylabel="-log10(padj)")
        sns_plot.set_xticklabels(sns_plot.get_xticklabels(), rotation=90)
        fig_out = os.path.join(fig_dir, id + "_DEG.pdf")
        sns_plot.figure.savefig(fig_out)
        plt.close()
