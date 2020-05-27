import datetime
# get time and date
start_time = datetime.datetime.now()
# make rest of imports
import os

from collections import defaultdict

import pdb

import pandas as pd
import numpy as np

import matplotlib as mpl
# set no display
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

from .. utils import other as uo

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'DejaVu Sans'
mpl.rcParams['pdf.fonttype'] = 42

def main_eval_random(args, argparser):

    file_dict = {
        'epcy' : "predictive_capability.xls",
        'deseq2' : "deseq2_genes.xls",
        'edger' : "edger_genes.xls",
        'limma' : "limma_voom_genes.xls",
        'trend' : "limma_trend_genes.xls",
    }

    # create dict with search params
    search_params_dict = {
        'methods' : args.METHODS,
        'designs_random' : args.DESIGN_RANDOM,
        'designs' : args.DESIGN
    }

    strl2fc_dict = {
        'epcy' : "L2FC",
        'deseq2' : "log2FoldChange",
        'edger' : "logFC",
        'limma' : "logFC",
        'trend' : "logFC",
    }

    strvalue_dict = {
        'epcy' : "KERNEL_MCC",
        'deseq2' : "pvalue",
        'edger' : "PValue",
        'limma' : "P.Value",
        'trend' : "P.Value",
    }

    strvalueadj_dict = {
        'epcy' : "AUC",
        'deseq2' : "padj",
        'edger' : "FDR",
        'limma' : "adj.P.Val",
        'trend' : "adj.P.Val",
    }

    top = 50

    df_biotype = None
    if args.BIOTYPE is not None:
        df_biotype = pd.read_csv(args.BF, sep="\t")
        selected_biotype = args.BIOTYPE.split(",")
        df_biotype = df_biotype.loc[df_biotype["gene_biotype"].isin(selected_biotype)]

    df_all_res = []
    for design in search_params_dict['designs_random']:
        dir_design = os.path.join(args.RANDOM_PATH, design)
        df_design = uo.get_design(args, "Query", dir_design)

        for method in search_params_dict['methods']:
            print('RANDOM design: {}, method: {}'.format(design, method))

            df_tmp = uo.read_diff_table(args, file_dict[method], method, dir_design, 1)

            if method != "epcy":
                df_tmp = df_tmp.reindex(df_tmp[strvalue_dict[method]].abs().sort_values(ascending=True).index)

            df_all_res.append(
                pd.DataFrame({
                    "ID": df_tmp["ID"][:top],
                    "METHOD": method,
                    "DESIGN": design,
                    "TYPE": "MCC" if method == "epcy" else "PVALUE",
                    "VALUE": df_tmp[strvalue_dict[method]][:top] if method == "epcy" else -np.log10(df_tmp[strvalue_dict[method]][:top]),
                    "L2FC": df_tmp[strl2fc_dict[method]][:top]
                })
            )

            if method != "epcy":
                df_tmp = df_tmp.reindex(df_tmp[strvalueadj_dict[method]].abs().sort_values(ascending=True).index)

                df_all_res.append(
                    pd.DataFrame({
                        "ID": df_tmp["ID"][:top],
                        "METHOD": method,
                        "DESIGN": design,
                        "TYPE": "PADJ",
                        "VALUE": -np.log10(df_tmp[strvalueadj_dict[method]][:top]),
                        "L2FC": df_tmp[strl2fc_dict[method]][:top]
                    })
                )

    for design in search_params_dict['designs']:
        dir_design = os.path.join(args.PATH, design)
        df_design = uo.get_design(args, "Query", dir_design)

        for method in search_params_dict['methods']:
            print('RANDOM design: {}, method: {}'.format(design, method))

            df_tmp = uo.read_diff_table(args, file_dict[method], method, dir_design, 1)

            if method != "epcy":
                df_tmp = df_tmp.reindex(df_tmp[strvalue_dict[method]].abs().sort_values(ascending=True).index)

            df_all_res.append(
                pd.DataFrame({
                    "ID": df_tmp["ID"][:top],
                    "METHOD": method,
                    "DESIGN": design,
                    "TYPE": "MCC" if method == "epcy" else "PVALUE",
                    "VALUE": df_tmp[strvalue_dict[method]][:top] if method == "epcy" else -np.log10(df_tmp[strvalue_dict[method]][:top]),
                    "L2FC": df_tmp[strl2fc_dict[method]][:top]
                })
            )

            if method != "epcy":
                df_tmp = df_tmp.reindex(df_tmp[strvalueadj_dict[method]].abs().sort_values(ascending=True).index)

                df_all_res.append(
                    pd.DataFrame({
                        "ID": df_tmp["ID"][:top],
                        "METHOD": method,
                        "DESIGN": design,
                        "TYPE": "PADJ",
                        "VALUE": -np.log10(df_tmp[strvalueadj_dict[method]][:top]),
                        "L2FC": df_tmp[strl2fc_dict[method]][:top]
                    })
                )

    print('MERGE ALL DATA')
    df_all_res = pd.concat(df_all_res)
    # create subfolder for log and res
    fig_dir = os.path.join(args.OUTDIR, "eval_random")

    if not os.path.exists(fig_dir):
        os.makedirs(fig_dir)

    print('SAVE TABLE')
    csv_out = os.path.join(fig_dir, "random_table.csv")
    df_all_res.to_csv(csv_out, index=False, sep="\t")

    print('PLOT FIG')

    for method in search_params_dict['methods']:
        df_tmp = df_all_res.loc[df_all_res["METHOD"] == method]
        if method != "epcy":
            # plot results pvalue
            df_tmp2 = df_tmp.loc[df_tmp["TYPE"] == "PVALUE"]
            plt_fig = sns.swarmplot(x="DESIGN", y="VALUE", data=df_tmp2, hue="METHOD")
            plt_fig.set(ylabel="-log10(pvalue)")
            #plt_fig.axhline(y=-np.log10(0.05), color='r', linestyle='--')
            plt_fig.set_xticklabels(plt_fig.get_xticklabels(), rotation=90)
            fig_out = os.path.join(fig_dir, method + "_pvalue.pdf")
            plt_fig.figure.savefig(fig_out)
            plt.close()

            # plot results pvalue adj
            df_tmp2 = df_tmp.loc[df_tmp["TYPE"] == "PADJ"]
            plt_fig = sns.swarmplot(x="DESIGN", y="VALUE", data=df_tmp2, hue="METHOD")
            plt_fig.set(ylabel="-log10(p_adjusted)")
            #plt_fig.axhline(y=-np.log10(0.05), color='r', linestyle='--')
            plt_fig.set_xticklabels(plt_fig.get_xticklabels(), rotation=90)
            fig_out = os.path.join(fig_dir, method + "_padj.pdf")
            plt_fig.figure.savefig(fig_out)
            plt.close()
        else:
            # plot results pvalue
            df_tmp2 = df_tmp.loc[df_tmp["TYPE"] == "MCC"]
            plt_fig = sns.swarmplot(x="DESIGN", y="VALUE", data=df_tmp2, hue="METHOD")
            plt_fig.set(ylabel="MCC")
            #plt_fig.axhline(y=0.2, color='r', linestyle='--')
            plt_fig.set_xticklabels(plt_fig.get_xticklabels(), rotation=90)
            fig_out = os.path.join(fig_dir, method + "_mcc.pdf")
            plt_fig.figure.savefig(fig_out)
            plt.close()
