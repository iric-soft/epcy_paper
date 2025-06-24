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


def main_eval_auc(args, argparser):

    file_dict = {
        'epcy' : "predictive_capability.tsv",
        'deseq2' : "deseq2_genes.xls",
        'edger' : "edger_genes.xls",
        'limma' : "limma_voom_genes.xls",
        'deseq2_pvalue': "deseq2_genes.xls",
        'edger_pvalue': "edger_genes.xls",
        'limma_pvalue': "limma_voom_genes.xls",
    }

    # create dict with search params
    search_params_dict = {
        'tops' : args.TOP_VALUES,
        'methods' : args.METHODS,
        'designs' : args.DESIGN,
        'padj' : args.PADJ
    }

    df_biotype = None
    if args.BIOTYPE is not None:
        df_biotype = pd.read_csv(args.BF, sep="\t")
        selected_biotype = args.BIOTYPE.split(",")
        df_biotype = df_biotype.loc[df_biotype["gene_biotype"].isin(selected_biotype)]

    dict_diff = defaultdict(list)
    dict_pred = defaultdict(list)
    dict_res = defaultdict(list)

    for design in search_params_dict['designs']:
        dir_design = os.path.join(args.PATH, design)
        df_design = uo.get_design(args, "Query", dir_design)

        for padj in search_params_dict['padj']:
            for method in search_params_dict['methods']:
                print('TRAINING/TESTING design: {}, method: {}, padj: {}'.format(design, method, padj))

                df_tmp = uo.read_diff_table(args, file_dict[method], method, dir_design, padj)
                if df_biotype is not None:
                    dict_diff[method] = df_tmp.loc[df_tmp["ID"].isin(df_biotype["ensembl_gene_id"])]
                else:
                    dict_diff[method] = df_tmp

            df_auc = dict_diff["epcy"]
            df_auc = df_auc.sort_values(["AUC", "abs_L2FC"], ascending=False)
            for top in search_params_dict['tops']:
                for method in search_params_dict['methods']:
                    top_ids = dict_diff[method]["ID"][:top]
                    df_top_auc = df_auc.loc[df_auc["ID"].isin(top_ids)]

                    dict_res['mean_auc'].append(np.mean(df_top_auc["AUC"]))
                    dict_res['method'].append(method)
                    dict_res['top'].append(top)
                    dict_res['padj'].append(padj)
                    dict_res['design'].append(design)

                dict_res['mean_auc'].append(np.mean(df_auc["AUC"][:top]))
                dict_res['method'].append("auc")
                dict_res['top'].append(top)
                dict_res['padj'].append(padj)
                dict_res['design'].append(design)
    # tranform dict of results to df
    df_auc = pd.DataFrame(dict_res)

    # create subfolder for log and res
    fig_dir = os.path.join(args.OUTDIR, "eval_auc")

    if not os.path.exists(fig_dir):
        os.makedirs(fig_dir)

    csv_out = os.path.join(fig_dir, "auc_table.csv")
    df_auc.to_csv(csv_out, index=False, sep="\t")

    # plot results
    plt_fig = sns.catplot(data=df_auc, x="top", y="mean_auc", hue="method", col="design", row = "padj", kind="point", facet_kws=dict(subplot_kws=dict(ylim=[0.45,1.05])))
    fig_out = os.path.join(fig_dir, "pred_auc.pdf")
    plt_fig.savefig(fig_out)
    plt.close()
