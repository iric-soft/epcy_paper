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

from statsmodels.nonparametric.bandwidths import bw_scott

from .. utils import other as uo
from sklearn.utils import shuffle as skl_shuffle

import logging

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'DejaVu Sans'
mpl.rcParams['pdf.fonttype'] = 42

#COPY/PASTE from EPCY
def bw_nrd0(x):
    #TODO need to improve speed of this part
    hi = np.std(x)
    iqr = np.subtract(*np.percentile(x, [75, 25]))
    lo = min(hi, iqr/1.34)
    if (lo == 0):
        lo = hi
        if (lo == 0):
            lo = abs(x[1])
            if (lo == 0):
                lo = 1

    # this function can be run by ne.evaluate, with all lo pre-computed
    return(0.9 * lo * len(x)**(-0.2))

def main_gene_density(args, argparser):

    col_pal = [
        mpl.colors.hex2color('#E62528'),
        mpl.colors.hex2color('#0D6BAC')
    ]

    file_dict = {
        'epcy' : "prediction_capability.xls",
        'deseq2' : "deseq2_genes.xls",
        'edger' : "edger_genes.xls",
        'limma' : "limma_voom_genes.xls",
    }

    dir_design = os.path.join(args.PATH, args.DESIGN)
    df_design = uo.get_design(args, args.QUERY, path=dir_design)

    df_exp = uo.get_exp(args, args.MATRIX)
    df_exp = df_exp.reindex(df_design["sample"], axis=1)
    list_genes = df_exp.index.tolist()

    num_query = len(np.where(df_design[args.SUBGROUP] == 1)[0])

    list_ids = [args.GENE_ID]

    f = plt.figure(figsize=(5, 5))
    gs = plt.GridSpec(4, 1)

    bw = bw_nrd0(df_exp.loc[args.GENE_ID])
    query_exp = df_exp.loc[args.GENE_ID][:num_query]
    ref_exp = df_exp.loc[args.GENE_ID][num_query:]

    df_swarn = pd.DataFrame(
        data={
            'log2(x+1)': np.append(query_exp, ref_exp),
            'subgroup': np.append(
                np.repeat(args.DESIGN, len(query_exp)),
                np.repeat("Other", len(ref_exp))
            )
        }
    )

    ax_kde = f.add_subplot(gs[0, 0])
    ax_swarm = f.add_subplot(gs[1:, 0], sharex=ax_kde)

    f.subplots_adjust(hspace=0)
    sns.despine(ax=ax_kde, top=True, right=True, left=True, bottom=True)
    sns.despine(ax=ax_swarm, top=True, right=True, left=True, bottom=False)

    # Turn off kde axis visibility
    plt.setp(ax_kde.get_xticklabels(), visible=False)
    plt.setp(ax_kde.get_yticklabels(), visible=False)
    plt.setp(ax_kde.yaxis.get_majorticklines(), visible=False)
    plt.setp(ax_kde.yaxis.get_minorticklines(), visible=False)
    plt.setp(ax_kde.xaxis.get_majorticklines(), visible=False)
    plt.setp(ax_kde.xaxis.get_minorticklines(), visible=False)
    ax_kde.yaxis.grid(False)
    ax_kde.xaxis.grid(False)

    # Turn off swarm axis visibility
    plt.setp(ax_swarm.get_yticklabels(), visible=False)
    plt.setp(ax_swarm.yaxis.get_majorticklines(), visible=False)
    plt.setp(ax_swarm.yaxis.get_minorticklines(), visible=False)
    ax_swarm.yaxis.grid(False)

    sns_plot = sns.kdeplot(query_exp, shade=True, bw=bw, color = col_pal[0], label=args.DESIGN, ax=ax_kde)
    #sns_plot = sns.rugplot(query_exp, color = "r")
    sns_plot = sns.kdeplot(ref_exp, shade=True, bw=bw, color = col_pal[1], label="Other", ax=ax_kde)
    #sns_plot = sns.rugplot(ref_exp, color = "b")
    sns_plot.set_title(str(args.GENE_ID) + " " + args.DESIGN + "\nbw=" + str(bw))

    sns_plot = sns.swarmplot(
        x="log2(x+1)", y="subgroup", data=df_swarn, ax=ax_swarm,
        palette=sns.color_palette([col_pal[0], col_pal[1]])
    )

    fig_dir = os.path.join(args.OUTDIR)
    if not os.path.exists(fig_dir):
        os.makedirs(fig_dir)
    fig_file = os.path.join(fig_dir, args.GENE_ID + "_density.pdf")
    sns_plot.figure.savefig(fig_file)
    plt.close()
