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

def main_density(args, argparser):

    dir_design = os.path.join(args.PATH, args.DESIGN)
    df_design = uo.get_design(args, dir_design)

    df_exp = uo.get_exp(args, args.MATRIX)
    df_exp = df_exp.reindex(df_design["sample"], axis=1)

    num_query = len(np.where(df_design[args.SUBGROUP] == 1)[0])
    query_exp = df_exp.loc[args.ID][:num_query]
    ref_exp = df_exp.loc[args.ID][num_query:]

    bw = bw_nrd0(df_exp.loc[args.ID])
    bw_s = bw_scott(df_exp.loc[args.ID])
    print(bw)
    print(bw_s)
    sns_plot = sns.kdeplot(query_exp, shade=True, bw=bw, color = "r", label=args.DESIGN);
    sns_plot = sns.rugplot(query_exp, color = "r");
    sns_plot = sns.kdeplot(ref_exp, shade=True, bw=bw, color = "b", label="Other");
    sns_plot = sns.rugplot(ref_exp, color = "b");

    sns_plot.set_title(str(args.ID) + " " + args.DESIGN + "\nbw=" + str(bw))

    fig_dir = os.path.join(args.OUTDIR, "density", args.DESIGN)
    if not os.path.exists(fig_dir):
        os.makedirs(fig_dir)
    fig_file = os.path.join(fig_dir, args.ID + "_density.pdf")
    sns_plot.figure.savefig(fig_file)
    plt.close()
