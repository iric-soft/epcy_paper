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
from sklearn.utils import shuffle as skl_shuffle

from sklearn import linear_model
from scipy.special import expit

import logging

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'DejaVu Sans'
mpl.rcParams['pdf.fonttype'] = 42

def main_log_reg(args, argparser):

    dir_design = os.path.join(args.PATH, args.DESIGN)
    df_design = uo.get_design(args, dir_design)

    df_feature = uo.get_exp(args, args.MATRIX)
    df_feature = df_feature.reindex(df_design["sample"], axis=1)
    df_feature = df_feature[df_design["sample"]]

    num_query = len(np.where(df_design[args.SUBGROUP] == 1)[0])
    query_exp = df_feature.loc[args.ID][:num_query]
    ref_exp = df_feature.loc[args.ID][num_query:]

    X = df_feature.loc[[args.ID]]
    X = X.T
    y = df_design[args.SUBGROUP]

    ########################################################
    #FROM:
    #https://scikit-learn.org/stable/auto_examples/linear_model/plot_logistic.html#sphx-glr-auto-examples-linear-model-plot-logistic-py
    ########################################################
    # Fit the classifier
    clf = linear_model.LogisticRegression(solver="liblinear", C = 1e10, penalty = 'l2', max_iter=1000)
    clf.fit(X, y)

    # and plot the result
    #plt.figure(1, figsize=(4, 3))
    #plt.clf()
    #plt.scatter(X.values.ravel(), y, color='black', zorder=20)
    sns_plot = sns.scatterplot(X.values.ravel(), y, color='black')

    min_x = min(X.values.ravel())
    max_x = max(X.values.ravel())
    print(min_x)
    print(max_x)
    X_test = np.linspace(min_x, max_x, 300)

    loss = expit(X_test * clf.coef_ + clf.intercept_).ravel()
    #plt.plot(X_test, loss, color='red', linewidth=3)
    sns_plot = sns.lineplot(X_test, loss, color='r')#, kwargs=dict(linewidth=3))

    sns_plot.axhline(.5, color='.5')

    #plt.xticks(range(-5, 10))
    #plt.yticks([0, 0.5, 1])
    #plt.ylim(-.25, 1.25)
    #plt.xlim(-4, 10)
    #plt.legend(('Logistic Regression Model', 'Linear Regression Model'),
    #           loc="lower right", fontsize='small')
    #plt.tight_layout()

    fig_dir = os.path.join(args.OUTDIR, "log_reg", args.DESIGN)
    if not os.path.exists(fig_dir):
        os.makedirs(fig_dir)
    fig_file = os.path.join(fig_dir, args.ID + "_log_reg.pdf")
    sns_plot.figure.savefig(fig_file)
    plt.close()
