
import os

import numpy as np
import pandas as pd

from scipy.stats import mannwhitneyu, ttest_ind

import matplotlib as mpl
# set no display
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'DejaVu Sans'
mpl.rcParams['pdf.fonttype'] = 42

col_pal = [
    mpl.colors.hex2color('#D21417'),
    mpl.colors.hex2color('#1483D2')
]


def plot_exp(fig, gs, i, query_exp, ref_exp, args, square_mk, circle_mk):
    (t_value, p_value_5) = ttest_ind(query_exp[1], ref_exp[1], equal_var=False)
    (t_value, p_value_100) = ttest_ind(query_exp[2], ref_exp[2], equal_var=False)

    df_swarn = pd.DataFrame(
        data={
            'x': np.append(query_exp[0], ref_exp[0]),
            'subgroup': np.append(
                np.repeat("Query", len(query_exp[0])),
                np.repeat("Other", len(ref_exp[0]))
            )
        }
    )

    ax_kde = fig.add_subplot(gs[0, i])
    ax_swarm1 = fig.add_subplot(gs[1:3, i], sharex=ax_kde)
    ax_swarm2 = fig.add_subplot(gs[4:7, i], sharex=ax_kde)

    fig.subplots_adjust(hspace=0)
    sns.despine(ax=ax_kde, top=True, right=True, left=True, bottom=True)
    sns.despine(ax=ax_swarm1, top=True, right=True, left=True, bottom=False)
    sns.despine(ax=ax_swarm2, top=True, right=True, left=True, bottom=False)

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
    plt.setp(ax_swarm1.get_xticklabels(), visible=False)
    plt.setp(ax_swarm1.get_yticklabels(), visible=False)
    plt.setp(ax_swarm1.yaxis.get_majorticklines(), visible=False)
    plt.setp(ax_swarm1.yaxis.get_minorticklines(), visible=False)
    plt.setp(ax_swarm1.xaxis.get_majorticklines(), visible=False)
    plt.setp(ax_swarm1.xaxis.get_minorticklines(), visible=False)
    ax_swarm1.xaxis.grid(False)
    ax_swarm1.yaxis.grid(False)

    plt.setp(ax_swarm2.get_yticklabels(), visible=False)
    plt.setp(ax_swarm2.yaxis.get_majorticklines(), visible=False)
    plt.setp(ax_swarm2.yaxis.get_minorticklines(), visible=False)
    ax_swarm1.yaxis.grid(False)

    sns_plot = sns.kdeplot(query_exp[0], shade=True,
                           color=col_pal[0], label="Query", ax=ax_kde)
    # sns_plot = sns.rugplot(query_exp, color = "r")
    sns_plot = sns.kdeplot(ref_exp[0], shade=True, color=col_pal[1],
                           label="Other", ax=ax_kde)
    # sns_plot = sns.rugplot(ref_exp, color = "b")
    sns_plot.set_title("t_test pvalue 5 =" + str(p_value_5) + "\nt_test pvalue 100 =" + str(p_value_100))

    df_swarn = pd.DataFrame(
        data={
            'x': np.append(query_exp[1], ref_exp[1]),
            'subgroup': np.append(
                np.repeat("Query", len(query_exp[1])),
                np.repeat("Other", len(ref_exp[1]))
            )
        }
    )

    if args.STRIP:
        sns_plot = sns.stripplot(
            x="x", y="subgroup", data=df_swarn, ax=ax_swarm1,
            size=args.SIZE, jitter=0.4,
            palette=sns.color_palette([col_pal[0], col_pal[1]])
        )
    else:
        sns_plot = sns.swarmplot(
            x="x", y="subgroup", data=df_swarn, ax=ax_swarm1, size=args.SIZE,
            palette=sns.color_palette([col_pal[0], col_pal[1]])
        )

    # Change shape in function of subgroup
    collections = sns_plot.collections
    unique_colors = [list(col_pal[0]) + [1], list(col_pal[1]) + [1]]
    markers = [circle_mk, square_mk]
    for collection in collections:
        paths = []
        for current_color in collection.get_facecolors():
            for possible_marker, possible_color in zip(markers, unique_colors):
                if np.array_equal(current_color, possible_color):
                    paths.append(possible_marker)
                    break
        collection.set_paths(paths)


    df_swarn = pd.DataFrame(
        data={
            'x': np.append(query_exp[2], ref_exp[2]),
            'subgroup': np.append(
                np.repeat("Query", len(query_exp[2])),
                np.repeat("Other", len(ref_exp[2]))
            )
        }
    )

    if args.STRIP:
        sns_plot = sns.stripplot(
            x="x", y="subgroup", data=df_swarn, ax=ax_swarm2,
            size=args.SIZE, jitter=0.4,
            palette=sns.color_palette([col_pal[0], col_pal[1]])
        )
    else:
        sns_plot = sns.swarmplot(
            x="x", y="subgroup", data=df_swarn, ax=ax_swarm2, size=args.SIZE,
            palette=sns.color_palette([col_pal[0], col_pal[1]])
        )

    # Change shape in function of subgroup
    collections = sns_plot.collections
    unique_colors = [list(col_pal[0]) + [1], list(col_pal[1]) + [1]]
    markers = [circle_mk, square_mk]
    for collection in collections:
        paths = []
        for current_color in collection.get_facecolors():
            for possible_marker, possible_color in zip(markers, unique_colors):
                if np.array_equal(current_color, possible_color):
                    paths.append(possible_marker)
                    break
        collection.set_paths(paths)

    x_label = "x"

    sns_plot.set(xlabel=x_label)

    return(sns_plot)


def main_diff_pred(args, argparser):

    fig_dir = os.path.join(args.OUTDIR, "diff_pred")

    if not os.path.exists(fig_dir):
        os.makedirs(fig_dir)


    # dummy plots, just to get the Path objects
    fig, ax = plt.subplots(1, 1)
    a = ax.scatter([1, 2], [3, 4], marker='s')
    b = ax.scatter([1, 2], [3, 4])
    square_mk, = a.get_paths()
    circle_mk, = b.get_paths()

    num_plot = 3
    fig = plt.figure(figsize=(5*num_plot, 7))
    gs = plt.GridSpec(7, num_plot)

    query_exp = [[] for x in range(0, num_plot)]
    ref_exp = [[] for x in range(0, num_plot)]

    query_exp[0] = np.random.normal(10, 1, size=10000)
    query_exp[1] = np.random.normal(10, 1, size=6)
    query_exp[2] = np.random.normal(10, 1, size=100)
    ref_exp[0] = np.random.normal(18, 1, size=10000)
    ref_exp[1] = np.random.normal(18, 1, size=6)
    ref_exp[2] = np.random.normal(18, 1, size=100)
    sns_plot = plot_exp(fig, gs, 0, query_exp, ref_exp, args, square_mk, circle_mk)

    query_exp[0] = np.random.normal(10, 4, size=10000)
    query_exp[1] = np.random.normal(10, 4, size=6)
    query_exp[2] = np.random.normal(10, 4, size=100)
    ref_exp[0] = np.random.normal(18, 4, size=10000)
    ref_exp[1] = np.random.normal(18, 4, size=6)
    ref_exp[2] = np.random.normal(18, 4, size=100)
    sns_plot = plot_exp(fig, gs, 1, query_exp, ref_exp, args, square_mk, circle_mk)

    query_exp[0] = np.append(np.random.normal(3, 0.2, size=5000), np.random.normal(7, 0.2, size=5000))
    query_exp[1] = np.append(np.random.normal(3, 0.2, size=3), np.random.normal(7, 0.2, size=3))
    query_exp[2] = np.append(np.random.normal(3, 0.2, size=50), np.random.normal(7, 0.2, size=50))
    ref_exp[0] = np.random.normal(5, 0.2, size=10000)
    ref_exp[1] = np.random.normal(5, 0.2, size=6)
    ref_exp[2] = np.random.normal(5, 0.2, size=100)
    sns_plot = plot_exp(fig, gs, 2, query_exp, ref_exp, args, square_mk, circle_mk)



    fig_out = os.path.join(fig_dir, "diff_pred.pdf")
    sns_plot.figure.savefig(fig_out)
    plt.close()
