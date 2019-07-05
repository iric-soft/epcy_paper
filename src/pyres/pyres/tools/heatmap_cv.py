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
from sklearn.utils import shuffle as skl_shuffle

import logging

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'DejaVu Sans'
mpl.rcParams['pdf.fonttype'] = 42

def main_heatmap_cv(args, argparser):

    file_dict = {
        'epcy' : "prediction_capability.xls",
        'deseq2' : "deseq2_genes.xls",
        'edger' : "edger_genes.xls",
        'limma' : "limma_voom_genes.xls",
    }

    criteria_dict = {
        'epcy' : "KERNEL_MCC",
        'deseq2' : "pvalue",
        'edger' : "PValue",
        'limma' : "P.Value",
    }

    # set fold argument depending of input param for fold analysis
    fold = "fold" + str(args.FOLD) if not args.LOO else 'foldloo'

    # create dict with search params
    search_params_dict = {
        'designs' : args.DESIGN
    }
    method = args.METHOD
    top = args.TOP

    df_biotype = pd.read_csv(args.BF, sep="\t")
    if (args.BIOTYPE is not None):
        selected_biotype = args.BIOTYPE.split(",")
        df_biotype = df_biotype.loc[df_biotype["gene_biotype"].isin(selected_biotype)]

    for design in search_params_dict['designs']:
        dir_design = os.path.join(args.PATH, design, "cv", fold)
        datasets = os.listdir(dir_design)
        n_datasets =  int(args.N_DATASETS) if args.N_DATASETS != 'all' else len(datasets)

        result = None
        for dataset in datasets[:n_datasets]:
            path_dataset_train = os.path.join(dir_design, dataset, "train")
            data = uo.read_diff_table(args, file_dict[method], method, df_biotype["ensembl_gene_id"], path_dir=path_dataset_train)

            tmp = data[["ID",criteria_dict[method]]][:top]
            tmp.columns = ["ID", dataset]
            tmp[dataset] = list(range(1, top + 1))
            #print(tmp)
            #if method != "epcy":
            #    tmp[dataset] = -log10(tmp[dataset])

            if result is None:
                result = pd.DataFrame(tmp)
            else:
                result = result.merge(tmp, how='outer', on="ID")

        #replace missing values (features not present in all dataset)
        result.fillna(top + 1 , inplace=True)
        result = result.set_index("ID")
        max_value = np.amax(result.values)
        min_value = np.amin(result.values)

        sns_plot = sns.clustermap(
            result, linewidths=0, metric='euclidean',
            xticklabels=False, yticklabels=False,
            vmin=min_value, vmax=max_value#, cmap="Blues"
        )
        sns_plot.fig.suptitle(method + " " + design + ": " + str(len(result.index)) + " genes")

        fig_dir = os.path.join(args.OUTDIR, "figures", design)
        if not os.path.exists(fig_dir):
            os.makedirs(fig_dir)
        fig_file = os.path.join(fig_dir, method + "_" + str(top) + "_heatmap_cv.pdf")
        sns_plot.savefig(fig_file)
        plt.close()
