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

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'DejaVu Sans'
mpl.rcParams['pdf.fonttype'] = 42

def main_eval_cv(args, argparser):

    file_dict = {
        'epcy' : "predictive_capability.tsv",
        'deseq2' : "deseq2_genes.xls",
        'edger' : "edger_genes.xls",
        'limma' : "limma_voom_genes.xls",
    }

    # set fold argument depending of input param for fold analysis
    fold = "fold" + str(args.FOLD) if not args.LOO else 'foldloo'

    # create dict with search params
    search_params_dict = {
        'tops' : args.TOP_VALUES,
        'methods' : args.METHODS,
        'designs' : args.DESIGN,
        'pvalues' : args.PVALUES
    }


    df_exp = uo.get_exp(args, args.MATRIX)
    list_genes = df_exp.index.tolist()

    dict_diff = defaultdict(list)
    dict_pred = defaultdict(list)
    dict_res = defaultdict(list)

    for design in search_params_dict['designs']:
        dir_design = os.path.join(args.PATH, design, "cv", fold)
        datasets = os.listdir(dir_design)
        n_datasets =  int(args.N_DATASETS) if args.N_DATASETS != 'all' else len(datasets)

        for method in search_params_dict['methods']:

            for pvalue in search_params_dict['pvalues']:
                print('TRAINING/TESTING design: {}, method: {}, pvalue: {}, n_datasets: {} / {}'.format(design, method, pvalue, n_datasets, len(datasets)))

                for dataset in datasets[:n_datasets]:
                    path_dataset_train = os.path.join(dir_design, dataset, "train")
                    df_design_train = uo.get_design(args, "Query", path_dataset_train)

                    path_dataset_test = os.path.join(dir_design, dataset, "test")
                    df_design_test = uo.get_design(args, "Query", path_dataset_test)

                    dict_diff[method] = uo.read_diff_table(args, file_dict[method], method, path_dataset_train, pvalue, list_genes)

                    Nmax = len(dict_diff[method])
                    print("##{}, {}, {} ## NMAX: {} ".format(design, dataset.upper(), method.upper(), Nmax))

                    for top in search_params_dict['tops']:
                        # define our features
                        df_feature = df_exp.loc[df_exp.index.isin(dict_diff[method]["ID"][:top])]

                        df_exp_train = df_feature[df_design_train["sample"]]
                        df_label_train = df_design_train[args.SUBGROUP]

                        df_exp_test = df_feature[df_design_test["sample"]]
                        df_label_test = df_design_test[args.SUBGROUP]

                        if args.USE_LR:
                            dict_pred_tmp = uo.run_LR(df_exp_train, df_label_train, df_exp_test, df_label_test, design, method, dataset, top, pvalue)
                            for key in dict_pred_tmp.keys():
                                dict_pred[key] = dict_pred[key] + dict_pred_tmp[key]

                        if args.USE_RANDF:
                            dict_randF_tmp = uo.run_rand_forest(df_exp_train, df_label_train, df_exp_test, df_label_test, design, method, dataset, top, pvalue)
                            for key in dict_randF_tmp.keys():
                               dict_pred[key] = dict_pred[key] + dict_randF_tmp[key]

                        # add a seed column to be consistent with shuffle part
                        dict_pred['seed'] = dict_pred['seed'] + [float('NaN')]

        # tranform dict of results to df
        df_LR = pd.DataFrame(dict_pred)

    ##############################################################################################################
    ###################### COMPUTE PERF AND PLOT ##################################################################
    ##############################################################################################################

    #init contengency table for mcc computing
    ct = [0,0,0,0]
    for design in search_params_dict['designs']:
        for method in search_params_dict['methods']:
            for top in search_params_dict['tops']:
                for pvalue in search_params_dict['pvalues']:
                    for learn in df_LR.learn.unique() :
                        df = df_LR.loc[df_LR.design == design]
                        df = df.loc[df.top == top]
                        df = df.loc[df.method == method]
                        df = df.loc[df.learn == learn]
                        df = df.loc[df.pvalue == pvalue]

                        print(df)

                        # continue only if df is non-empty
                        if len(df) < 2 : break

                        df = df.sort_values(by=['dataset'], ascending=[True])

                        # convert probs to array
                        probas = df.proba.values
                        # same for labels
                        labels = np.array([label for x in df.label.values for label in x])
                        ids_query = np.where(labels == 1)
                        ids_ref = np.where(labels == 0)

                        # compute mcc
                        ct[0] = np.sum(probas[ids_query] > 0.5)
                        ct[1] = np.sum(probas[ids_query] <= 0.5)
                        ct[3] = np.sum(probas[ids_ref] <= 0.5)
                        ct[2] = np.sum(probas[ids_ref] > 0.5)
                        mcc = uo.get_mcc(ct)

                        auc = uo.auc_u_test(probas, ids_query, ids_ref)

                        dict_res['auc'].append(auc[0])
                        dict_res['mcc'].append(mcc)
                        dict_res['method'].append(method)
                        dict_res['top'].append(top)
                        dict_res['design'].append(design)
                        dict_res['learn'].append(learn)
                        dict_res['pvalue'].append(pvalue)


    # create subfolder for log and res
    fig_dir = os.path.join(args.OUTDIR, "eval_cv", "L_Reg")
    if args.USE_RANDF:
        fig_dir = os.path.join(args.OUTDIR, "eval_cv", "rand_forest")

    if not os.path.exists(fig_dir):
        os.makedirs(fig_dir)

    # convert dict to df for analysis
    df_res = pd.DataFrame(dict_res)
    csv_out = os.path.join(fig_dir, "pred_table.csv")
    df_res.to_csv(csv_out, index=False, sep="\t")

    # plot results
    plt_fig = sns.catplot(data=df_res, x="top", y="mcc", hue="method", col="design", row = "pvalue", kind="point", facet_kws=dict(subplot_kws=dict(ylim=[-0.05,1.05])))
    fig_out =  os.path.join(fig_dir, "pred_mcc.pdf")
    plt_fig.savefig(fig_out)
    plt.close()

    #plt_fig = sns.catplot(data=df_res, x="top", y="auc", hue="method", col="design", row = "pvalue", kind="point", facet_kws=dict(subplot_kws=dict(ylim=[-0.05,1.05])))
    #fig_out =  os.path.join(fig_dir, "pred_auc.pdf")
    #plt_fig.savefig(fig_out)
    #plt.close()
