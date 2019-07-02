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

def main_eval_cv(args, argparser):

    # log name format: DD_MM_YY.log
    logname = '{}_{}_{}.log'.format(start_time.day, start_time.month, start_time.year)

    # create subfolder for log and res
    fig_dir = os.path.join(args.OUTDIR, 'figures')
    if not os.path.exists(fig_dir):
        os.makedirs(fig_dir)

    log_dir = os.path.join(args.OUTDIR, 'logs')
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)

    csv_dir = os.path.join(args.OUTDIR, 'tables')
    if not os.path.exists(csv_dir):
        os.makedirs(csv_dir)

    # init log
    logfile = os.path.join(log_dir, logname)
    logging.basicConfig(filename=logfile,level=logging.DEBUG)
    logging.info('###### BEGIN RUN - TIME: {} \t#####'.format(start_time.strftime("%Y-%m-%d %H:%M")))
    logging.info(str(args))

    file_dict = {
        'epcy' : "prediction_capability.xls",
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
        'designs' : args.DESIGN
    }


    df_exp = uo.get_exp(args, args.MATRIX)
    list_genes = df_exp.index.tolist()

    warmup = datetime.datetime.now()
    logging.info('Search parameters: {}'.format(search_params_dict))
    logging.info('WARM UP ELAPSED TIME : {} seconds\t'.format((warmup - start_time).seconds))

    dict_diff = defaultdict(list)
    dict_pred = defaultdict(list)
    dict_res = defaultdict(list)

    dict_diff_SHUFFLED = defaultdict(list)
    dict_pred_SHUFFLED = defaultdict(list)
    dict_res_SHUFFLED = defaultdict(list)

    for design in search_params_dict['designs']:
        dir_design = os.path.join(args.PATH, "design", design, "cv", fold)
        datasets = os.listdir(dir_design)
        n_datasets =  int(args.N_DATASETS) if args.N_DATASETS != 'all' else len(datasets)

        for method in search_params_dict['methods']:
            logging.info('TRAINING/TESTING design: {}, method: {}, n_datasets: {} / {}'.format(design, method, n_datasets, len(datasets)))
            print('TRAINING/TESTING design: {}, method: {}, n_datasets: {} / {}'.format(design, method, n_datasets, len(datasets)))
            checkpoint =  datetime.datetime.now()

            for dataset in datasets[:n_datasets]:
                path_dataset_train = os.path.join(dir_design, dataset, "train")
                df_design_train = uo.get_design(args, path_dataset_train)

                path_dataset_test = os.path.join(dir_design, dataset, "test")
                df_design_test = uo.get_design(args, path_dataset_test)

                dict_diff[method] = uo.read_diff_table(args, file_dict[method], method, list_genes, path_dir=path_dataset_train)

                Nmax = len(dict_diff[method])
                logging.info("##{}, {}, {} ## NMAX: {} ".format(design, dataset.upper(),method.upper(), Nmax))
                print("##{}, {}, {} ## NMAX: {} ".format(design, dataset.upper(), method.upper(), Nmax))

                for top in search_params_dict['tops']:
                    # define our features
                    df_feature = df_exp.loc[df_exp.index.isin(dict_diff[method]["ID"][:top])]

                    df_exp_train = df_feature[df_design_train["sample"]]
                    df_label_train = df_design_train[args.SUBGROUP]

                    df_exp_test = df_feature[df_design_test["sample"]]
                    df_label_test = df_design_test[args.SUBGROUP]

                    if args.USE_LR:
                        dict_pred_tmp = uo.run_LR(df_exp_train, df_label_train, df_exp_test, df_label_test, design, method, dataset, top)
                        for key in dict_pred_tmp.keys():
                            dict_pred[key] = dict_pred[key] + dict_pred_tmp[key]

                    if args.USE_RANDF:
                        dict_randF_tmp = uo.run_rand_forest(df_exp_train, df_label_train, df_exp_test, df_label_test, design, method, dataset, top)
                        for key in dict_randF_tmp.keys():
                           dict_pred[key] = dict_pred[key] + dict_randF_tmp[key]

                    # add a seed column to be consistent with shuffle part
                    dict_pred['seed'] = dict_pred['seed'] + [float('NaN')]

                # only if specified
                if args.SHUFFLE_SEEDS != None:
                    for seed in args.SHUFFLE_SEEDS:
                        logging.info("## {} SHUFFLED-{}, SEED: {}".format(design, dataset.upper(), seed))
                        print("## {} SHUFFLED-{}, SEED: {}".format(design, dataset.upper(), seed))

                        # set random seed for repeatabiblity
                        np.random.seed(seed)
                        # shuffle dataset
                        dict_diff_SHUFFLED[method] = dict_diff[method].iloc[np.random.permutation(np.arange(len(dict_diff[method])))]

                        for top in search_params_dict['tops']:
                            df_feature = df_exp.loc[df_exp.index.isin(dict_diff_SHUFFLED[method]["ID"][:top])]

                            df_exp_train = df_feature[df_design_train["sample"]]
                            df_label_train = df_design_train[args.SUBGROUP]

                            df_exp_test = df_feature[df_design_test["sample"]]
                            df_label_test = df_design_test[args.SUBGROUP]

                            if args.USE_LR:
                                dict_pred_tmp = uo.run_LR(df_exp_train, df_label_train, df_exp_test, df_label_test, design, method, dataset, top)
                                for key in dict_pred_tmp.keys():
                                    dict_pred_SHUFFLED[key] = dict_pred_SHUFFLED[key] + dict_pred_tmp[key]

                            if args.USE_RANDF:
                                dict_randF_tmp = uo.run_rand_forest(df_exp_train, df_label_train, df_exp_test, df_label_test, design, method, dataset, top)
                                for key in dict_randF_tmp.keys():
                                   dict_pred_SHUFFLED[key] = dict_pred_SHUFFLED[key] + dict_randF_tmp[key]

                            # add a seed column
                            dict_pred_SHUFFLED['seed'] = dict_pred_SHUFFLED['seed'] + [seed]

            pred_time = datetime.datetime.now()
            logging.info('PREDICTION ELAPSED {} design={}: {} seconds \t'.format(method.upper(), design, (pred_time - checkpoint).seconds ))

        # tranform dict of results to df
        df_LR = pd.DataFrame(dict_pred)
        df_LR_path = os.path.join(args.OUTDIR, 'tables', 'LR_{}.csv'.format(start_time.isoformat()))
        df_LR.to_csv(df_LR_path, index=False, sep="\t")
        logging.info('FILEPATH[df_LR]: {}'.format(df_LR_path))

        if args.SHUFFLE_SEEDS != None:
            df_LR_SHUFFLED = pd.DataFrame(dict_pred_SHUFFLED)
            df_LR_SHUFFLED_path = os.path.join(args.OUTDIR, 'tables', 'LR_SHUFFLED_{}.csv'.format(start_time.isoformat()))
            df_LR_SHUFFLED.to_csv(df_LR_SHUFFLED_path, index=False, sep="\t")
            logging.info('FILEPATH[df_LR_SHUFFLED]: {}'.format(df_LR_SHUFFLED_path))
            print('FILEPATH[df_LR_SHUFFLED]: {}'.format(df_LR_SHUFFLED_path))

    ##############################################################################################################
    ###################### COMPUTE PERF AND PLOT ##################################################################
    ##############################################################################################################

    #init contengency table for mcc computing
    ct = [0,0,0,0]
    for design in search_params_dict['designs']:
        for method in search_params_dict['methods']:
            for top in search_params_dict['tops']:
                for learn in df_LR.learn.unique() :
                    df = df_LR.loc[df_LR.design == design]
                    df = df.loc[df.top == top]
                    df = df.loc[df.method == method]
                    df = df.loc[df.learn == learn]

                    # continue only if df is non-empty
                    if len(df) < 2 : break

                    df = df.sort_values(by=['dataset'], ascending=[True])

                    # convert probs to array
                    probas = np.array([y for x in df.proba.values for y in x[:,0]])
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

                    dict_res['mcc'].append(mcc)
                    dict_res['method'].append(method)
                    dict_res['top'].append(top)
                    dict_res['design'].append(design)
                    dict_res['learn'].append(learn)

    # convert dict to df for analysis
    df_res = pd.DataFrame(dict_res)
    csv_out = os.path.join(args.OUTDIR, "tables", "res_{}.csv".format(start_time.isoformat()))
    df_res.to_csv(csv_out, index=False, sep="\t")

    # plot results
    plt_fig = sns.catplot(data=df_res, x="top", y="mcc", hue="method", col="design", row = "learn", kind="point", facet_kws=dict(subplot_kws=dict(ylim=[-0.05,1.05])))
    fig_out =  os.path.join(args.OUTDIR, "figures", "mcc_{}.pdf".format(start_time.isoformat()))
    plt_fig.savefig(fig_out)
    plt.close()

    end_time = datetime.datetime.now()
    logging.info('###### END RUN - TIME: {}\tRUNTIME: {} \t#####\n\n'.format(end_time.strftime("%Y-%m-%d %H:%M"),(end_time - start_time).seconds))
