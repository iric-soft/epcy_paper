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
    # create path for log file
    logfile = os.path.join(log_dir, logname)
    # set config for logging (all info will go to logfile)
    if args.EVENT_LOG : logging.basicConfig(filename=logfile,level=logging.DEBUG)
    # log start time of run of module
    logging.info('###### BEGIN RUN - TIME: {} \t#####'.format(start_time.strftime("%Y-%m-%d %H:%M")))
    # log arguments
    logging.info(str(args))

    # if analysis is done on genes
    if args.GENE:
        file_dict = {
            # results filenames
            'kt' : "genes_prediction.xls",
            'kt_log2FC' : "genes_prediction.xls",
            'deseq2' :  "deseq2_genes_h5.xls",
            'deseq2_lrt' :"deseq2_genes_h5.xls",
            'edger_lrt' :  "edger_lrt_genes_h5.xls",
            'edger' : "edger_qlf_genes_h5.xls",
            'limma' : "limma_voom_genes_h5.xls",
        }

        col_ids_name = "gene_id"
        #file_exp = "/u/eaudemard/project/kt_paper/data/all_exp_404.csv"
        file_exp = "/u/eaudemard/project/kt_paper/data/all_exp_tcga_brca.csv"

    # if analysis is done on transcripts
    else:
        file_dict = {
            'kt' : "trans_prediction.xls",
            'kt_log2FC' : "trans_prediction.xls",
            'deseq2' : "deseq2_trans_h5.xls",
            'deseq2_lrt' : "deseq2_trans_h5.xls",
            'edger_lrt' : "edger_lrt_trans_h5.xls",
            'edger' : "edger_qlf_trans_h5.xls",
            'limma' : "limma_voom_trans_h5.xls",
        }

        col_ids_name = "trans_id"
        file_exp = "/u/eaudemard/project/kt_paper/data/"

    # set fold argument depending of input param for fold analysis
    fold = "fold" + str(args.FOLD) if not args.LOO else 'foldloo'
    # create dict with search params
    search_params_dict = {
        'top_values' : args.TOP_VALUES,
        'top_range' : args.TOP_RANGE, # top R best ranked genes to use for prediction
        'methods' : ["deseq2", "edger", "limma", "kt"],#, 'kt_log2FC'], ### tested methods, generic names, used for legending
        'subgroups' : args.SUBGROUPS ### hard code, but we can add argument
    }

    ## get all df for all samples & filters on biotype (ex coding)
    df_exp = uo.get_exp(args, file_exp)
    ## convert to list
    list_genes = df_exp.index.tolist()

    # log warmup time
    warmup = datetime.datetime.now()
    # log search params
    logging.info('Search parameters: {}'.format(search_params_dict))
    logging.info('WARM UP ELAPSED TIME : {} seconds\t'.format((warmup - start_time).seconds))
    # set tops
    tops = np.arange(args.TOP_RANGE[0], args.TOP_RANGE[1], args.TOP_RANGE[2]) if args.TOP_RANGE != None else args.TOP_VALUES
    # init a default dict of gene rankings
    dict_diff = defaultdict(list)
    # init empty dict of predictor result
    dict_logR = defaultdict(list)
    # init default dict of aucs (results)
    dict_auc = defaultdict(list)

    # init a default dict of gene rankings
    dict_diff_SHUFFLED = defaultdict(list)
    # init empty dict of predictor result SHUFFLED
    dict_logR_SHUFFLED = defaultdict(list)
    # init default dict of aucs (results)
    dict_auc_SHUFFLED = defaultdict(list)

    # cycle through subgroups
    for subgroup in search_params_dict['subgroups']:
        # updates path wrt subgroup and cv
        dir_dataset = os.path.join(args.PATH, subgroup, "cv", fold)
        # list all folds
        datasets = os.listdir(dir_dataset)
        # set number of datasets to include DEBUG ONLY
        n_datasets =  int(args.N_DATASETS) if args.N_DATASETS != 'all' else len(datasets) # len(datasets)
        if args.N_DATASETS != 'all': np.random.shuffle(datasets)
        # failed attempts to compute predictions
        failed = 0
        # cycle through methods
        for method in search_params_dict['methods']:
            logging.info('TRAINING/TESTING subgroup: {}, method: {}, n_datasets: {} / {}'.format(subgroup, method, n_datasets, len(datasets)))
            print('TRAINING/TESTING subgroup: {}, method: {}, n_datasets: {} / {}'.format(subgroup, method, n_datasets, len(datasets)))
            # log checkpoint
            checkpoint =  datetime.datetime.now()
            # cycle through datasets
            for dataset in datasets[:n_datasets]: # cycle through folds n = 404
                # cycle through tops

                path_dataset_train = os.path.join(dir_dataset, dataset, "train") # get train dataset path
                df_design_train = uo.get_design(args, path_dataset_train) # design of train data (project.csv) converted to df

                # same for test
                path_dataset_test = os.path.join(dir_dataset, dataset, "test")
                df_design_test = uo.get_design(args, path_dataset_test)

                # set top to all genes !
                allGenes = 25000
                # import dataframes of gene rankings for each method
                if method == "kt_log2FC":
                    df_ext_fc = None
                    if args.EXTFC:
                        df_ext_fc = uo.read_diff_table(args, file_dict["kt"], "kt", list_genes, path_dir=os.path.join(args.PATH, subgroup, "c"))
                        df_ext_fc = df_ext_fc[["ID", "log2_FC"]]
                        df_ext_fc.columns = ["ID", "log2_FC_ext"]
                    dict_diff[method] = uo.read_diff_table(args, file_dict[method], method, list_genes, path_dir=path_dataset_train, df_ext=df_ext_fc)
                else:
                    dict_diff[method] = uo.read_diff_table(args, file_dict[method], method, list_genes, path_dir=path_dataset_train)
                # compute Nmax
                Nmax = len(dict_diff[method])

                # log Nmax
                logging.info("##{}, {}, {} ## NMAX: {} ".format(subgroup, dataset.upper(),method.upper(), Nmax))
                print ("##{}, {}, {} ## NMAX: {} ".format(subgroup, dataset.upper(), method.upper(), Nmax))
                # log top 10 ranking for dataset
                logging.info("## RANKING {}, dataset {}, top 10 : {}".format(method,dataset, np.array(dict_diff[method]["ID"])[:10]))
                # cycle topX through 1 to Nmax with increment =  top_increment
                for top in tops:
                    #try :
                    ##############################################################################################################
                    ###################### TRAIN PREDICTORS AND COMPUTE PREDICTIONS ON ALL TEST SETS #############################
                    ##############################################################################################################
                    # define our gene signature
                    df_signature = df_exp.loc[df_exp.index.isin(dict_diff[method]["ID"][:top])]
                    # define our (X TRAIN)
                    df_exp_train = df_signature[df_design_train["sample"]]
                    # set our (Y TRAIN)
                    df_label_train = df_design_train[args.GROUP]

                    # define (X TEST)
                    df_exp_test = df_signature[df_design_test["sample"]]
                    # define (Y TEST)
                    df_label_test = df_design_test[args.GROUP]

                    if args.USE_LOGR:
                        # predict with logistic regression (may be change name)
                        dict_logR_tmp = uo.run_logR(df_exp_train, df_label_train, df_exp_test, df_label_test, subgroup, method, dataset, top)

                        # cleanup
                        for key in dict_logR_tmp.keys():
                            # append new result to entry's list
                            dict_logR[key] = dict_logR[key] + dict_logR_tmp[key]

                    if args.USE_RANDF:
                        # predict with logistic regression (may be change name)
                        dict_randF_tmp = uo.run_rand_forest(df_exp_train, df_label_train, df_exp_test, df_label_test, subgroup, method, dataset, top)

                        # cleanup
                        for key in dict_randF_tmp.keys():
                            # append new result to entry's list
                           dict_logR[key] = dict_logR[key] + dict_randF_tmp[key]
                    # add a seed column
                    dict_logR['seed'] = dict_logR['seed'] + [float('NaN')]



                    # except Exception as e:
                    #   logging.info("##TOP {}, {}, {} :LOGR COMPUTATION ERROR ==> Exception: {} ".format(top, method, dataset, e))
                    #   failed += 1
                ####### END (NO SHUFFLE) TOP FOR LOOP    ##############

                ####### BEGIN SHUFFLED TOP FOR LOOP       #############
                # only if specified
                if args.SHUFFLE_SEEDS != None:
                    # cycle through seeds
                    for seed in args.SHUFFLE_SEEDS:
                        # verbose

                        logging.info("## {} SHUFFLED-{}, SEED: {}".format(subgroup, dataset.upper(), seed))
                        print ("## {} SHUFFLED-{}, SEED: {}".format(subgroup, dataset.upper(), seed))
                       # set random seed for repeatabiblity
                        np.random.seed(seed)
                        # shuffle dataset
                        dict_diff_SHUFFLED[method] = dict_diff[method].iloc[np.random.permutation(np.arange(len(dict_diff[method])))]
                        # cycle through tops as usual
                        for top in tops:
                            # define our gene signature
                            df_signature = df_exp.loc[df_exp.index.isin(dict_diff_SHUFFLED[method]["ID"][:top])]
                            # define our (X TRAIN)
                            df_exp_train = df_signature[df_design_train["sample"]]
                            # set our (Y TRAIN)
                            df_label_train = df_design_train[args.GROUP]

                            # define (X TEST)
                            df_exp_test = df_signature[df_design_test["sample"]]
                            # define (Y TEST)
                            df_label_test = df_design_test[args.GROUP]

                            if args.USE_LOGR:
                                # predict with logistic regression (may be change name)
                                dict_logR_tmp = uo.run_logR(df_exp_train, df_label_train, df_exp_test, df_label_test, subgroup, method, dataset, top)

                                # cleanup
                                for key in dict_logR_tmp.keys():
                                    # append new result to entry's list
                                    dict_logR_SHUFFLED[key] = dict_logR_SHUFFLED[key] + dict_logR_tmp[key]

                            if args.USE_RANDF:
                                # predict with logistic regression (may be change name)
                                dict_randF_tmp = uo.run_rand_forest(df_exp_train, df_label_train, df_exp_test, df_label_test, subgroup, method, dataset, top)

                                # cleanup
                                for key in dict_randF_tmp.keys():
                                    # append new result to entry's list
                                   dict_logR_SHUFFLED[key] = dict_logR_SHUFFLED[key] + dict_randF_tmp[key]
                            # add a seed column
                            dict_logR_SHUFFLED['seed'] = dict_logR_SHUFFLED['seed'] + [seed]

                ####### END SHUFFLED DATASET TOP FORLOOP ##############

            ####### END DATASET FORLOOP ###########
            logging.info("GOT {} FAILED DATASET ANALYSES".format(failed))
            # get time
            pred_time = datetime.datetime.now()
            # log runtime for Prediction computing
            logging.info('PREDICTION ELAPSED {} subgroup={}: {} seconds \t'.format(method.upper(), subgroup,(pred_time - checkpoint).seconds ))

        # tranform dict of results to df
        df_logR = pd.DataFrame(dict_logR)
        # set path for output tables
        df_logR_path = os.path.join(args.OUTDIR, 'tables', 'logR_output_{}.csv'.format(start_time.isoformat()))
        # save logR df
        df_logR.to_csv(df_logR_path, index=False, sep="\t")
        # set path for table
        logging.info('FILEPATH[df_logR]: {}'.format(df_logR_path))

        if args.SHUFFLE_SEEDS != None:
            # tranform dict of results to df
            df_logR_SHUFFLED = pd.DataFrame(dict_logR_SHUFFLED)
            # set path for output tables
            df_logR_SHUFFLED_path = os.path.join(args.OUTDIR, 'tables', 'logR_output_SHUFFLED_{}.csv'.format(start_time.isoformat()))
            # save logR df
            df_logR_SHUFFLED.to_csv(df_logR_SHUFFLED_path, index=False, sep="\t")
            # log path to table
            logging.info('FILEPATH[df_logR_SHUFFLED]: {}'.format(df_logR_SHUFFLED_path))

        ####### END METHOD FORLOOP ###############
    ####### END SUBGROUP FORLOOP ############
    print('FILEPATH[df_logR_SHUFFLED]: {}'.format(df_logR_SHUFFLED_path))
    pdb.set_trace()
    ##############################################################################################################
    ###################### COMPUTE AUC AND PLOT ##################################################################
    ##############################################################################################################
    #init contengency table for mcc computing
    ct = [0,0,0,0]
    for subgroup in search_params_dict['subgroups']:
        for method in search_params_dict['methods']:
            # cycle thru tops
            for top in tops:
                for learn in ["L-Reg", "R-Forest"]:
                    # filter by top
                    df = df_logR.loc[df_logR.subgroup == subgroup]
                    # filter by top
                    df = df.loc[df.top == top]
                    # filter on method
                    df = df.loc[df.method == method]
                    # filter on learn
                    df = df.loc[df.learn == learn]
                    # continue only if df is non-empty
                    if len(df) < 2 : break
                    ####### sort ??
                    df = df.sort_values(by=['dataset'], ascending=[True])
                    # convert probs to array
                    # labels = np.array([label.strip('[').strip(']') for label in df.label]).astype(np.int)
                    probas = np.array([y for x in df.proba.values for y in x[:,0]])
                    # same for labels
                    labels = np.array([label for x in df.label.values for label in x])
                    # get query ids array
                    ids_query = np.where(labels == 1)
                    # get ref ids array
                    ids_ref = np.where(labels == 0)

                    # compute mcc
                    ct[0] = np.sum(probas[ids_query] > 0.5)
                    ct[1] = np.sum(probas[ids_query] <= 0.5)
                    ct[3] = np.sum(probas[ids_ref] <= 0.5)
                    ct[2] = np.sum(probas[ids_ref] > 0.5)
                    mcc = uo.get_mcc(ct)

                    acc = (ct[0]+ct[3]) / (len(df.label))
                    # compute auc for method's proba.
                    auc = uo.auc_u_test(probas, ids_query, ids_ref)
                    #
                    dict_auc['auc'].append(auc[0])
                    dict_auc['mcc'].append(mcc)
                    dict_auc['acc'].append(acc)
                    dict_auc['method'].append(method)
                    dict_auc['top'].append(top)
                    dict_auc['subgroup'].append(subgroup)
                    dict_auc['learn'].append(learn)

    # convert dict to df for analysis
    df_auc = pd.DataFrame(dict_auc)
    # create path to csv file
    csv_out = os.path.join(args.OUTDIR, "tables", "auc_output_{}.csv".format(start_time.isoformat()))
    # save csv to path
    df_auc.to_csv(csv_out, index=False, sep="\t")

    # plot results
    plt_auc = sns.catplot(data=df_auc, x="top", y="auc", hue="method", col="subgroup", row = "learn", kind="point")#, facet_kws=dict(subplot_kws=dict(alpha=.1)))
    # create path to figure output
    fig_out =  os.path.join(args.OUTDIR, "figures", "auc_output_{}.pdf".format(start_time.isoformat()))
    # save
    plt_auc.savefig(fig_out)
    # close plot
    plt.close()

    # plot results
    plt_auc = sns.catplot(data=df_auc, x="top", y="mcc", hue="method", col="subgroup", row = "learn", kind="point", facet_kws=dict(subplot_kws=dict(ylim=[-0.05,1.05])))
    # create path to figure output
    fig_out =  os.path.join(args.OUTDIR, "figures", "mcc_output_{}.pdf".format(start_time.isoformat()))
    # save
    plt_auc.savefig(fig_out)
    # close plot
    plt.close()

    # get time and date
    end_time = datetime.datetime.now()
    # log runtime
    logging.info('###### END RUN - TIME: {}\tRUNTIME: {} \t#####\n\n'.format(end_time.strftime("%Y-%m-%d %H:%M"),(end_time - start_time).seconds))
