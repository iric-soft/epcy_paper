import datetime
# get time and date
start_time = datetime.datetime.now()
# make rest of imports
import os

from collections import defaultdict
import pdb
import seaborn as sns
import pandas as pd
import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt

from .. utils import other as uo

import logging

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'DejaVu Sans'
mpl.rcParams['pdf.fonttype'] = 42

def main_train_pred(args, argparser):

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

    # if analysis is done on genes
    if args.GENE:
        file_dict = {
            # results filenames
            'kt' : "genes_prediction.xls",
            'deseq2' :  "deseq2_genes_h5.xls",
            'deseq2_lrt' :"deseq2_genes_h5.xls",
            'edger_lrt' :  "edger_lrt_genes_h5.xls",
            'edger_qlf' : "edger_qlf_genes_h5.xls",
            'limma' : "limma_voom_genes_h5.xls",
        }

        col_ids_name = "gene_id"
        file_exp = "/u/eaudemard/project/kt_paper/data/all_exp_404.csv"

    # if analysis is done on transcripts
    else:
        file_dict = {
            'kt' : "trans_prediction.xls",
            'deseq2' : "deseq2_trans_h5.xls",
            'deseq2_lrt' : "deseq2_trans_h5.xls",
            'edger_lrt' : "edger_lrt_trans_h5.xls",
            'edger_qlf' : "edger_qlf_trans_h5.xls",
            'limma' : "limma_voom_trans_h5.xls",
        }

        col_ids_name = "trans_id"
        file_exp = "/u/eaudemard/project/kt_paper/data/"

    search_params_dict = {
        'tops' : [1, 3, 10, 20, 50, 100, 200], # top R best ranked genes to use for prediction
        'methods' : ["deseq2", "edger", "limma", "kt"], ### tested methods, generic names, used for legending
        'subgroups' : ["28_inv16_vs28", "28_inv16", "18_t8_21", "132_FLT3-ITD", "139_NPM1_mut"] ### hard code, but we can add argument
    }

    ## get all df for all samples & filters on biotype (ex coding)
    df_exp = uo.get_exp(args, file_exp)
    df_exp_pred = uo.get_exp(args, args.EXP)

    df_exp = df_exp[pd.notnull(df_exp)]
    df_exp_pred = df_exp_pred[pd.notnull(df_exp_pred)]

    df_exp.drop_duplicates(inplace=True)
    df_exp_pred.drop_duplicates(inplace=True)

    df_exp = pd.merge(df_exp, df_exp_pred, left_index=True, right_index=True)
    list_genes = df_exp.index.tolist()

    ##############################################################################################################
    ###################### TRAIN PREDICTORS AND COMPUTE PREDICTIONS ON ALL TEST SETS #############################
    ##############################################################################################################
    dict_supervised = defaultdict(list)
    dict_unsupervised = defaultdict(list)
    for subgroup in search_params_dict['subgroups']:
        print("##########" + subgroup + "##########")
        dir_dataset = os.path.join(args.PATH, subgroup)
        df_design_train = uo.get_design(args, path=dir_dataset) # design of train data (project.csv) converted to df

        df_design_test = uo.get_design(args, path=dir_dataset, file_name="project_aml_tcga.csv")

        dict_diff = defaultdict(list)
        for top in search_params_dict['tops']: # n = 5
            print(top)
            for method in search_params_dict['methods']: # n = 4
######          # add kt_binom_mcc , kt_auc, kt_kde_mcc (default)

                dict_diff[method] = uo.read_diff_table(args, file_dict[method], method, list_genes, dir_dataset)

                # define our gene signature
                df_signature = pd.DataFrame(df_exp.loc[df_exp.index.isin(dict_diff[method]["ID"][:top])])
                # define our (X TRAIN)
                df_exp_train = df_signature[df_design_train["sample"]]
                # set our (Y TRAIN)
                df_label_train = df_design_train[args.GROUP]

                # define (X TEST)
                df_exp_test = df_signature[df_design_test["sample"]]
                # define (Y TEST)
                df_label_test = df_design_test[args.GROUP]
                # predict with logistic regression (may be change name)
                dict_tmp = uo.run_logR(df_exp_train, df_label_train, df_exp_test, df_label_test, subgroup, method, "dataset_0", top)

                for key in dict_tmp.keys():
                    dict_supervised[key] = dict_supervised[key] + dict_tmp[key]

                dict_tmp = uo.run_rand_forest(df_exp_train, df_label_train, df_exp_test, df_label_test, subgroup, method, "dataset_0", top)

                for key in dict_tmp.keys():
                    dict_supervised[key] = dict_supervised[key] + dict_tmp[key]

                dict_tmp = uo.run_kmeans(df_exp_test, df_label_test, subgroup, method, "dataset_0", top)

                for key in dict_tmp.keys():
                    dict_unsupervised[key] = dict_unsupervised[key] + dict_tmp[key]


    df_supervised = pd.DataFrame(dict_supervised)
    df_unsupervised = pd.DataFrame(dict_unsupervised)
    ##############################################################################################################
    ############################## COMPUTE PREDICTION SCORES FOR TOPS, METHODS  ##################################
    ##############################################################################################################
    ct = [0,0,0,0]
    dict_auc = defaultdict(list)

    for subgroup in search_params_dict['subgroups']:
        for top in search_params_dict['tops']:
            for method in search_params_dict['methods']:
                for learn in ['L-Reg', 'R-Forest']:
                    df = df_supervised.loc[df_supervised.top == top]
                    df = df.loc[df.method == method]
                    df = df.loc[df.learn == learn]
                    df = df.loc[df.subgroup == subgroup]
                    # convert probs to array
                    proba = [y for x in df.proba.values for y in x[:,0] ]
                    # same for labels
                    label = [label for x in df.label.values for label in x]
                    label = np.array(label)
                    proba = np.array(proba)
                    ids_query = np.where(label == 1)
                    ids_ref = np.where(label == 0)
                    # compute mcc
                    ct[0] = np.sum(proba[ids_query] > 0.5)
                    ct[1] = np.sum(proba[ids_query] <= 0.5)
                    ct[3] = np.sum(proba[ids_ref] < 0.5)
                    ct[2] = np.sum(proba[ids_ref] >= 0.5)
                    mcc = uo.get_mcc(ct)
                    # compute auc for method's proba.
                    auc = uo.auc_u_test(proba, ids_query, ids_ref)
                    dict_auc['auc'].append(auc[0])
                    dict_auc['mcc'].append(mcc)
                    dict_auc['method'].append(method)
                    dict_auc['top'].append(top)
                    dict_auc['subgroup'].append(subgroup)
                    dict_auc['learn'].append(learn)
                    dict_auc['type'].append("test")

                    proba = [y for x in df.proba_train.values for y in x[:,0] ]
                    # same for labels
                    label = [label for x in df.label_train.values for label in x]
                    label = np.array(label)
                    proba = np.array(proba)
                    ids_query = np.where(label == 1)
                    ids_ref = np.where(label == 0)
                    # compute mcc
                    ct[0] = np.sum(proba[ids_query] > 0.5)
                    ct[1] = np.sum(proba[ids_query] <= 0.5)
                    ct[3] = np.sum(proba[ids_ref] < 0.5)
                    ct[2] = np.sum(proba[ids_ref] >= 0.5)
                    mcc = uo.get_mcc(ct)
                    #print(proba)
                    #print(label)
                    # compute auc for method's proba.
                    auc = uo.auc_u_test(proba, ids_query, ids_ref)
                    dict_auc['auc'].append(auc[0])
                    dict_auc['mcc'].append(mcc)
                    dict_auc['method'].append(method)
                    dict_auc['top'].append(top)
                    dict_auc['subgroup'].append(subgroup)
                    dict_auc['learn'].append(learn)
                    dict_auc['type'].append("train")


    df_auc = pd.DataFrame(dict_auc)
    df_auc_test =  df_auc.loc[df_auc.type == "test"]
    plt_auc = sns.catplot(data=df_auc_test, x="top", y="mcc", hue="method", col="subgroup", row="learn", kind="point")#, facet_kws=dict(linestyles=["-", "--"])
    fig_out = os.path.join(fig_dir, "mcc_tcga_aml_test.pdf")
    plt_auc.savefig(fig_out)
    plt.close()

    df_auc = pd.DataFrame(dict_auc)
    df_auc_test =  df_auc.loc[df_auc.type == "test"]
    plt_auc = sns.catplot(data=df_auc_test, x="top", y="auc", hue="method", col="subgroup", row="learn", kind="point")#, facet_kws=dict(linestyles=["-", "--"])
    fig_out = os.path.join(fig_dir, "auc_tcga_aml_test.pdf")
    plt_auc.savefig(fig_out)
    plt.close()

    plt_pred = sns.catplot(data=df_unsupervised, x="top", y="pred", hue="method", col="subgroup", row="learn", kind="point")#, facet_kws=dict(linestyles=["-", "--"])
    fig_out = os.path.join(fig_dir, "pred_tcga_aml_test.pdf")
    plt_pred.savefig(fig_out)
    plt.close()

    df_auc_train = df_auc.loc[df_auc.type == "train"]
    plt_auc = sns.catplot(data=df_auc_train, x="top", y="auc", hue="method", col="subgroup", row="learn", kind="point")#, facet_kws=dict(linestyles=["-", "--"])
    fig_out = os.path.join(fig_dir, "auc_tcga_aml_train.pdf")
    plt_auc.savefig(fig_out)
    plt.close()

    df_auc_train = df_auc.loc[df_auc.type == "train"]
    plt_auc = sns.catplot(data=df_auc_train, x="top", y="mcc", hue="method", col="subgroup", row="learn", kind="point")#, facet_kws=dict(linestyles=["-", "--"])
    fig_out = os.path.join(fig_dir, "mcc_tcga_aml_train.pdf".format(start_time.isoformat()))
    plt_auc.savefig(fig_out)
    plt.close()

    df_auc_path = os.path.join(csv_dir, "tcga_aml.xls")
    df_auc.to_csv(df_auc_path, index=False, sep="\t")
