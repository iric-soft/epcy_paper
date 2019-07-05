import os

from collections import defaultdict

import seaborn as sns
import pandas as pd
import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'DejaVu Sans'
mpl.rcParams['pdf.fonttype'] = 42

from .. utils import other as uo

def main_clust_exp(args, argparser):
    path_file_out = os.path.join(args.PATH, "res", "deg")
    if not os.path.exists(path_file_out):
        os.makedirs(path_file_out)

    file_dict = {
        'epcy' : "prediction_capability.xls",
        'deseq2' : "deseq2_genes.xls",
        'edger' : "edger_genes.xls",
        'limma' : "limma_voom_genes.xls",
    }

    # create dict with search params
    search_params_dict = {
        'tops' : args.TOP_VALUES,
        'methods' : args.METHODS,
        'designs' : args.DESIGN
    }

    prefix_file_out += "_" + args.PREFIX
    prefix_file_out += "_lfc" + str(args.LOG_FC)
    prefix_file_out += "_MCC" + str(args.MCC)

    designs = args.design.split(",")

    dict_perf = None
    dict_silh = None
    dict_diff = defaultdict(list)
    for design in search_params_dict['designs']:
        print(design)
        path_design_file_out = os.path.join(args.PATH, design, "res", "deg")

        if not os.path.exists(path_design_file_out):
            os.makedirs(path_design_file_out)

        path_out = os.path.join(path_design_file_out, prefix_file_out)

        dir_design = os.path.join(args.PATH, "design", design)
        df_design = uo.get_design(args, dir_design)
        df_exp = uo.get_exp(args, args.MATRIX)
        max_exp = np.amax(df_exp.values)

        df_exp = df_exp[df_design["sample"]]
        list_genes = df_exp.index.tolist()

        dict_dist = None
        for method in search_params_dict['methods']:
        #for top in search_params_dict['tops']:
            print(method)
            dict_diff[method] = uo.read_diff_table(args, file_dict[method], method, list_genes, path_dir=path_dataset_train)

            #for method in search_params_dict['methods']:
            for top in search_params_dict['tops']:
                print(top)
                silh_score = uo.get_exp_cluster(args, df_exp, df_design, dict_diff[method], path_out, method, design, top, max_exp)
