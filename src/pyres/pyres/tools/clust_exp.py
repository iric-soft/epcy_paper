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
        'num_cluster' : args.NUM_CLUSTER,
        'tops' : args.TOP_VALUES,
        'methods' : args.METHODS,
        'designs' : args.DESIGN
    }

    dict_diff = defaultdict(list)
    dict_res = defaultdict(list)
    for design in search_params_dict['designs']:
        path_out = os.path.join(args.OUTDIR, "clust_exp", "kmean", design)

        if not os.path.exists(path_out):
            os.makedirs(path_out)


        dir_design = os.path.join(args.PATH, design)
        df_design = uo.get_design(args, dir_design)
        df_exp = uo.get_exp(args, args.MATRIX)
        max_exp = np.amax(df_exp.values)

        df_exp = df_exp[df_design["sample"]]
        list_genes = df_exp.index.tolist()

        dict_dist = None
        for method in search_params_dict['methods']:
        #for top in search_params_dict['tops']:
            dict_diff[method] = uo.read_diff_table(args, file_dict[method], method, dir_design, args.PVALUE, list_genes)

            #for method in search_params_dict['methods']:
            for top in search_params_dict['tops']:
                for num_c in search_params_dict['num_cluster']:
                    print("##{}, {}, TOP {}, CLUSETER {}".format(design, method.upper(), str(top), str(num_c)))
                    silh_score = uo.get_exp_cluster(args, df_exp, df_design, dict_diff[method], path_out, method, design, top, num_c, max_exp)

                    dict_res['silh_score'].append(silh_score)
                    dict_res['method'].append(method)
                    dict_res['top'].append(top)
                    dict_res['design'].append(design)
                    dict_res['num_cluster'].append(num_c)

    # create subfolder for log and res
    fig_dir = os.path.join(args.OUTDIR, "clust_exp", "kmean")
    if not os.path.exists(fig_dir):
        os.makedirs(fig_dir)

    # convert dict to df for analysis
    # df_res = pd.DataFrame(dict_res)
    # csv_out = os.path.join(fig_dir, "clust_table.csv")
    # df_res.to_csv(csv_out, index=False, sep="\t")

    # plot results
    # plt_fig = sns.catplot(data=df_res, x="top", y="silh_score", hue="method", row="num_cluster", col="design", kind="point", facet_kws=dict(subplot_kws=dict(ylim=[-0.05,1.05])))
    # fig_out =  os.path.join(fig_dir, "clust_silh.pdf")
    # plt_fig.savefig(fig_out)
    # plt.close()
