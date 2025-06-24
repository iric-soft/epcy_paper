

# get time and date
# make rest of imports
import os
import sys
import gc

from collections import defaultdict
from collections import Counter

import pdb

import pandas as pd
import numpy as np
from scipy import spatial

import matplotlib as mpl
# set no display
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import LogNorm

import seaborn as sns

#from utils import other as uo
from .. utils import other as uo
import random

import datetime
from scipy.stats import pearsonr
from scipy.stats import linregress

start_time = datetime.datetime.now()

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'DejaVu Sans'
mpl.rcParams['pdf.fonttype'] = 42

# for test
# from types import SimpleNamespace
# args = SimpleNamespace(**{ 'N_GENES': 21953, 'MCC': -2, 'LOG_FC':0, 'TYPE_QUANT': 'readcounts', 'QUANT': 'cellranger', 'cellNumber':[3000, 5000, 8000, 10000], 'E_FPR': [0.00001, 0.0001, 0.001], 'RANDOMS': [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20], 'REPS': [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20], 'CELLTYPES': ['cd14', 'cd56_nk', 'cd34', 'cytotoxic_t', 'b_cells', 'naive_t', 'memory_t', 'regulatory_t', 'cd4', 'naive_cytotoxic'], 'METHODS': ['epcy', 'trend', 'mast', 'wilcox'], 'PATH': '/u/eaudemard/project/epcy_paper/data/design', 'OUTDIR': '/u/eaudemard/project/epcy_paper/data/res', 'SUBGROUP': 'subgroup' })


def main_eval_sc(args, argparser):
    df_id2name = pd.read_csv("../../../data/other/Homo_sapiens.GRCh37.82.id2name.tsv", sep="\t")

    file_dict = {
        'epcy': "predictive_capability.tsv",
        'epcy_bagging': "predictive_capability.tsv",
        'deseq2': "deseq2_genes.xls",
        'edger': "edger_genes.xls",
        'voom': "limma_voom_genes.xls",
        'trend': "limma_trend_genes.xls",
        'mast': "mast_genes.xls",
        'wilcox': "presto_genes.xls",
    }

    # create dict with search params
    search_params_dict = {
        'cellNumber': args.cellNumber,
        'reps': args.REPS,
        'randoms': args.RANDOMS,
        'e_fpr': args.E_FPR,
        'cellTypes': args.CELLTYPES,
        'methods': args.METHODS
    }

    col_pal = sns.color_palette("colorblind")
    col_pal = [col_pal[i] for i in [9, 2, 4, 7]]
    col_pal_cmap = sns.color_palette("colorblind", as_cmap=True)
    col_pal_cmap = [col_pal_cmap[i] for i in [9, 2, 4, 7]]
    col_pal_hex = sns.color_palette("colorblind").as_hex()
    col_pal_hex = [col_pal_hex[i] for i in [9, 2, 4]]

    fig_dir = os.path.join(args.OUTDIR, "eval_sc")
    if not os.path.exists(fig_dir):
        os.makedirs(fig_dir)

    log_c = sys.float_info.min
    dict_analyses = defaultdict(list)
    for method in search_params_dict['methods']:
        dict_diff = defaultdict(list)
        value_adj_colname = ""
        value_colname = ""
        l2fc_colname = ""
        if "epcy" in method:
            value_adj_colname = "KERNEL_MCC"
            value_colname = "KERNEL_MCC"
            l2fc_colname = "L2FC"
        elif method == "deseq2":
            value_adj_colname = "padj"
            value_colname = ""
            l2fc_colname = "log2FoldChange"
        elif method == "edger":
            value_adj_colname = "FDR"
            value_colname = "PValue"
            l2fc_colname = "logFC"
        elif method == "voom":
            value_adj_colname = "adj.P.Val"
            value_colname = "P.Value"
            l2fc_colname = "logFC"
        elif method == "trend":
            value_adj_colname = "adj.P.Val"
            value_colname = "P.Value"
            l2fc_colname = "logFC"
        elif method == "mast":
            value_adj_colname = "adj.P.Val"
            value_colname = "pval"
            l2fc_colname = "logFC"
        elif method == "wilcox":
            value_adj_colname = "padj"
            value_colname = "pval"
            l2fc_colname = "logFC"

        for cellNumber in search_params_dict['cellNumber']:
            for cellType in search_params_dict['cellTypes']:
                for rep in search_params_dict['reps']:
                    print('SC SUBSAMPLING: cellNumber: {}, cellType: {}, rep: {}, method: {}'.format(str(cellNumber), cellType, rep, method))
                    path_sc = os.path.join(args.PATH, '10X_FACS_reduce_{}_{}'.format(str(cellNumber), str(rep)))
                    design = uo.find_folder(path_sc, cellType)
                    path_design = os.path.join(path_sc, design)
                    df_design = uo.get_design(args, "Query", path_design)
                    num_ref = df_design[df_design[args.SUBGROUP] == 0].count().iloc[0]
                    num_query = df_design[df_design[args.SUBGROUP] == 1].count().iloc[0]
                    df_diff = uo.read_diff_table(
                        args, file_dict[method], method,
                        path_design, 1
                    )

                    if "epcy" in method:
                        dict_analyses_tmp = defaultdict(list)

                        dict_analyses_tmp['cellNumber'] = [cellNumber] 
                        dict_analyses_tmp['cellType'] = [cellType] 
                        dict_analyses_tmp['rep'] = [rep] 

                        dict_analyses_tmp['num_query'] = [num_query] 
                        dict_analyses_tmp['num_ref'] = [num_ref] 

                        for key in dict_analyses_tmp.keys():
                            dict_analyses[key] = dict_analyses[key] + dict_analyses_tmp[key]
                    else:
                        df_diff[value_colname] = -np.log10(df_diff[value_colname] + log_c)
                        df_diff[value_adj_colname] = -np.log10(df_diff[value_adj_colname] + log_c)

                    dict_diff_tmp = defaultdict(list)
                    dict_diff_tmp['ID'] = df_diff["ID"].tolist()
                    dict_diff_tmp['cellNumber'] = [cellNumber] * df_diff["ID"].shape[0]
                    dict_diff_tmp['cellType'] = [cellType] * df_diff["ID"].shape[0]
                    dict_diff_tmp['rep'] = [rep] * df_diff["ID"].shape[0]
                    dict_diff_tmp['method'] = [method] * df_diff["ID"].shape[0]
                    dict_diff_tmp['value'] = (df_diff[value_adj_colname]).tolist()
                    dict_diff_tmp['l2fc'] = df_diff[l2fc_colname].tolist()

                    for key in dict_diff_tmp.keys():
                        dict_diff[key] = dict_diff[key] + dict_diff_tmp[key]
        df_diff = pd.DataFrame(dict_diff)  
        del dict_diff
        fig_dir = os.path.join(args.OUTDIR, "eval_sc")
        csv_out = os.path.join(fig_dir, method + "_diff_table.csv")
        df_diff.to_csv(csv_out, index=False, sep="\t")

    list_df_diff = []
    for method in search_params_dict['methods']:
        fig_dir = os.path.join(args.OUTDIR, "eval_sc")
        csv_out = os.path.join(fig_dir, method + "_diff_table.csv")
        list_df_diff.append(pd.read_csv(csv_out, sep="\t"))
    
    df_diff = pd.concat(list_df_diff, axis=0, ignore_index=True)
    del list_df_diff
    
    print('SAVE TABLES')    
    fig_dir = os.path.join(args.OUTDIR, "eval_sc")
    csv_out = os.path.join(fig_dir, "diff_table.csv")
    df_diff.to_csv(csv_out, index=False, sep="\t")
    # df_diff = pd.read_csv(csv_out, sep="\t")
    df_diff.loc[df_diff["method"] == "trend", "l2fc"] = -df_diff.loc[df_diff["method"] == "trend", "l2fc"]
    
    df_analyses = pd.DataFrame(dict_analyses) 
    del dict_analyses 
    fig_dir = os.path.join(args.OUTDIR, "eval_sc")
    csv_out = os.path.join(fig_dir, "analyses_table.csv")
    df_analyses.to_csv(csv_out, index=False, sep="\t")
    # df_analyses = pd.read_csv(csv_out, sep="\t")

    df_analyses.sort_values(["cellNumber", "cellType", "rep"], inplace=True)
    df_diff.sort_values(by='value', ascending=False, inplace=True)
    
    ####### cuttoff on all cellType

    quantiles = [0.9999, 0.9995, 0.999, 0.995, 0.990]
    df_cutoffs = []
    for quantile in quantiles:
        index_cutoff = int(args.N_GENES * len(search_params_dict['cellTypes']) * (1-quantile)) * len(search_params_dict['reps'])

        df_cutoff = df_diff.groupby(["cellNumber","method"]).apply(lambda x: x.iloc[index_cutoff], include_groups=False).reset_index()
        df_cutoff["quantile"] = quantile
        df_cutoff["index_cutoff_rep"] = index_cutoff
        df_cutoff["index_cutoff"] = index_cutoff / len(search_params_dict['reps'])
        df_cutoffs.append(df_cutoff)
    
    df_cutoff = pd.concat(df_cutoffs, axis=0, ignore_index=True)
    del df_cutoffs

    df_diff_counts = []
    dict_uniq_id = defaultdict(list)
    for row in df_cutoff.iterrows():
        cellNumber = row[1]["cellNumber"]
        method = row[1]["method"]
        cutoff = row[1]["value"]
        quantile = row[1]["quantile"]
        
        df_diff_selected = df_diff[(df_diff["cellNumber"] == cellNumber) & (df_diff["method"] == method) & (df_diff["value"] >= cutoff)].copy()
        df_diff_count = pd.DataFrame(df_diff_selected.groupby(["cellType"])['ID'].value_counts()).reset_index()
        df_diff_count["cellNumber"] = cellNumber
        df_diff_count["method"] = method
        df_diff_count["quantile"] = quantile
        df_diff_counts.append(df_diff_count)

        dict_uniq_id_tmp = defaultdict(list)
        dict_uniq_id_tmp['cellNumber'] = [cellNumber] 
        dict_uniq_id_tmp['method'] = [method] 
        dict_uniq_id_tmp['style'] = ["At leat one replicat"]
        dict_uniq_id_tmp['quantile'] = [quantile] 
        dict_uniq_id_tmp['expected'] = [row[1]["index_cutoff"]] 
        dict_uniq_id_tmp['selected'] = [df_diff_count.shape[0]] 

        for key in dict_uniq_id_tmp.keys():
            dict_uniq_id[key] = dict_uniq_id[key] + dict_uniq_id_tmp[key]
        
        dict_uniq_id_tmp = defaultdict(list)
        dict_uniq_id_tmp['cellNumber'] = [cellNumber] 
        dict_uniq_id_tmp['method'] = [method]
        dict_uniq_id_tmp['style'] = ["All replicats"]
        dict_uniq_id_tmp['quantile'] = [quantile] 
        dict_uniq_id_tmp['expected'] = [row[1]["index_cutoff"]] 
        dict_uniq_id_tmp['selected'] = [df_diff_count[df_diff_count["count"] == 20].shape[0]] 
        
        for key in dict_uniq_id_tmp.keys():
            dict_uniq_id[key] = dict_uniq_id[key] + dict_uniq_id_tmp[key]

    df_diff_count = pd.concat(df_diff_counts, axis=0, ignore_index=True)
    df_uniq_id = pd.DataFrame(dict_uniq_id)

    plt_fig = sns.relplot(
        data=df_uniq_id, 
        x="cellNumber", y="selected", hue="method", style="style",
        markers=True, col="quantile", kind="line", palette=col_pal,
        col_wrap=6, col_order=quantiles, facet_kws=dict(sharey=False)
    )
    plt_fig.add_legend()
    sns.move_legend(plt_fig, "upper right", bbox_to_anchor=(1, 1))
    for ax in plt_fig.axes:
        quantile = float(ax.get_title().replace("quantile = ", ""))
        df_tmp = df_uniq_id[df_uniq_id['quantile'] == quantile].copy()
        ax.axhline(y=df_tmp["expected"].iloc[0], color='r', linestyle='-.')
        for idx,row in df_tmp.iterrows():
            x = row["cellNumber"]
            y = row["selected"]
            text = row["selected"]
            ax.text(x+.05,y,text, horizontalalignment='left')
    fig_dir = os.path.join(args.OUTDIR, "eval_sc")            
    fig_out = os.path.join(fig_dir, "num_unique_id_selected_zero_pvalue.pdf")
    plt_fig.savefig(fig_out)
    plt.close()
    
    plt_fig = sns.relplot(
        data=df_uniq_id, 
        x="cellNumber", y="selected", hue="method", style="style",
        markers=True, col="quantile", kind="line", palette=col_pal,
        col_wrap=6, col_order=quantiles, facet_kws=dict(sharey=False)
    )
    plt_fig.add_legend()
    #sns.move_legend(plt_fig, "upper left", bbox_to_anchor=(0, 0))
    for ax in plt_fig.axes:
        quantile = float(ax.get_title().replace("quantile = ", ""))
        ax.set_yscale("log")
        df_tmp = df_uniq_id[df_uniq_id['quantile'] == quantile].copy()
        ax.axhline(y=df_tmp["expected"].iloc[0], color='r', linestyle='-.')
        for idx,row in df_tmp.iterrows():
            x = row["cellNumber"]
            y = row["selected"]
            text = row["selected"]
            ax.text(x+.05,y,text, horizontalalignment='left')
    fig_dir = os.path.join(args.OUTDIR, "eval_sc")            
    fig_out = os.path.join(fig_dir, "num_unique_id_selected_zero_pvalue_log.pdf")
    plt_fig.savefig(fig_out)
    plt.close()
    
    for cellNumber in search_params_dict['cellNumber']:
        for quantile in quantiles:
            df_diff_count_tmp = df_diff_count[(df_diff_count["cellNumber"] == cellNumber) & (df_diff_count["quantile"] == quantile)].copy()

            if df_diff_count_tmp.shape[0] != 0:
                fig_dir = os.path.join(args.OUTDIR, "eval_sc", str(cellNumber))
                if not os.path.exists(fig_dir):
                    os.makedirs(fig_dir)
                
                plt_fig = sns.FacetGrid(df_diff_count_tmp, col="method", row="cellType")
                plt_fig.map_dataframe(
                    sns.histplot, 
                    x="count",bins=20, 
                )
                fig_out = os.path.join(fig_dir, str(quantile) + "_count_selected.pdf")
                plt_fig.savefig(fig_out)
                plt.close()

    df_detected_highly = df_diff_count[df_diff_count["count"] > 15].groupby(["method", "cellNumber", "cellType", "quantile"]).size().reset_index(name='num_gene')
    df_detected_highly["label"] = "highly replicated"
    df_detected_lowly = df_diff_count[df_diff_count["count"] <= 15].groupby(["method", "cellNumber", "cellType", "quantile"]).size().reset_index(name='num_gene')
    df_detected_lowly["label"] = "lowly replicated"

    df_dected = pd.concat([df_detected_highly, df_detected_lowly], axis=0, ignore_index=True)
    df_dected = df_dected.fillna(0)

    for cellNumber in search_params_dict['cellNumber']:
        df_dected_tmp = df_dected[(df_dected["cellNumber"] == cellNumber)].copy()
        fig_dir = os.path.join(args.OUTDIR, "eval_sc", str(cellNumber))
        if not os.path.exists(fig_dir):
            os.makedirs(fig_dir)
        
        plt_fig = sns.catplot(kind='bar', data=df_dected_tmp, x='cellType', y='num_gene', hue='label', col='method', row='quantile', sharey=False)
        fig_out = os.path.join(fig_dir, "count_selected.pdf")
        plt_fig.savefig(fig_out)
        plt.close()


###### volcano plot

    df_diff_mean = df_diff.groupby(["method", "cellNumber", "cellType", "ID"])["value"].mean()
    df_diff_mean = df_diff_mean.to_frame()
    df_diff_mean = df_diff_mean.rename({"value": "mean_value"}, axis='columns')
    df_diff_mean.reset_index(inplace=True)

    df_diff_std = df_diff.groupby(["method", "cellNumber", "cellType", "ID"])["value"].std(ddof=0)
    df_diff_std = df_diff_std.to_frame()
    df_diff_std = df_diff_std.rename({"value": "y_std"}, axis='columns')
    df_diff_std.reset_index(inplace=True)

    df_diff_l2fc = df_diff.groupby(["method", "cellNumber", "cellType", "ID"])["l2fc"].mean()
    df_diff_l2fc = df_diff_l2fc.to_frame()
    df_diff_l2fc = df_diff_l2fc.rename({"l2fc": "mean_l2fc"}, axis='columns')
    df_diff_l2fc.reset_index(inplace=True)

    df_diff_resume = pd.merge(df_diff_mean, df_diff_std, on=['method', 'cellNumber', "cellType", 'ID'])
    df_diff_resume = pd.merge(df_diff_resume, df_diff_l2fc, on=['method', 'cellNumber', "cellType", 'ID'])

    #df_selected = df_diff_resume[(df_diff_resume["method"] == "epcy") & (df_diff_resume["cellNumber"] == 10000) & (df_diff_resume["cellType"] == "cd34")].copy()
    #high_std_ids = df_selected[df_selected["y_std"] > 0.05]['ID'].values

    #df_diff_selected = df_diff[(df_diff["method"] == "epcy") & (df_diff["cellNumber"] == 10000) & (df_diff["cellType"] == "cd34")].copy()
    #df_diff_selected["std"] = "low"
    #df_diff_selected.loc[df_diff_selected['ID'].isin(high_std_ids), "std"] = "high"

    #path_sc = os.path.join(args.PATH, '10X_FACS_reduce_{}_{}'.format(str(10000), str(1)))
    #design = uo.find_folder(path_sc, "cd34")
    #path_design = os.path.join(path_sc, design)
    #df_epcy = uo.read_diff_table(
    #    args, file_dict["epcy"], "epcy",
    #    path_design, 1
    #)
    #df_epcy["std"] = "low"
    #df_epcy.loc[df_epcy['ID'].isin(high_std_ids), "std"] = "high"
    #df_epcy[df_epcy["std"] == "high"]


    for method in search_params_dict['methods']:
        df_diff_method = df_diff_resume[(df_diff_resume["method"] == method)].copy()
        df_diff_method.sort_values(by='y_std', ascending=True, inplace=True)

        if method == "wilcox":
            col_dot = [col_pal[3]]
        if method == "mast":
            col_dot = [col_pal[2]]
        if method == "trend":
            col_dot = [col_pal[1]]
        if method == "epcy":
            col_dot = [col_pal[0]]

        df_diff_method['y_std_binned'] = pd.cut(df_diff_method['y_std'], bins=4).astype(str)
        intervals = df_diff_method['y_std_binned'].unique()
        size_mapping = {intervals[0]: 20, intervals[1]: 80, intervals[2]: 140, intervals[3]: 200}

        plt_fig = sns.relplot(
            x="mean_l2fc", y="mean_value", hue="method",
            row="cellType", col="cellNumber", size="y_std_binned", sizes=size_mapping,
            col_order=[3000, 5000, 8000, 10000],
            palette=col_dot,
            data=df_diff_method
        )

        fig_out = os.path.join(fig_dir, method + "_volcano_4.pdf")
        plt_fig.savefig(fig_out)
        plt.close()

        fig_out = os.path.join(fig_dir, method + "_volcano_4.png")
        plt_fig.savefig(fig_out)
        plt.close()

    for method in search_params_dict['methods']:
        for cellType in search_params_dict['cellTypes']:
            df_diff_method = df_diff_resume[(df_diff_resume["method"] == method) &  (df_diff_resume["cellType"] == cellType) ].copy()
            df_diff_method.sort_values(by='y_std', ascending=True, inplace=True)
            
            if method == "wilcox":
                col_dot = [col_pal[3]]
            if method == "mast":
                col_dot = [col_pal[2]]
            if method == "trend": 
                col_dot = [col_pal[1]]
            if method == "epcy":
                col_dot = [col_pal[0]]

            df_diff_method['y_std_binned'] = pd.cut(df_diff_method['y_std'], bins=4).astype(str)
            intervals = df_diff_method['y_std_binned'].unique()
            size_mapping = {intervals[0]: 20, intervals[1]: 80, intervals[2]: 140, intervals[3]: 200}

            plt_fig = sns.relplot(
                x="mean_l2fc", y="mean_value", hue="method",
                row="cellType", col="cellNumber", size="y_std_binned", sizes=size_mapping,
                col_order=[3000, 5000, 8000, 10000],
                palette=col_dot,
                data=df_diff_method
            )

            fig_out = os.path.join(fig_dir, method + "_volcano_4.pdf")
            plt_fig.savefig(fig_out)
            plt.close()

    col_cutoff = sns.color_palette("flare")

    for cellNumber in search_params_dict['cellNumber']:
        fig_dir = os.path.join(args.OUTDIR, "eval_sc", str(cellNumber))
        if not os.path.exists(fig_dir):
            os.makedirs(fig_dir)

        for method in search_params_dict['methods']:
            for cellType in search_params_dict['cellTypes']:
                print('Volcano: cellNumber: {}, cellType: {}, method: {}'.format(str(cellNumber), cellType, method))
                df_diff_method = df_diff[(df_diff["cellNumber"] == cellNumber) & (df_diff["cellType"] == cellType) & (df_diff["method"] == method)].copy()
                df_diff_method.sort_values(by='value', ascending=True, inplace=True)

                if method == "wilcox":
                    col_dot = [col_pal[3]]
                if method == "mast":
                    col_dot = [col_pal[2]]
                if method == "trend":
                    col_dot = [col_pal[1]]
                if method == "epcy":
                    col_dot = [col_pal[0]]

                plt_fig = sns.relplot(
                    x="l2fc", y="value", hue="method",
                    col="rep", col_wrap=4, 
                    col_order=search_params_dict['reps'],
                    palette=col_dot,
                    data=df_diff_method
                )

                df_cutoff_tmp = df_cutoff[(df_cutoff["cellNumber"] == cellNumber) & (df_cutoff["method"] == method)].copy()
                
                for ax in plt_fig.axes:
                    for row in df_cutoff_tmp.iterrows():
                        cutoff = row[1]["value"]
                        quantile = row[1]["quantile"]
                        match quantile:
                            case 0.9999:
                                color = col_cutoff[5]
                            case 0.9995:
                                color = col_cutoff[3]
                            case 0.999:
                                color = col_cutoff[4]
                            case 0.995:
                                color = col_cutoff[2]
                            case 0.990:
                                color = col_cutoff[1]
                        ax.axhline(y=cutoff, color=color, linestyle='--')
                    rep = float(ax.get_title().replace("rep = ", ""))
                    df_tmp = df_diff_method[df_diff_method['rep'] == rep].copy()
                    df_tmp.sort_values(by=['value','l2fc'], ascending=True, inplace=True)
                    
                    i = 0
                    for idx,row in df_tmp.tail(5).iterrows():
                        x = row["l2fc"]
                        y = row["value"]
                        text = row["ID"]
                        text = df_id2name[df_id2name["gene_id"] == text]["gene_name"].iloc[0]
                        random.seed(text)
                        if 'epcy' in method:
                            offsets_x = np.arange(-0.4, 0.4, 0.01).tolist()
                            offsets_y = np.arange(-0.2, 0.2, 0.01).tolist()
                        else:
                            offsets_x = np.arange(-0.5, 0.5, 0.01).tolist()
                            if method == "wilcox":
                                offsets_x = np.arange(-2, 2, 0.05).tolist()
                            if method == "mast":
                                offsets_x = np.arange(-1, 1, 0.02).tolist()
                            offsets_y = np.arange(-40, -10, 1).tolist()

                        offset_x = offsets_x[int(len(offsets_x) / 10)*i]
                        offset_y = random.choice(offsets_y)

                        ax.plot([x, x + offset_x], [y, y+offset_y], color=col_dot[0], zorder=4)
                        ax.text(
                            x+offset_x ,y+offset_y,text, 
                            horizontalalignment='center',
                            verticalalignment='center',
                        )
                        i += 1
                    del df_tmp

                fig_out = os.path.join(fig_dir, method + "_" + cellType + "_volcano.pdf")
                plt_fig.savefig(fig_out)
                fig_out = os.path.join(fig_dir, method + "_" + cellType + "_volcano.png")
                plt_fig.savefig(fig_out)
                plt.close()

                del df_diff_method, df_cutoff_tmp

    def plot_volcano(ax, df_diff, df_cutoff, cellNumber, cellType, method, rep, num_gene_display=5):
        df_diff_method = df_diff[(df_diff["cellNumber"] == cellNumber) & (df_diff["cellType"] == cellType) & (df_diff["method"] == method) & (df_diff["rep"] == rep)].copy()
        df_diff_method.sort_values(by=['value','l2fc'], ascending=True, inplace=True)
        print('Volcano: cellNumber: {}, cellType: {}, method: {}, rep: {}'.format(str(cellNumber), cellType, method, rep))

        if method == "wilcox":
            col_dot = [col_pal[3]]
        if method == "mast":
            col_dot = [col_pal[2]]
        if method == "trend":
            col_dot = [col_pal[1]]
        if method == "epcy":
            col_dot = [col_pal[0]]

        sns.scatterplot(
            x="l2fc", y="value", hue="method",
            palette=col_dot,
            data=df_diff_method, ax=ax
        )

        if "epcy" in method:
            ax.set(ylim=(-0.1, 1))
        else:
            ax.set(ylim=(0, 315))

        i = 0
        for idx,row in df_diff_method.tail(num_gene_display).iterrows():
            x = row["l2fc"]
            y = row["value"]
            text = row["ID"]
            text = df_id2name[df_id2name["gene_id"] == text]["gene_name"].iloc[0]
            random.seed(text)
            if 'epcy' in method:
                offsets_x = np.arange(-0.4, 0.4, 0.01).tolist()
                offsets_y = np.arange(-0.2, 0.2, 0.01).tolist()
            else:
                offsets_x = np.arange(-0.5, 0.5, 0.01).tolist()
                if method == "wilcox":
                    offsets_x = np.arange(-2, 2, 0.05).tolist()
                if method == "mast":
                    offsets_x = np.arange(-1, 1, 0.02).tolist()
                offsets_y = np.arange(-40, -10, 1).tolist()

            offset_x = offsets_x[int(len(offsets_x) / 10)*i]
            offset_y = random.choice(offsets_y)

            ax.plot([x, x + offset_x], [y, y+offset_y], color=col_dot[0], zorder=4)
            ax.text(
                x+offset_x ,y+offset_y,text, 
                horizontalalignment='center',
                verticalalignment='center',
            )
            i += 1
    
        df_cutoff_tmp = df_cutoff[(df_cutoff["cellNumber"] == cellNumber) & (df_cutoff["method"] == method)].copy()
        
        for row in df_cutoff_tmp.iterrows():
            cutoff = row[1]["value"]
            quantile = row[1]["quantile"]
            match quantile:
                case 0.9999:
                    color = col_cutoff[5]
                case 0.9995:
                    color = col_cutoff[3]
                case 0.999:
                    color = col_cutoff[4]
                case 0.995:
                    color = col_cutoff[2]
                case 0.990:
                    color = col_cutoff[1]
            ax.axhline(y=cutoff, color=color, linestyle='--')
        del df_cutoff_tmp, df_diff_method

    def fig_volcano(cellType, df_diff, df_cutoff):
        print('*** Fig Volcano on cellType: {}'.format(cellType))
        fig, axs = plt.subplots(ncols=4, nrows=4, sharey='row', figsize=(30, 18), dpi=600 )

        cellNumber = 3000
        for i, method in enumerate(search_params_dict['methods']):
            for j, rep in enumerate([1,12]):
                ax = axs[i,j]
                plot_volcano(ax, df_diff, df_cutoff, cellNumber, cellType, method, rep)

        cellNumber = 5000
        rep = 1
        for i, method in enumerate(search_params_dict['methods']):
            ax=axs[i, 2]
            plot_volcano(ax, df_diff, df_cutoff, cellNumber, cellType, method, rep)
        
        cellNumber = 10000
        rep = 1
        for i, method in enumerate(search_params_dict['methods']):
            ax =axs[i, 3]
            plot_volcano(ax, df_diff, df_cutoff, cellNumber, cellType, method, rep)

        fig_dir = os.path.join(args.OUTDIR, "eval_sc")
        fig_out = os.path.join(fig_dir, cellType+ "_fig_volcano.pdf")
        print(fig_out)
        fig.savefig(fig_out)
        fig_out = os.path.join(fig_dir, cellType+ "_fig_volcano.png")
        fig.savefig(fig_out)
        print(fig_out)
        plt.close()

        cellNumber = 3000
        for i, method in enumerate(search_params_dict['methods']):
            fig, axs = plt.subplots(ncols=2, nrows=1, sharey='row', figsize=(8, 4), dpi=600)
            for j, rep in enumerate([1,12]):
                plot_volcano(axs[j], df_diff, df_cutoff, cellNumber, cellType, method, rep)
            fig_out = os.path.join(fig_dir, method + "_" + cellNumber + "_" + cellType + "1_12_fig_volcano.pdf")
            fig.savefig(fig_out)
            plt.close()

        rep = 1
        for i, method in enumerate(search_params_dict['methods']):
            fig, axs = plt.subplots(ncols=2, nrows=1, sharey='row', figsize=(8, 4), dpi=600 )
            for j, cellNumber in enumerate([5000,10000]):
                plot_volcano(axs[j], df_diff, df_cutoff, cellNumber, cellType, method, rep)
            fig_out = os.path.join(fig_dir, method + "_" + rep + "_" + cellType + "5000_10000_fig_volcano.pdf")
            fig.savefig(fig_out)
            plt.close()
        
    for cellType in search_params_dict['cellTypes']:
        fig_volcano(cellType, df_diff, df_cutoff)

    ###### Pairwise plot

    def hexbin(x, y, color, **kwargs):
        cmap = sns.light_palette(color, as_cmap=True)
        plt.hexbin(x, y, gridsize=15, cmap=cmap, norm=LogNorm(),**kwargs)
        plt.colorbar()

    all_df_long = []
    for cellType in search_params_dict['cellTypes']:
        list_df = []
        for cellNumber in search_params_dict['cellNumber']:
            fig_dir = os.path.join(args.OUTDIR, "eval_sc", str(cellNumber))
            if not os.path.exists(fig_dir):
                os.makedirs(fig_dir)
            for method in search_params_dict['methods']:
                print('Pairwise: cellNumber: {}, cellType: {}, method: {}'.format(str(cellNumber), cellType, method))
                df = df_diff[(df_diff["cellNumber"] == cellNumber) & (df_diff["cellType"] == cellType) & (df_diff["method"] == method)].copy()
                df = df.pivot(index='ID', columns='rep', values='value')
                df = df.reset_index()
                df.dropna(inplace=True)
                df_long = df.copy()

                df["method"] = [method] * df.shape[0]

                df_long['ID'] = range(1, df_long.shape[0]+1)
                df_long.columns = ['x' + str(col) for col in df_long.columns]
                df_long = df_long.rename(columns={df_long.columns[0]: 'ID'})
                df_long = df_long.rename(columns={df_long.columns[1]: 'y'})
                df_long = pd.wide_to_long(df_long, stubnames='x', i='ID', j='y')
                df_long['method'] = [method] * df_long.shape[0]
                df_long['cellNumber'] = [cellNumber] * df_long.shape[0]
                list_df.append(df_long)
                
                lim_value = (0, 310)
                label_value = "-log10(padj)"
                if method == "epcy":
                    lim_value = (-0.20, 1)
                    label_value = "MCC"

                if method == "wilcox":
                    col_dot = [col_pal[3]]
                if method == "mast":
                    col_dot = [col_pal[2]]
                if method == "trend":
                    col_dot = [col_pal[1]]
                if method == "epcy":
                    col_dot = [col_pal[0]]

                for pair in  [[1,12],[3,14],[11,18]]:
                    df.sort_values(by=pair[0], ascending=True, inplace=True)

                    # to build square scatterplot
                    plt.figure(figsize=(5, 5))

                    plt_fig = sns.scatterplot(
                        data=df, x=pair[1], y=pair[0], hue="method",
                        palette=col_dot
                    )

                    xordinal = 'th'
                    if pair[1] == 1:
                        xordinal = 'st'
                    if pair[1] == 2:
                        xordinal = 'nd'
                    if pair[1] == 3:
                        xordinal = 'rd'
                    yordinal = 'th'
                    if pair[0] == 1:
                        yordinal = 'st'
                    if pair[0] == 2:
                        yordinal = 'nd'
                    if pair[0] == 3:
                        yordinal = 'rd'
                    plt.xlabel(str(pair[1]) + xordinal + " replicate's " + label_value)
                    plt.ylabel(str(pair[0]) + yordinal + " replicate's " + label_value)

                    plt.gca().set_aspect('equal', adjustable='box')
                    #plt.tight_layout()

                    plt_fig.set(xlim=lim_value, ylim=lim_value)
                    lim_diag = [lim_value[1], lim_value[0]]
                    plt_fig.plot(lim_diag, lim_diag, '-r')
                    fig_out = os.path.join(fig_dir, method + "_" + cellType + "_pairwise_" + str(pair[0]) + "_" + str(pair[1]) + ".png")
                    plt_fig.get_figure().savefig(fig_out)
                    fig_out = os.path.join(fig_dir, method + "_" + cellType + "_pairwise_" + str(pair[0]) + "_" + str(pair[1]) + ".pdf")
                    plt_fig.get_figure().savefig(fig_out)
                    plt.close()

                del df, df_long
                gc.collect()
        df_long = pd.concat(list_df, axis=0, ignore_index=True)
        del list_df
        gc.collect()

        plt_fig = sns.FacetGrid(
            df_long, col="cellNumber", row="method", 
            hue="method", palette=col_pal_cmap,
            sharex=False, sharey=False,
            gridspec_kws={"hspace":0.4, "wspace":0.5},
            height=3, aspect=1.4
        )
        plt_fig.map(hexbin, "x", "y")
        plt_fig.set_titles("{row_name} using {col_name} cells")
        plt.subplots_adjust(top=0.92) # adjust the Figure to display title
        plt.suptitle("Hexbins to compare analyses reproducibility on " + cellType + " cells", fontsize=14)
        for ax in plt_fig.axes.flatten():
            if 'epcy' in ax.get_title():
                lim_value = (-0.20, 1)
                ax.set_ylabel("1st rep's MCC")
                ax.set_xlabel("Other rep's MCC")
            else:
                lim_value = (0, 310)
                ax.set_ylabel("1st rep's -log10(padj)")
                ax.set_xlabel("Other rep's -log10(padj)")
            ax.set(xlim=lim_value, ylim=lim_value)
            lim_diag = [lim_value[1], lim_value[0]]
            ax.plot(lim_diag, lim_diag, '-r')

        fig_dir = os.path.join(args.OUTDIR, "eval_sc")
        fig_out = os.path.join(fig_dir, cellType + "_pairwise.png")
        plt_fig.savefig(fig_out)
        fig_out = os.path.join(fig_dir, cellType + "_pairwise.pdf")
        plt_fig.savefig(fig_out, bbox_inches='tight')
        plt.close()
        all_df_long.append(df_long)
        del df_long
        gc.collect()
    
    df_long = pd.concat(all_df_long, axis=0, ignore_index=True)
    del all_df_long
    gc.collect()

    plt_fig = sns.FacetGrid(
        df_long, col="cellNumber", row="method", 
        hue="method", palette=col_pal_cmap,
        sharex=False, sharey=False,
        gridspec_kws={"hspace":0.4, "wspace":0.5},
        height=3, aspect=1.4
    )
    plt_fig.map(hexbin, "x", "y")
    plt_fig.set_titles("{row_name} using {col_name} cells")
    plt.subplots_adjust(top=0.92) # adjust the Figure to display title
    plt.suptitle("Hexbins to compare analyses reproducibility of all cells type", fontsize=14)
    for ax in plt_fig.axes.flatten():
        if 'epcy' in ax.get_title():
            lim_value = (-0.20, 1)
            ax.set_ylabel("1st rep's MCC")
            ax.set_xlabel("Other rep's MCC")
        else:
            lim_value = (0, 310)
            ax.set_ylabel("1st rep's -log10(padj)")
            ax.set_xlabel("Other rep's -log10(padj)")
        ax.set(xlim=lim_value, ylim=lim_value)
        lim_diag = [lim_value[1], lim_value[0]]
        ax.plot(lim_diag, lim_diag, '-r')

    fig_dir = os.path.join(args.OUTDIR, "eval_sc")
    fig_out = os.path.join(fig_dir, "pairwise.png")
    plt_fig.savefig(fig_out)
    fig_out = os.path.join(fig_dir, "pairwise.pdf")
    plt_fig.savefig(fig_out)
    plt.close()
    del df_long
    gc.collect()


