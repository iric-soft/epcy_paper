import datetime
# get time and date
start_time = datetime.datetime.now()
# make rest of imports
import os
from functools import reduce
from collections import defaultdict

import pdb

import pandas as pd
import numpy as np
from scipy import stats
import upsetplot as up

import matplotlib as mpl
# set no display
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from decimal import Decimal
import subprocess

import umap
from sklearn.decomposition import PCA

from .. utils import other as uo
#from utils import other as uo

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'DejaVu Sans'
mpl.rcParams['pdf.fonttype'] = 42

# for test
# from types import SimpleNamespace
# args = SimpleNamespace(**{ 'N_GENES': 60564, 'CPM': 1, 'MCC': -2, 'LOG_FC':0, 'TYPE_QUANT': 'readcounts', 'QUANT': 'STAR_RSEM', 'DESIGNS':["30_t15_17", "30_inv16"], 'QUANTILES': [0.9999, 0.9995, 0.999, 0.995, .99], 'METHODS': ['deseq2', 'edger', 'voom', 'epcy'], 'PATH': '/u/eaudemard/project/epcy_paper/data/design/leucegene3/', 'OUTDIR': '/u/eaudemard/project/epcy_paper/data/res', 'SUBGROUP': 'subgroup', 'BF': '/u/eaudemard/project/epcy_paper/data/other/GRCh38_84_genes_biotype.tsv', 'MATRIX': '/u/eaudemard/project/epcy_paper/data/leucegene3/STAR_RSEM/readcounts.xls' })
# args.BIOTYPE = None

def main_eval_bulk(args, argparser):

    file_dict = {
        'epcy' : "predictive_capability.tsv",
        'deseq2' : "deseq2_genes.xls",
        'edger' : "edger_genes.xls",
        'voom' : "limma_voom_genes.xls",
        'trend' : "limma_trend_genes.xls",
        'mast' : "mast_genes.xls",
    }

    # create dict with search params
    search_params_dict = {
        'methods' : args.METHODS,
        'designs' : args.DESIGNS
    }

    strl2fc_dict = {
        'epcy' : "L2FC",
        'deseq2' : "log2FoldChange",
        'edger' : "logFC",
        'voom' : "logFC",
        'trend' : "logFC",
    }

    strvalue_dict = {
        'epcy' : "KERNEL_MCC",
        'deseq2' : "pvalue",
        'edger' : "PValue",
        'voom' : "P.Value",
        'trend' : "P.Value",
        'mast' : "pval",
    }

    strvalueadj_dict = {
        'epcy' : "KERNEL_MCC",
        'deseq2' : "padj",
        'edger' : "FDR",
        'voom' : "adj.P.Val",
        'trend' : "adj.P.Val",
        'mast' : "pval",
    }

    cohort_order = [*search_params_dict['designs']]
    print(cohort_order)

    df_biotype = pd.read_csv(args.BF, sep="\t")
    df_biotype['ensembl_gene_id'] = df_biotype['ensembl_gene_id'].str.split('.').str[0]

    if args.BIOTYPE is not None:
        selected_biotype = args.BIOTYPE.split(",")
        df_biotype = df_biotype.loc[df_biotype["gene_biotype"].isin(selected_biotype)]

    df_exp = uo.get_exp(args, args.MATRIX)

    df_all_res = []
    for design in search_params_dict['designs']:
        dir_design = os.path.join(args.PATH, design)
        for method in search_params_dict['methods']:
            print('RANDOM design: {}, method: {}'.format(design, method))

            df_tmp = uo.read_diff_table(args, file_dict[method], method, dir_design, 1)

            if method != "epcy":
                df_tmp = df_tmp.reindex(df_tmp[strvalue_dict[method]].abs().sort_values(ascending=True).index)

            df_all_res.append(
                pd.DataFrame({
                    "ID": df_tmp["ID"],
                    "METHOD": method,
                    "DESIGN": design,
                    "TYPE": "MCC" if method == "epcy" else "PVALUE",
                    "VALUE": df_tmp[strvalue_dict[method]] if method == "epcy" else -np.log10(df_tmp[strvalue_dict[method]]),
                    "TOP": [x for x in range(1, df_tmp.shape[0]+1, 1)]
                })
            )

            if method != "epcy":
                df_tmp = df_tmp.reindex(df_tmp[strvalueadj_dict[method]].abs().sort_values(ascending=True).index)

                df_all_res.append(
                    pd.DataFrame({
                        "ID": df_tmp["ID"],
                        "METHOD": method,
                        "DESIGN": design,
                        "TYPE": "PADJ",
                        "VALUE": -np.log10(df_tmp[strvalueadj_dict[method]]),
                        "TOP": [x for x in range(1, df_tmp.shape[0]+1, 1)]
                    })
                )

    print('MERGE ALL DATA')
    df_all_res = pd.concat(df_all_res)
    df_all_res = df_all_res.replace([np.inf, -np.inf], 340)

    # create subfolder for log and res
    fig_dir = os.path.join(args.OUTDIR, "eval_bulk")

    if not os.path.exists(fig_dir):
        os.makedirs(fig_dir)

    print('SAVE TABLE')
    csv_out = os.path.join(fig_dir, "bulk_table.csv")
    df_all_res.to_csv(csv_out, index=False, sep="\t")
    # df_all_res = pd.read_csv(csv_out, sep="\t")
    df_all_res['ID'] = df_all_res['ID'].str.split('.').str[0]

    col_pal = sns.color_palette("colorblind")
    col_pal = [col_pal[i] for i in [4, 1, 2, 9]]

    print('PLOT FIG')
    res_cutoff = []
    for method in search_params_dict['methods']:
        df_tmp = df_all_res.loc[df_all_res["TYPE"] != "PVALUE"]
        df_tmp = df_tmp.loc[df_tmp["METHOD"] == method]
        for design in search_params_dict['designs']:
            df_tmp_design = df_tmp.loc[df_tmp["DESIGN"] == design]
            df_tmp_design = df_tmp_design.sort_values(by='VALUE', ascending=False)
            for quantile in args.QUANTILES:
                index_cutoff = int(args.N_GENES * (1-quantile))
                cutoff = df_tmp_design["VALUE"].iloc[index_cutoff]
                d = {'quantile': [quantile], 'design': [design], 'num_gene': [index_cutoff], 'METHOD': [method], 'cutoff': [cutoff]}
                res_cutoff.append(pd.DataFrame(data=d))
    res_cutoff = pd.concat(res_cutoff)


    csv_out = os.path.join(fig_dir, "tab_cutoff.csv")
    res_cutoff = pd.pivot_table(res_cutoff, values='cutoff', index=['quantile', 'num_gene', 'design'], columns=['METHOD'])
    res_cutoff.reset_index(drop=False, inplace=True)
    res_cutoff.to_csv(csv_out, index=False, sep="\t")
    # res_cutoff = pd.read_csv(csv_out, sep="\t")
    
    genes_highlight = [
        "PDPN", "ARHGAP4", "RASL12", 
        'HRASLS5', #PLAAT5
        "UNCX", "GABRE", "CDK18", "NCR2", "HOXA9",
        "AC025811.3"
        
    ]
    markers = {
        "default": "o", 
        "PDPN": "v", "ARHGAP4": "<", "RASL12": "D",
        'HRASLS5': "d", 
        "UNCX": "^",
        "GABRE": "s", "CDK18": "X",
        "NCR2": ">", "HOXA9": "P", 
        "AC025811.3": "*"
    }

    linestyles=[ # (0, (1, 1)), #dotted
        (0, (3, 5, 1, 5)), (0, (3, 5, 1, 5, 1, 5)), # 'dot dashed
        (0, (5, 10)), (0, (3, 1, 1, 1, 1, 1)), # 'densely dot dashed'
        (0, (5, 1)) # 'densely dashed'
    ]
    col_pal_condition = [
        mpl.colors.hex2color('#D21417'),
        mpl.colors.hex2color('#1483D2')
    ]

    for design in search_params_dict['designs']:
        df_all_tmp = df_all_res.loc[df_all_res["DESIGN"] == design]
        df_epcy = df_all_tmp.loc[df_all_tmp["METHOD"] == "epcy"]
        df_epcy = df_epcy.rename(columns={"VALUE": "epcy"})
        df_deg = df_all_tmp.loc[df_all_tmp["METHOD"] != "epcy"]
        deg_methods = pd.unique(df_deg.METHOD)

        for deg_m in deg_methods:
            df_deg_m = df_all_tmp.loc[df_all_tmp["METHOD"] == deg_m]
            df_deg_m = df_deg_m.loc[df_deg_m["TYPE"] == "PADJ"]
            df_deg_m = df_deg_m.rename(columns={"VALUE": deg_m})

            df_merge = pd.merge(df_epcy, df_deg_m, on='ID', how='outer')
            df_merge = pd.merge(df_merge, df_biotype, left_on='ID', right_on='ensembl_gene_id', how='left')
            df_merge["shape"] = "default"
            if design == "30_t15_17":
                for gene in genes_highlight:
                    df_merge.loc[df_merge["external_gene_name"] == gene, "shape"] = gene
            sns_plot = sns.scatterplot(
                x="epcy", y=deg_m, style="shape", linewidth=0,
                data=df_merge[df_merge["shape"] == "default"], 
                color="black", markers=markers, style_order=["default"],   
            )
            if df_merge[df_merge["shape"] != "default"].shape[0] > 0:
                sns_plot = sns.scatterplot(
                    x="epcy", y=deg_m, style="shape", linewidth=1, edgecolor="red",
                    data=df_merge[df_merge["shape"] != "default"], 
                    color="black", markers=markers, style_order=genes_highlight,   
                )
            sns.move_legend(sns_plot, "upper left", bbox_to_anchor=(1, 1))

            i_quant = 0
            for quantile in args.QUANTILES:
                epcy_cutoff = res_cutoff.loc[(res_cutoff["design"] == design) & (res_cutoff["quantile"] == quantile)]["epcy"].values[0]
                deg_cutoff = res_cutoff.loc[(res_cutoff["design"] == design) & (res_cutoff["quantile"] == quantile)][deg_m].values[0]  
                sns_plot.axes.axvline(epcy_cutoff, ls=linestyles[i_quant], color='black')
                sns_plot.axes.axhline(deg_cutoff, ls=linestyles[i_quant], color='black')
                i_quant += 1

            fig_out = os.path.join(fig_dir, design + "_epcy_vs_" + deg_m + ".pdf")
            sns_plot.figure.savefig(fig_out, bbox_inches='tight')
            plt.close()
            print("DONE: " + design + "_epcy_vs_" + deg_m)

    for quantile in args.QUANTILES:
        for design in search_params_dict['designs'][::-1]:
            dir_design = os.path.join(args.PATH, design)
            df_design = uo.get_design(args, "Query", dir_design)
            df_upset_raw = []
            for method in search_params_dict['methods']:
                print('RANDOM design: {}, method: {}'.format(design, method))
                df_tmp = uo.read_diff_table(args, file_dict[method], method, dir_design, 1)
                if method != "epcy":
                    df_tmp[strvalueadj_dict[method]] = -np.log10(df_tmp[strvalueadj_dict[method]])
                num_selected = res_cutoff.loc[(res_cutoff["design"] == design) & (res_cutoff["quantile"] == quantile)]["num_gene"].values[0]
                df_tmp.sort_values(by=strvalueadj_dict[method], ascending=False, inplace=True)
                df_tmp_selected = df_tmp.iloc[0:num_selected]
                df_upset_raw.append(
                    pd.DataFrame({
                        method: [True] * num_selected,
                        'ID': df_tmp['ID'][0:num_selected].values,
                    })
                )
            df_upset_raw = reduce(lambda left,right: pd.merge(left,right,on='ID', how='outer'), df_upset_raw)
            df_upset_raw = df_upset_raw.fillna(False).infer_objects(copy=False)
            df_upset = df_upset_raw.groupby(search_params_dict['methods']).size()
            up.plot(df_upset)
            fig_out = os.path.join(fig_dir, design + "_upset_" + str(quantile) + ".pdf")
            plt.savefig(fig_out)
            plt.close()

            def display_profile(df, method, design, fig_dir, quantile):
                for id in df["ID"]:
                    command = 'epcy profile_rna -d {} --condition subgroup -m {} --log --cpm -o {} --no_density --ids {}'.format(
                        os.path.join(args.PATH, design, "design.tsv"),
                        args.MATRIX,
                        os.path.join(fig_dir, "profile", method, str(quantile)),
                        id
                    )
                    print(command)
                    process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
                    output, error = process.communicate()

            def compute_result (df, df_design, design, df_exp):
                df_method_exp = df_exp.loc[df["ID"]]
                df_mehod_exp_query = df_method_exp[df_design[df_design["subgroup"] == 1]["sample"]]
                query_variance = df_mehod_exp_query.var(axis=1).to_frame(name='variance')
                query_variance["subgroup"] = design
                query_mean = df_mehod_exp_query.mean(axis=1).to_frame(name='mean')
                query_mean["subgroup"] = design
                df_method_exp_ref = df_method_exp[df_design[df_design["subgroup"] == 0]["sample"]]
                ref_variance = df_method_exp_ref.var(axis=1).to_frame(name='variance')
                ref_variance["subgroup"] = "Other AML"
                ref_mean = df_method_exp_ref.mean(axis=1).to_frame(name='mean')
                ref_mean["subgroup"] = "Other AML"
                df_fc = (query_mean['mean'] - ref_mean['mean']).to_frame(name='log2FC')

                df_var = pd.concat([query_variance, ref_variance])
                df_var.reset_index(inplace=True)
                df_mean = pd.concat([query_mean, ref_mean]) 
                df_mean.reset_index(inplace=True)
                df_result = df_var.merge(df_mean, on=['ID', 'subgroup'])
                return(df_result, df_fc)

            df_epcy = df_upset_raw[
                (df_upset_raw["epcy"] == True) &
                (df_upset_raw["deseq2"] == False) &
                (df_upset_raw["edger"] == False) &
                (df_upset_raw["voom"] == False)
            ]
            df_epcy_res, df_epcy_fc = compute_result(df_epcy, df_design, design, df_exp)
            df_epcy_res['method'] = 'epcy'
            df_epcy_fc['method'] = 'epcy'
            df_exp_uniq = df_exp.loc[df_epcy_fc.index]
            df_exp_uniq["ID"] = df_exp_uniq.index
            df_exp_uniq = pd.melt(df_exp_uniq, id_vars=["ID"], value_vars=df_exp_uniq.columns[0:-1], value_name="log2(CPM(x)+1)")
            df_exp_uniq["condition"] = df_exp_uniq["variable"].apply(lambda x: design if df_design[df_design["sample"] == x]["subgroup"].iloc[0] == 1 else "Other AML")

            if quantile >= 0.999:
                sns_plot = sns.catplot(
                    data=df_exp_uniq, kind="swarm",
                    x="log2(CPM(x)+1)", y="condition", hue="condition", col="ID", col_wrap=4,
                    palette=sns.color_palette([col_pal_condition[1], col_pal_condition[0]])
                )
                fig_file = os.path.join(fig_dir, design + "_epcy_unique_" + str(quantile) + ".pdf")
                sns_plot.savefig(fig_file)
                plt.close('all')


            df_deseq2 = df_upset_raw[
                (df_upset_raw["epcy"] == False) &
                (df_upset_raw["deseq2"] == True) &
                (df_upset_raw["edger"] == False) &
                (df_upset_raw["voom"] == False)
            ]
            df_deseq2_res, df_deseq2_fc  = compute_result(df_deseq2, df_design, design, df_exp)
            df_deseq2_res['method'] = 'deseq2'
            df_deseq2_fc['method'] = 'deseq2'
            df_exp_uniq = df_exp.loc[df_deseq2_fc.index]
            df_exp_uniq["ID"] = df_exp_uniq.index
            df_exp_uniq = pd.melt(df_exp_uniq, id_vars=["ID"], value_vars=df_exp_uniq.columns[0:-1], value_name="log2(CPM(x)+1)")
            df_exp_uniq["condition"] = df_exp_uniq["variable"].apply(lambda x: design if df_design[df_design["sample"] == x]["subgroup"].iloc[0] == 1 else "Other AML")

            if quantile >= 0.999:
                sns_plot = sns.catplot(
                    data=df_exp_uniq, kind="swarm",
                    x="log2(CPM(x)+1)", y="condition", hue="condition", col="ID", col_wrap=4,
                    palette=sns.color_palette([col_pal_condition[1], col_pal_condition[0]])
                )
                fig_file = os.path.join(fig_dir, design + "_deseq2_unique_" + str(quantile) + ".pdf")
                sns_plot.savefig(fig_file)
                plt.close('all')



            df_edger = df_upset_raw[
                (df_upset_raw["epcy"] == False) &
                (df_upset_raw["deseq2"] == False) &
                (df_upset_raw["edger"] == True) &
                (df_upset_raw["voom"] == False)
            ]
            df_edger_res, df_edger_fc  = compute_result(df_edger, df_design, design, df_exp)
            df_edger_res['method'] = 'edger'
            df_edger_fc['method'] = 'edger'
            df_exp_uniq = df_exp.loc[df_edger_fc.index]
            df_exp_uniq["ID"] = df_exp_uniq.index
            df_exp_uniq = pd.melt(df_exp_uniq, id_vars=["ID"], value_vars=df_exp_uniq.columns[0:-1], value_name="log2(CPM(x)+1)")
            df_exp_uniq["condition"] = df_exp_uniq["variable"].apply(lambda x: design if df_design[df_design["sample"] == x]["subgroup"].iloc[0] == 1 else "Other AML")

            if quantile >= 0.999:
                sns_plot = sns.catplot(
                    data=df_exp_uniq, kind="swarm",
                    x="log2(CPM(x)+1)", y="condition", hue="condition", col="ID", col_wrap=4,
                    palette=sns.color_palette([col_pal_condition[1], col_pal_condition[0]])
                )
                fig_file = os.path.join(fig_dir, design + "_edger_unique_" + str(quantile) + ".pdf")
                sns_plot.savefig(fig_file)
                plt.close('all')



            df_voom = df_upset_raw[
                (df_upset_raw["epcy"] == False) &
                (df_upset_raw["deseq2"] == False) &
                (df_upset_raw["edger"] == False) &
                (df_upset_raw["voom"] == True)
            ]
            df_voom_res, df_voom_fc  = compute_result(df_voom, df_design, design, df_exp)
            df_voom_res['method'] = 'voom'
            df_voom_fc['method'] = 'voom'
            df_exp_uniq = df_exp.loc[df_voom_fc.index]
            df_exp_uniq["ID"] = df_exp_uniq.index
            df_exp_uniq = pd.melt(df_exp_uniq, id_vars=["ID"], value_vars=df_exp_uniq.columns[0:-1], value_name="log2(CPM(x)+1)")
            df_exp_uniq["condition"] = df_exp_uniq["variable"].apply(lambda x: design if df_design[df_design["sample"] == x]["subgroup"].iloc[0] == 1 else "Other AML")

            if quantile >= 0.999:
                sns_plot = sns.catplot(
                    data=df_exp_uniq, kind="swarm",
                    x="log2(CPM(x)+1)", y="condition", hue="condition", col="ID", col_wrap=4,
                    palette=sns.color_palette([col_pal_condition[1], col_pal_condition[0]])
                )
                fig_file = os.path.join(fig_dir, design + "_voom_unique_" + str(quantile) + ".pdf")
                sns_plot.savefig(fig_file)
                plt.close('all')



            df_res = pd.concat([df_deseq2_res,  df_edger_res, df_voom_res, df_epcy_res])
            df_fc = pd.concat([df_deseq2_fc, df_edger_fc, df_voom_fc, df_epcy_fc])
            
            sns_plot = sns.boxplot(
                data=df_res, x="method", y="variance", hue="subgroup",
                linewidth=1, #palette=col_pal
            )
            fig_out = os.path.join(fig_dir, design + "_variance_unique_" + str(quantile) + ".pdf")
            sns_plot.figure.savefig(fig_out)
            plt.close('all')

            sns_plot = sns.boxplot(
                data=df_res, x="method", y="mean", hue="subgroup",
                linewidth=1, #palette=col_pal
            )
            fig_out = os.path.join(fig_dir, design + "_mean_unique_" + str(quantile) + ".pdf")
            sns_plot.figure.savefig(fig_out)
            plt.close('all')

            sns_plot = sns.boxplot(
                data=df_fc, x="method", y="log2FC", hue="method",
                linewidth=1, palette=col_pal
            )
            fig_out = os.path.join(fig_dir, design + "_fc_unique_" + str(quantile) + ".pdf")
            sns_plot.figure.savefig(fig_out)
            plt.close('all')

            sns_plot = sns.violinplot(
                data=df_res, x="method", y="variance", hue='subgroup', density_norm='width',
                split=True, cut=0, gap=.1, inner="quart" #palette=col_pal
            )
            fig_out = os.path.join(fig_dir, design + "_variance_unique_violin_" + str(quantile) + ".pdf")
            sns_plot.figure.savefig(fig_out)
            plt.close('all')

            sns_plot = sns.violinplot(
                data=df_res, x="method", y="mean", hue='subgroup', density_norm='width',
                split=True, cut=0, gap=.1, inner="quart" #palette=col_pal
            )
            fig_out = os.path.join(fig_dir, design + "_mean_unique_violin_" + str(quantile) + ".pdf")
            sns_plot.figure.savefig(fig_out)
            plt.close('all')
            

            df_epcy = df_upset_raw[
                (df_upset_raw["epcy"] == True)
            ]
            df_epcy_res, df_epcy_fc = compute_result(df_epcy, df_design, design, df_exp)
            df_epcy_res['method'] = 'epcy'
            df_epcy_fc['method'] = 'epcy'

            df_deseq2 = df_upset_raw[
                (df_upset_raw["deseq2"] == True)
            ]
            df_deseq2_res, df_deseq2_fc  = compute_result(df_deseq2, df_design, design, df_exp)
            df_deseq2_res['method'] = 'deseq2'
            df_deseq2_fc['method'] = 'deseq2'

            df_edger = df_upset_raw[
                (df_upset_raw["edger"] == True)
            ]
            df_edger_res, df_edger_fc  = compute_result(df_edger, df_design, design, df_exp)
            df_edger_res['method'] = 'edger'
            df_edger_fc['method'] = 'edger'

            df_voom = df_upset_raw[
                (df_upset_raw["voom"] == True)
            ]
            df_voom_res, df_voom_fc  = compute_result(df_voom, df_design, design, df_exp)
            df_voom_res['method'] = 'voom'
            df_voom_fc['method'] = 'voom'

            df_res = pd.concat([df_deseq2_res,  df_edger_res, df_voom_res, df_epcy_res])
            df_fc = pd.concat([df_deseq2_fc, df_edger_fc, df_voom_fc, df_epcy_fc])
            
            sns_plot = sns.boxplot(
                data=df_res, x="method", y="variance", hue="subgroup",
                linewidth=1, #palette=col_pal
            )
            fig_out = os.path.join(fig_dir, design + "_variance_" + str(quantile) + ".pdf")
            sns_plot.figure.savefig(fig_out)
            plt.close('all')

            sns_plot = sns.boxplot(
                data=df_res, x="method", y="mean", hue="subgroup",
                linewidth=1, #palette=col_pal
            )
            fig_out = os.path.join(fig_dir, design + "_mean_" + str(quantile) + ".pdf")
            sns_plot.figure.savefig(fig_out)
            plt.close('all')

            sns_plot = sns.boxplot(
                data=df_fc, x="method", y="log2FC", hue="method",
                linewidth=1, palette=col_pal
            )
            fig_out = os.path.join(fig_dir, design + "_fc_" + str(quantile) + ".pdf")
            sns_plot.figure.savefig(fig_out)
            plt.close('all')

            sns_plot = sns.violinplot(
                data=df_res, x="method", y="variance", hue='subgroup', density_norm='width',
                split=True, cut=0, gap=.1, inner="quart" #palette=col_pal
            )
            fig_out = os.path.join(fig_dir, design + "_variance_violin_" + str(quantile) + ".pdf")
            sns_plot.figure.savefig(fig_out)
            plt.close('all')

            sns_plot = sns.violinplot(
                data=df_res, x="method", y="mean", hue='subgroup', density_norm='width',
                split=True, cut=0, gap=.1, inner="quart" #palette=col_pal
            )
            fig_out = os.path.join(fig_dir, design + "_mean_violin_" + str(quantile) + ".pdf")
            sns_plot.figure.savefig(fig_out)
            plt.close('all')

    for design in search_params_dict['designs'][::-1]:
        df_all_l2fc = []
        dir_design = os.path.join(args.PATH, design)
        df_epcy = uo.read_diff_table(args, file_dict["epcy"], "epcy", dir_design, 1)
        df_epcy.set_index('ID', inplace=True)
        df_epcy['ID'] = df_epcy.index
        for method in search_params_dict['methods']:
            print('RANDOM design: {}, method: {}'.format(design, method))
            df_tmp = uo.read_diff_table(args, file_dict[method], method, dir_design, 1)
            df_tmp.set_index('ID', inplace=True)
            df_tmp['ID'] = df_tmp.index
            if method != "epcy":
                df_tmp[strvalueadj_dict[method]] = -np.log10(df_tmp[strvalueadj_dict[method]])

            if method in ['voom', 'edger']:
                df_tmp[strl2fc_dict[method]] = -df_tmp[strl2fc_dict[method]]

            for quantile in args.QUANTILES:
                num_selected = res_cutoff.loc[(res_cutoff["design"] == design) & (res_cutoff["quantile"] == quantile)]["num_gene"].values[0]
                    
                df_tmp.sort_values(by=strvalueadj_dict[method], ascending=False, inplace=True)
                df_tmp_selected = df_tmp.iloc[0:num_selected]

                df_tmp_epcy = df_epcy[df_epcy["ID"].isin(df_tmp_selected["ID"])]
                df_tmp_epcy = df_tmp_epcy.loc[df_tmp_selected.index]
                df_all_l2fc.append(
                    pd.DataFrame({
                        'method': [method] * num_selected,
                        'ID': df_tmp['ID'][0:num_selected].values,
                        'log2FC': df_tmp[strl2fc_dict[method]][0:num_selected].values,
                        'log2FC_epcy': df_tmp_epcy[strl2fc_dict['epcy']][0:num_selected].values,
                        'quantile': [quantile] * num_selected,
                    })
                )  
        df_all_l2fc = pd.concat(df_all_l2fc) 

        for quantile in args.QUANTILES:

            for fc_str in ['log2FC_epcy']: #, 'log2FC']:
                epcy_deseq = stats.mannwhitneyu(
                    df_all_l2fc.loc[(df_all_l2fc['method'] == "epcy") & (df_all_l2fc['quantile'] == quantile), fc_str],
                    df_all_l2fc.loc[(df_all_l2fc['method'] == "deseq2") & (df_all_l2fc['quantile'] == quantile), fc_str]
                )[1]
                epcy_voom = stats.mannwhitneyu(
                    df_all_l2fc.loc[(df_all_l2fc['method'] == "epcy") & (df_all_l2fc['quantile'] == quantile), fc_str],
                    df_all_l2fc.loc[(df_all_l2fc['method'] == "voom") & (df_all_l2fc['quantile'] == quantile), fc_str]
                )[1]
                epcy_edger = stats.mannwhitneyu(
                    df_all_l2fc.loc[(df_all_l2fc['method'] == "epcy") & (df_all_l2fc['quantile'] == quantile), fc_str],
                    df_all_l2fc.loc[(df_all_l2fc['method'] == "edger") & (df_all_l2fc['quantile'] == quantile), fc_str]
                )[1]
                deseq_voom = stats.mannwhitneyu(
                    df_all_l2fc.loc[(df_all_l2fc['method'] == "deseq2") & (df_all_l2fc['quantile'] == quantile), fc_str],
                    df_all_l2fc.loc[(df_all_l2fc['method'] == "voom") & (df_all_l2fc['quantile'] == quantile), fc_str]
                )[1]
                deseq_edger = stats.mannwhitneyu(
                    df_all_l2fc.loc[(df_all_l2fc['method'] == "deseq2") & (df_all_l2fc['quantile'] == quantile), fc_str],
                    df_all_l2fc.loc[(df_all_l2fc['method'] == "edger") & (df_all_l2fc['quantile'] == quantile), fc_str]
                )[1]
                voom_edger = stats.mannwhitneyu(
                    df_all_l2fc.loc[(df_all_l2fc['method'] == "voom") & (df_all_l2fc['quantile'] == quantile), fc_str],
                    df_all_l2fc.loc[(df_all_l2fc['method'] == "edger") & (df_all_l2fc['quantile'] == quantile), fc_str]
                )[1]

                sns_plot = sns.boxplot(
                    x="method", y=fc_str, data=df_all_l2fc[df_all_l2fc['quantile']==quantile],
                    hue="method", linewidth=1, palette=col_pal
                )
                sns_plot.set_title(
                    "log2FC distribution of quantile " + str(quantile) + "\n" +
                    "epcy_deseq=" + "{:.2E}".format(Decimal(epcy_deseq)) + " " + 
                    "epcy_voom=" + "{:.2E}".format(Decimal(epcy_voom)) + " " + 
                    "epcy_edger=" + "{:.2E}".format(Decimal(epcy_edger)) + "\n" + 
                    "deseq_voom=" + "{:.2E}".format(Decimal(deseq_voom)) + " " +
                    "deseq_edger=" + "{:.2E}".format(Decimal(deseq_edger)) + " " +
                    "voom_edger=" + "{:.2E}".format(Decimal(voom_edger)) 
                )
                fig_out = os.path.join(fig_dir, "all_" + fc_str + "_" + str(quantile) + "_" + design + ".pdf")
                sns_plot.figure.savefig(fig_out)
                plt.close('all')

