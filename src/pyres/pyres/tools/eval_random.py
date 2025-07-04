import datetime
# get time and date
start_time = datetime.datetime.now()
# make rest of imports
import os

from collections import defaultdict

import pdb

import pandas as pd
import numpy as np
import upsetplot as up

import matplotlib as mpl
# set no display
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

from .. utils import other as uo

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'DejaVu Sans'
mpl.rcParams['pdf.fonttype'] = 42

def main_eval_random(args, argparser):

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
        'designs_random' : args.DESIGN_RANDOM,
        'designs' : args.DESIGN
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
        'epcy' : "AUC",
        'deseq2' : "padj",
        'edger' : "FDR",
        'voom' : "adj.P.Val",
        'trend' : "adj.P.Val",
        'mast' : "pval",
    }

    cohort_order = [*search_params_dict['designs_random'], *search_params_dict['designs']]
    print(cohort_order)
    top = 50

    df_biotype = None
    if args.BIOTYPE is not None:
        df_biotype = pd.read_csv(args.BF, sep="\t")
        selected_biotype = args.BIOTYPE.split(",")
        df_biotype = df_biotype.loc[df_biotype["gene_biotype"].isin(selected_biotype)]

    df_all_res = []
    for design in search_params_dict['designs_random']:
        dir_design = os.path.join(args.RANDOM_PATH, design)
        df_design = uo.get_design(args, "Query", dir_design)

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

    for design in search_params_dict['designs']:
        dir_design = os.path.join(args.PATH, design)
        df_design = uo.get_design(args, "Query", dir_design)

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
    fig_dir = os.path.join(args.OUTDIR, "eval_random")

    if not os.path.exists(fig_dir):
        os.makedirs(fig_dir)

    print('SAVE TABLE')
    csv_out = os.path.join(fig_dir, "random_table.csv")
    df_all_res.to_csv(csv_out, index=False, sep="\t")

    col_pal = sns.color_palette("colorblind")
    col_pal = [col_pal[i] for i in [4, 1, 2, 9]]

    print('PLOT FIG')
    res_fdr = []
    res_ids_fdr = []
    res_cutoff = []
    res_top = []
    for method in search_params_dict['methods']:
        df_tmp = df_all_res.loc[df_all_res["METHOD"] == method]
        df_tmp = df_tmp.loc[df_tmp["TOP"] <= 150]

        if method == "deseq2":
            col_dot = [col_pal[0]]
        if method == "mast":
            col_dot = [col_pal[0]]
        if method == "edger":
            col_dot = [col_pal[1]]
        if method == "voom":
            col_dot = [col_pal[2]]
        if method == "trend":
            col_dot = [col_pal[2]]
        if method == "epcy":
            col_dot = [col_pal[3]]

        if method != "epcy":
            # plot results pvalue
            df_tmp2 = df_tmp.loc[df_tmp["TYPE"] == "PVALUE"]
            df_tmp2 = df_tmp2.replace(np.inf, 350)
            plt_fig = sns.swarmplot(
                x="DESIGN", y="VALUE",
                data=df_tmp2,
                order=cohort_order,
                palette=col_dot
            )
            plt_fig.set(ylabel="-log10(pvalue)")
            #plt_fig.axhline(y=-np.log10(0.05), color='r', linestyle='--')
            plt_fig.set_xticklabels(plt_fig.get_xticklabels(), rotation=90)
            fig_out = os.path.join(fig_dir, method + "_pvalue.pdf")
            plt_fig.figure.savefig(fig_out)
            plt.close()

            # plot results pvalue adj
            df_tmp2 = df_tmp.loc[df_tmp["TYPE"] == "PADJ"]
            df_tmp2 = df_tmp2.replace(np.inf, 350)
            plt_fig = sns.swarmplot(
                x="DESIGN", y="VALUE",
                data=df_tmp2,
                order=cohort_order,
                palette=col_dot
            )
            plt_fig.set(ylabel="-log10(p_adjusted)")
            #plt_fig.axhline(y=-np.log10(0.05), color='r', linestyle='--')
            plt_fig.set_xticklabels(plt_fig.get_xticklabels(), rotation=90)
            fig_out = os.path.join(fig_dir, method + "_padj.pdf")
            plt_fig.figure.savefig(fig_out)
            plt.close()
        else:
            # plot results pvalue
            df_tmp2 = df_tmp.loc[df_tmp["TYPE"] == "MCC"]
            plt_fig = sns.swarmplot(
                x="DESIGN", y="VALUE",
                data=df_tmp2,
                order=cohort_order,
                palette=col_dot
            )
            plt_fig.set(ylabel="MCC")
            #plt_fig.axhline(y=0.2, color='r', linestyle='--')
            plt_fig.set_xticklabels(plt_fig.get_xticklabels(), rotation=90)
            fig_out = os.path.join(fig_dir, method + "_mcc.pdf")
            plt_fig.figure.savefig(fig_out)
            plt.close()

        df_all_res = df_all_res.loc[df_all_res["TYPE"] != "PVALUE"]
        df_tmp_random = df_all_res.loc[df_all_res["METHOD"] == method]
        df_tmp_random = df_tmp_random.loc[df_tmp_random["DESIGN"].isin(search_params_dict['designs_random'])]

        df_tmp_random = df_tmp_random.sort_values(by='VALUE', ascending=False)
        num_random = len(search_params_dict['designs_random'])
        for e_fpr in args.E_FPR:
            num_gene = int(args.N_GENES * num_random * e_fpr)
            cutoff = df_tmp_random["VALUE"].iloc[num_gene]
            df_tmp_candidate = df_all_res.loc[df_all_res["METHOD"] == method]
            df_tmp_candidate = df_tmp_candidate.loc[df_tmp_candidate["DESIGN"].isin(search_params_dict['designs'])]
            df_tmp_candidate = df_tmp_candidate.loc[df_tmp_candidate["VALUE"] >= cutoff]

            res = df_tmp_candidate["DESIGN"].value_counts().rename_axis('DESIGN').reset_index(name='counts')
            res["eFPR"] = e_fpr
            res["METHOD"] = method
            res_fdr.append(res)
            df_tmp_candidate["eFPR"] = e_fpr
            res_ids_fdr.append(df_tmp_candidate)

            d = {'eFPR': [e_fpr], 'num_gene': [num_gene], 'METHOD': [method], 'cutoff': [cutoff]}
            res_cutoff.append(pd.DataFrame(data=d))

        for design in search_params_dict['designs']:
            df_tmp_candidate = df_all_res.loc[df_all_res["METHOD"] == method]
            df_tmp_candidate = df_tmp_candidate.loc[df_tmp_candidate["DESIGN"] == design]
            df_tmp_candidate = df_tmp_candidate.sort_values(["VALUE"], ascending=[False])
            top3 = df_tmp_candidate["VALUE"].iloc[3]
            d = {'DESIGN': [design], 'num_gene': [3-1], 'METHOD': [method], 'cutoff': [top3]}
            res_top.append(pd.DataFrame(data=d))
            top_x = df_tmp_candidate["VALUE"].iloc[args.TOP-1]
            d = {'DESIGN': [design], 'num_gene': [args.TOP], 'METHOD': [method], 'cutoff': [top_x]}
            res_top.append(pd.DataFrame(data=d))

    res_fdr = pd.concat(res_fdr)
    res_cutoff = pd.concat(res_cutoff)
    res_top = pd.concat(res_top)
    res_ids_fdr = pd.concat(res_ids_fdr)


    csv_out = os.path.join(fig_dir, "tab_cutoff.csv")
    res_cutoff = pd.pivot_table(res_cutoff, values='cutoff', index=['eFPR', 'num_gene'], columns=['METHOD'])
    res_cutoff.reset_index(drop=False, inplace=True)
    res_cutoff.to_csv(csv_out, index=False, sep="\t")

    #res_cutoff.reindex_axis(['eFPR', 'num_gene', 'deseq2', 'edger', 'voom', "epcy"], axis=1)
    csv_out = os.path.join(fig_dir, "tab_top.csv")
    res_top = pd.pivot_table(res_top, values='cutoff', index=['num_gene', 'DESIGN'], columns=['METHOD'])
    res_top.reset_index(drop=False, inplace=True)
    res_top.to_csv(csv_out, index=False, sep="\t")

    markers=['o', 'v', 's', 'x']
    linestyles=['-', '--', '-.', ':']
    if len(args.METHODS) == 3:
        markers=['o', 's', 'x']
        linestyles=['-', '-.', ':']
        col_pal=[col_pal[0], col_pal[2], col_pal[3]]

    for design in args.DESIGN:
        res_fdr_tmp = res_fdr[res_fdr["DESIGN"] == design]
        sns_plot = sns.pointplot(
            x="eFPR", y="counts", hue="METHOD",
            hue_order=args.METHODS, order=args.E_FPR,
            data=res_fdr_tmp,
            linestyles=linestyles,
            markers=markers,
            palette=col_pal
        )
        sns_plot.set_title("Number of gene candidates\n function to % of FDR")
        fig_out = os.path.join(fig_dir, design + "_num_gene_candidates_fdr.pdf")
        sns_plot.figure.savefig(fig_out)
        plt.close()

    sns_plot = sns.catplot(
        kind="point", col="DESIGN",
        x="eFPR", y="counts", hue="METHOD",
        hue_order=args.METHODS, order=args.E_FPR,
        data=res_fdr,
        linestyles=linestyles,
        markers=markers,
        palette=col_pal
    )
    fig_out = os.path.join(fig_dir, "num_gene_candidates_fdr.pdf")
    sns_plot.savefig(fig_out)
    plt.close()

    if (len(args.DESIGN) >= 2):
        res_common = []
        for e_fpr in args.E_FPR:
            df_tmp = res_ids_fdr.loc[res_ids_fdr["eFPR"] == e_fpr]
            for method in args.METHODS:
                df_tmp_method = df_tmp.loc[df_tmp["METHOD"] == method]
                df_tmp1 = df_tmp_method.loc[df_tmp_method["DESIGN"] == args.DESIGN[0]]
                df_tmp2 = df_tmp_method.loc[df_tmp_method["DESIGN"] == args.DESIGN[1]]
                df_common = df_tmp1.loc[df_tmp1["ID"].isin(df_tmp2['ID'])]
                d = {'eFPR': [e_fpr], 'METHOD': [method], 'shared_genes': [df_common.shape[0]]}
                res_common.append(pd.DataFrame(data=d))

        res_common = pd.concat(res_common)

        sns_plot = sns.pointplot(
            kind="point", col="DESIGN",
            x="eFPR", y="shared_genes", hue="METHOD",
            hue_order=args.METHODS, order=args.E_FPR,
            data=res_common,
            linestyles=linestyles,
            markers=markers,
            palette=col_pal
        )
        fig_out = os.path.join(fig_dir, "shared_genes_by_fdr.pdf")
        sns_plot.figure.savefig(fig_out)
        plt.close()

    for design in search_params_dict['designs']:
        df_all_tmp = df_all_res.loc[df_all_res["DESIGN"] == design]
        df_merge = df_all_tmp.loc[df_all_tmp["METHOD"] == "epcy"]
        df_merge = df_merge.rename(columns={"VALUE": "epcy"})
        df_deg = df_all_tmp.loc[df_all_tmp["METHOD"] != "epcy"]
        deg_methods = pd.unique(df_deg.METHOD)

        for deg_m in deg_methods:
            df_deg_m = df_all_tmp.loc[df_all_tmp["METHOD"] == deg_m]
            df_deg_m = df_deg_m.loc[df_deg_m["TYPE"] == "PADJ"]
            df_deg_m = df_deg_m.rename(columns={"VALUE": deg_m})

            df_merge = pd.merge(df_merge, df_deg_m, on='ID', how='outer')

        csv_out = os.path.join(fig_dir, design + "_merge_table.xls")
        df_merge.to_csv(csv_out, index=False, sep="\t")

        methods_b = ["epcy_b"]
        df_merge["epcy_b"] = True
        cutoff = res_cutoff["epcy"][1]
        df_merge.loc[df_merge["epcy"] < cutoff, "epcy_b"] = False
        df_merge.loc[pd.isna(df_merge["epcy"]), "epcy_b"] = False

        if "deseq2" in search_params_dict["methods"]:
            df_merge["deseq2_b"] = True
            cutoff = res_cutoff["deseq2"][1]
            df_merge.loc[df_merge["deseq2"] < cutoff, "deseq2_b"] = False
            df_merge.loc[pd.isna(df_merge["deseq2"]), "deseq2_b"] = False
            methods_b.append("deseq2_b")

        if "edger" in search_params_dict["methods"]:
            df_merge["edger_b"] = True
            cutoff = res_cutoff["edger"][1]
            df_merge.loc[df_merge["edger"] < cutoff, "edger_b"] = False
            df_merge.loc[pd.isna(df_merge["edger"]), "edger_b"] = False
            methods_b.append("edger_b")

        if "voom" in search_params_dict["methods"]:
            df_merge["voom_b"] = True
            cutoff = res_cutoff["voom"][1]
            df_merge.loc[df_merge["voom"] < cutoff, "voom_b"] = False
            df_merge.loc[pd.isna(df_merge["voom"]), "voom_b"] = False
            methods_b.append("voom_b")

        if "mast" in search_params_dict["methods"]:
            df_merge["mast_b"] = True
            cutoff = res_cutoff["mast"][1]
            df_merge.loc[df_merge["mast"] < cutoff, "mast_b"] = False
            df_merge.loc[pd.isna(df_merge["mast"]), "mast_b"] = False
            methods_b.append("mast_b")

        if "trend" in search_params_dict["methods"]:
            df_merge["trend_b"] = True
            cutoff = res_cutoff["trend"][1]
            df_merge.loc[df_merge["trend"] < cutoff, "trend_b"] = False
            df_merge.loc[pd.isna(df_merge["trend"]), "trend_b"] = False
            methods_b.append("trend_b")

        res = df_merge.groupby(methods_b).size()
        res = res[1:]
        up.plot(res)
        fig_out = os.path.join(fig_dir, design + "_upset_1.pdf")
        plt.savefig(fig_out)
        plt.close()

        df_merge["epcy_b"] = True
        cutoff = res_cutoff["epcy"][0]
        df_merge.loc[df_merge["epcy"] < cutoff, "epcy_b"] = False
        df_merge.loc[pd.isna(df_merge["epcy"]), "epcy_b"] = False

        if "deseq2" in search_params_dict["methods"]:
            df_merge["deseq2_b"] = True
            cutoff = res_cutoff["deseq2"][0]
            df_merge.loc[df_merge["deseq2"] < cutoff, "deseq2_b"] = False
            df_merge.loc[pd.isna(df_merge["deseq2"]), "deseq2_b"] = False

        if "edger" in search_params_dict["methods"]:
            df_merge["edger_b"] = True
            cutoff = res_cutoff["edger"][0]
            df_merge.loc[df_merge["edger"] < cutoff, "edger_b"] = False
            df_merge.loc[pd.isna(df_merge["edger"]), "edger_b"] = False

        if "voom" in search_params_dict["methods"]:
            df_merge["voom_b"] = True
            cutoff = res_cutoff["voom"][0]
            df_merge.loc[df_merge["voom"] < cutoff, "voom_b"] = False
            df_merge.loc[pd.isna(df_merge["voom"]), "voom_b"] = False

        if "mast" in search_params_dict["methods"]:
            df_merge["mast_b"] = True
            cutoff = res_cutoff["mast"][0]
            df_merge.loc[df_merge["mast"] < cutoff, "mast_b"] = False
            df_merge.loc[pd.isna(df_merge["mast"]), "mast_b"] = False

        if "trend" in search_params_dict["methods"]:
            df_merge["trend_b"] = True
            cutoff = res_cutoff["trend"][0]
            df_merge.loc[df_merge["trend"] < cutoff, "trend_b"] = False
            df_merge.loc[pd.isna(df_merge["trend"]), "trend_b"] = False

        res = df_merge.groupby(methods_b).size()
        res = res[1:]
        up.plot(res)
        fig_out = os.path.join(fig_dir, design + "_upset_0.pdf")
        plt.savefig(fig_out)
        plt.close()

        res_top_tmp = res_top.loc[res_top["DESIGN"] == design]
        res_top_tmp = res_top_tmp.loc[res_top_tmp["num_gene"] == args.TOP]
        df_merge["epcy_b"] = True
        cutoff = res_top_tmp["epcy"].values[0]
        df_merge.loc[df_merge["epcy"] < cutoff, "epcy_b"] = False
        df_merge.loc[pd.isna(df_merge["epcy"]), "epcy_b"] = False

        if "deseq2" in search_params_dict["methods"]:
            df_merge["deseq2_b"] = True
            cutoff = res_top_tmp["deseq2"].values[0]
            df_merge.loc[df_merge["deseq2"] < cutoff, "deseq2_b"] = False
            df_merge.loc[pd.isna(df_merge["deseq2"]), "deseq2_b"] = False

        if "edger" in search_params_dict["methods"]:
            df_merge["edger_b"] = True
            cutoff = res_top_tmp["edger"].values[0]
            df_merge.loc[df_merge["edger"] < cutoff, "edger_b"] = False
            df_merge.loc[pd.isna(df_merge["edger"]), "edger_b"] = False

        if "voom" in search_params_dict["methods"]:
            df_merge["voom_b"] = True
            cutoff = res_top_tmp["voom"].values[0]
            df_merge.loc[df_merge["voom"] < cutoff, "voom_b"] = False
            df_merge.loc[pd.isna(df_merge["voom"]), "voom_b"] = False

        if "mast" in search_params_dict["methods"]:
            df_merge["mast_b"] = True
            cutoff = res_top_tmp["mast"].values[0]
            df_merge.loc[df_merge["mast"] < cutoff, "mast_b"] = False
            df_merge.loc[pd.isna(df_merge["mast"]), "mast_b"] = False

        if "trend" in search_params_dict["methods"]:
            df_merge["trend_b"] = True
            cutoff = res_top_tmp["trend"].values[0]
            df_merge.loc[df_merge["trend"] < cutoff, "trend_b"] = False
            df_merge.loc[pd.isna(df_merge["trend"]), "trend_b"] = False

        res = df_merge.groupby(methods_b).size()
        res = res[1:]
        up.plot(res)
        fig_out = os.path.join(fig_dir, design + "_upset_top_" + str(args.TOP) + ".pdf")
        plt.savefig(fig_out)
        plt.close()

        for method in search_params_dict["methods"]:
            if method != "epcy":
                sns_plot = sns.scatterplot(
                    x="epcy", y=method, size=2, linewidth=0,
                    data=df_merge, color="black"
                )
                sns_plot.axes.axvline(res_cutoff["epcy"][0], ls='--', color='black')
                sns_plot.axes.axvline(res_top_tmp["epcy"].values[0], ls=':', color='black')
                sns_plot.axes.axhline(res_cutoff[method][0], ls='--', color='black')
                sns_plot.axes.axhline(res_top_tmp[method].values[0], ls=':', color='black')

                if method == "deseq2":
                    col_dot = [col_pal[0]]
                if method == "mast":
                    col_dot = [col_pal[0]]
                if method == "edger":
                    col_dot = [col_pal[1]]
                if method == "voom":
                    col_dot = [col_pal[2]]
                if method == "trend":
                    col_dot = [col_pal[1]]

                sns_plot.fill_between(
                    df_merge["epcy"],
                    0, res_cutoff[method][0],
                    where=(df_merge["epcy"] > res_cutoff["epcy"][0]) & (df_merge["epcy"] < res_top_tmp["epcy"].values[0]),
                    color=col_pal[len(col_pal)-1], alpha=0.25)
                sns_plot.fill_between(
                    df_merge["epcy"],
                    0, res_cutoff[method][0],
                    where=df_merge["epcy"] > res_top_tmp["epcy"].values[0],
                    color=col_pal[len(col_pal)-1], alpha=0.5)
                sns_plot.fill_between(
                    df_merge["epcy"],
                    res_cutoff[method][0], res_top_tmp[method].values[0],
                    where=df_merge["epcy"] > res_top_tmp["epcy"].values[0],
                    color=col_pal[len(col_pal)-1], alpha=0.25)

                sns_plot.fill_between(
                    df_merge["epcy"],
                    res_cutoff[method][0], res_top_tmp[method].values[0],
                    where=df_merge["epcy"] < res_cutoff["epcy"][0],
                    color=col_dot, alpha=0.25)
                sns_plot.fill_between(
                    df_merge["epcy"],
                    res_top_tmp[method].values[0], df_merge[method].max(),
                    where=df_merge["epcy"] < res_cutoff["epcy"][0],
                    color=col_dot, alpha=0.5)
                sns_plot.fill_between(
                    df_merge["epcy"],
                    res_top_tmp[method].values[0], df_merge[method].max(),
                    where=(df_merge["epcy"] > res_cutoff["epcy"][0]) & (df_merge["epcy"] < res_top_tmp["epcy"].values[0]),
                    color=col_dot, alpha=0.25)

                for line in range(0, df_merge.shape[0]):
                    if df_merge.epcy[line] > res_top_tmp["epcy"].values[0] or df_merge[method][line] > res_top_tmp[method].values[0]:
                        sns_plot.text(
                            df_merge.epcy[line]+0.02, df_merge[method][line],
                            df_merge.ID[line], horizontalalignment='left',
                            size='medium', color='black', weight='semibold')

                    #if df_merge.epcy[line] < 0.2 and df_merge[method][line] > res_cutoff[method][0]:
                    #    sns_plot.text(
                    #        df_merge.epcy[line]+0.02, df_merge[method][line],
                    #        df_merge.ID[line], horizontalalignment='left',
                    #        size='medium', color='black', weight='semibold')

                    #if df_merge.epcy[line] > 0.45 and df_merge[method][line] < res_cutoff[method][0]:
                    #    sns_plot.text(
                    #        df_merge.epcy[line]+0.02, df_merge[method][line],
                    #        df_merge.ID[line], horizontalalignment='left',
                    #        size='medium', color='black', weight='semibold')

                fig_out = os.path.join(fig_dir, design + "_epcy_vs_" + method + "_top_" + str(args.TOP) + ".pdf")
                sns_plot.figure.savefig(fig_out)
                plt.close()
