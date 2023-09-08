import argparse

from .argparser.clust_exp import *
from .argparser.clust_exp_umap import *
from .argparser.clust_all_umap import *
from .argparser.gen_cv import *
from .argparser.eval_cv import *
from .argparser.heatmap_cv import *
from .argparser.density import *
from .argparser.gene_density import *
from .argparser.log_reg import *
from .argparser.eval_tt import *
from .argparser.eval_auc import *
from .argparser.eval_random import *
from .argparser.diff_pred import *
from .argparser.eval_ss import *

from .tools.clust_exp import main_clust_exp
from .tools.clust_exp_umap import main_clust_exp_umap
from .tools.clust_all_umap import main_clust_all_umap
from .tools.gen_cv import main_gen_cv
from .tools.eval_cv import main_eval_cv
from .tools.heatmap_cv import main_heatmap_cv
from .tools.density import main_density
from .tools.gene_density import main_gene_density
from .tools.log_reg import main_log_reg
from .tools.eval_tt import main_eval_tt
from .tools.eval_auc import main_eval_auc
from .tools.eval_random import main_eval_random
from .tools.diff_pred import main_diff_pred
from .tools.eval_ss import main_eval_ss



# ###########################################################################
# Main function
def main():

    argparser = argparse.ArgumentParser(prog='PROG')
    subparsers = argparser.add_subparsers(help='sub-command help')

    # create the argparser for the "diff_pred" command
    diff_pred = subparsers.add_parser(
        'diff_pred',
        help='Display plot to compare diff gene vs pred gene.'
    )
    diff_pred.set_defaults(func=main_diff_pred)
    get_argparser_diff_pred(diff_pred)

    # create the argparser for the "eval_ss" command
    eval_ss = subparsers.add_parser(
        'eval_ss',
        help='Evaluate subsampling results.'
    )
    eval_ss.set_defaults(func=main_eval_ss)
    get_argparser_eval_ss(eval_ss)

    # create the argparser for the "density" command
    density = subparsers.add_parser(
        'density',
        help='Display a density plot of top x features.'
    )
    density.set_defaults(func=main_density)
    get_argparser_density(density)

    # create the argparser for the "gene_density" command
    gene_density = subparsers.add_parser(
        'gene_density',
        help='Display a density plot of one gene.'
    )
    gene_density.set_defaults(func=main_gene_density)
    get_argparser_gene_density(gene_density)


    # create the argparser for the "clust_exp" command
    log_reg = subparsers.add_parser(
        'log_reg',
        help='Display a logistic regression plot build on a feature.'
    )
    log_reg.set_defaults(func=main_log_reg)
    get_argparser_log_reg(log_reg)

    # create the argparser for the "clust_exp" command
    clust_exp = subparsers.add_parser(
        'clust_exp',
        help='Clustring compute on expression data.'
    )
    clust_exp.set_defaults(func=main_clust_exp)
    get_argparser_clust_exp(clust_exp)

    # create the argparser for the "clust_exp_umap" command
    clust_exp_umap = subparsers.add_parser(
        'clust_exp_umap',
        help='UMAP clustring compute on expression data.'
    )
    clust_exp_umap.set_defaults(func=main_clust_exp_umap)
    get_argparser_clust_exp_umap(clust_exp_umap)

    # create the argparser for the "clust_all_umap" command
    clust_all_umap = subparsers.add_parser(
        'clust_all_umap',
        help='UMAP clustring using all design together.'
    )
    clust_all_umap.set_defaults(func=main_clust_all_umap)
    get_argparser_clust_all_umap(clust_all_umap)

    # create the argparser for the "gen_cv" command
    gen_cv = subparsers.add_parser(
        'gen_cv',
        help='Generate folder and project file for cross validation'
    )
    gen_cv.set_defaults(func=main_gen_cv)
    get_argparser_gen_cv(gen_cv)

    # create the argparser for the "eval_random" command
    eval_random = subparsers.add_parser(
        'eval_random',
        help='Evalutate performance of kt, LDE (Limma + EdgeR + DEseq) on random design'
    )
    eval_random.set_defaults(func=main_eval_random)
    get_argparser_eval_random(eval_random)

    # create the argparser for the "eval_auc" command
    eval_auc = subparsers.add_parser(
        'eval_auc',
        help='Evalutate AUC performance of kt, LDE (Limma + EdgeR + DEseq)'
    )
    eval_auc.set_defaults(func=main_eval_auc)
    get_argparser_eval_auc(eval_auc)

    # create the argparser for the "eval_cv" command
    eval_cv = subparsers.add_parser(
        'eval_cv',
        help='Evalutate performance of kt, LDE (Limma + EdgeR + DEseq)'
    )
    eval_cv.set_defaults(func=main_eval_cv)
    get_argparser_eval_cv(eval_cv)

    # create the argparser for the "eval_tt" command
    eval_tt = subparsers.add_parser(
        'eval_tt',
        help='Evalutate performance of kt, LDE (Limma + EdgeR + DEseq)'
    )
    eval_tt.set_defaults(func=main_eval_tt)
    get_argparser_eval_cv(eval_tt)


    # create the argparser for the "heatmap_cv" command
    heatmap_cv = subparsers.add_parser(
        'heatmap_cv',
        help='display heatmap of top x features a cross all datasets of cross validation.'
    )
    heatmap_cv.set_defaults(func=main_heatmap_cv)
    get_argparser_heatmap_cv(heatmap_cv)

    # recover arguments
    args = argparser.parse_args()

    # execute the command
    args.func(args, argparser)
