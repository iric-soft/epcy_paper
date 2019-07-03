import argparse

from .argparser.clust_exp import *
from .argparser.gen_cv import *
from .argparser.eval_cv import *
from .argparser.heatmap_cv import *

from .tools.clust_exp import main_clust_exp
from .tools.gen_cv import main_gen_cv
from .tools.eval_cv import main_eval_cv
from .tools.heatmap_cv import main_heatmap_cv



# ###########################################################################
# Main function
def main():

    argparser = argparse.ArgumentParser(prog='PROG')
    subparsers = argparser.add_subparsers(help='sub-command help')

    # create the argparser for the "clust_exp" command
    clust_exp = subparsers.add_parser(
        'clust_exp',
        help='Clustring compute on expression data.'
    )
    clust_exp.set_defaults(func=main_clust_exp)
    get_argparser_compare_deg(clust_exp)

    # create the argparser for the "gen_cv" command
    gen_cv = subparsers.add_parser(
        'gen_cv',
        help='Generate folder and project file for cross validation'
    )
    gen_cv.set_defaults(func=main_gen_cv)
    get_argparser_gen_cv(gen_cv)

    # create the argparser for the "eval_cv" command
    eval_cv = subparsers.add_parser(
        'eval_cv',
        help='Evalutate performance of kt, LDE (Limma + EdgeR + DEseq)'
    )
    eval_cv.set_defaults(func=main_eval_cv)
    get_argparser_eval_cv(eval_cv)

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
