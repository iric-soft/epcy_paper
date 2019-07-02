import argparse

from .argparser.compare_deg import *
from .argparser.compare_mcc import *
from .argparser.gen_cv import *
from .argparser.eval_cv import *

from .tools.compare_deg import main_compare_deg
from .tools.compare_mcc import main_compare_mcc
from .tools.gen_cv import main_gen_cv
from .tools.eval_cv import main_eval_cv



# ###########################################################################
# Main function
def main():

    argparser = argparse.ArgumentParser(prog='PROG')
    subparsers = argparser.add_subparsers(help='sub-command help')

    # create the argparser for the "deg" command
    deg = subparsers.add_parser(
        'compare_deg',
        help='Compare with deg analysis.'
    )
    deg.set_defaults(func=main_compare_deg)
    get_argparser_compare_deg(deg)

    # create the argparser for the "mcc" command
    mcc = subparsers.add_parser(
        'compare_mcc',
        help='Compare with mcc.'
    )
    mcc.set_defaults(func=main_compare_mcc)
    get_argparser_compare_mcc(mcc)

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

    # recover arguments
    args = argparser.parse_args()

    # execute the command
    args.func(args, argparser)
