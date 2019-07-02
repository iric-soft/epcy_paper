from .common import *

def get_argparser_eval_cv(parser):

    parser.add_argument("-f",
                        dest="FOLD",
                        help="fold value (Default: 3).",
                        type=int,
                        default=3)

    parser.add_argument("-m",
                        dest="MATRIX",
                        help="path matrix file.",
                        type=str,
                        default=None)

    parser.add_argument("-p",
                        dest="PATH",
                        help="path to data folder ",
                        type=lambda x: is_valid_path(parser, x))

    parser.add_argument("--cpm",
                        dest="CPM",
                        help="To normalize the matrix, as Count Par Million (CPM)",
                        action='store_true')
                        
    parser.add_argument("--loo",
                        dest="LOO",
                        help="Activate leave one out (Default: False)",
                        action='store_true')

    parser.add_argument("--biotype",
                        dest="BIOTYPE",
                        help="List of filtred biotype (ex: protein_coding,antisense,pseudogene).",
                        type=str,
                        default=None)

    parser.add_argument("--bf",
                        dest="BF",
                        help="path to biotype annotation file.",
                        type=str)

    parser.add_argument("--design",
                        dest="DESIGN",
                        help="Name of designs to evaluate",
                        nargs='+' ,
                        type=str,
                        default=None)

    parser.add_argument("--lfc",
                        dest="LOG_FC",
                        help="abs(LOG_FC) filter value (Default: 0).",
                        type=float,
                        default=0)

    parser.add_argument("--mcc",
                        dest="MCC",
                        help="MCC filter value(Default: 0).",
                        type=float,
                        default=0.0)

    parser.add_argument("--methods",
                        dest = "METHODS",
                        help = 'methods to be tested',
                        type = str,
                        nargs = '+',
                        default = ['deseq2', 'edger', 'limma', 'epcy'])

    parser.add_argument("--n_datasets",
                        dest="N_DATASETS",
                        help="Number of datasets to use, DEBUG ONLY!",
                        type=str,
                        default='all')

    parser.add_argument("--outdir",
                        dest="OUTDIR",
                        help="Output directory for evaluation metrics tables and log file",
                        type=str,
                        default='OUT')

    parser.add_argument("--pvalue",
                        dest="PVALUE",
                        help="PValue filter value(Default: 0.05).",
                        type=float,
                        default=0.05)

    parser.add_argument("--query",
                        dest="QUERY",
                        help="Query value in the class column (Default: Query)",
                        type=str,
                        default="Query")

    parser.add_argument("--shuffle_seeds",
                        dest = "SHUFFLE_SEEDS",
                        help = 'set seeds for dataset rankings shuffles before predictions',
                        nargs = '+',
                        type = int,
                        default = None)

    parser.add_argument("--subgroup",
                        dest="SUBGROUP",
                        help="Name column of group sample in the design file (Default: subgroup)",
                        type=str,
                        default="subgroup")

    parser.add_argument("--top_values",
                        dest = "TOP_VALUES",
                        type = int,
                        help = 'Explicitly set values of top to search',
                        nargs= '+',
                        default = [1, 3, 5, 10, 50, 100, 200])

    parser.add_argument("--use_LR",
                        dest = "USE_LR",
                        help = 'Use a logistic regression (without regulariser) classifier for prediction',
                        action = 'store_true')

    parser.add_argument("--use_randF",
                        dest = "USE_RANDF",
                        help = 'Use a random forest classifier for prediction',
                        action = 'store_true')


    parser.set_defaults(LOO=False)
    parser.set_defaults(USE_RANDF=False)
    parser.set_defaults(USE_LR=False)
