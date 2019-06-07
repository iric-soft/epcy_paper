from .common import *

def get_argparser_eval_cv(parser):
    parser.add_argument("-p",
                        dest="PATH",
                        help="path to the fold foder ",
                        type=lambda x: is_valid_path(parser, x))

    parser.add_argument("--extfc",
                        dest="EXTFC",
                        help="activate corrected foldchange available in subfolder c.",
                        action='store_true')

    parser.add_argument("--loo",
                        dest="LOO",
                        help="Activate leave one out (Default: False)",
                        action='store_true')

    parser.add_argument("-f",
                        dest="FOLD",
                        help="fold value (Default: 3).",
                        type=int,
                        default=3)

    parser.add_argument("--group",
                        dest="GROUP",
                        help="Name column of group sample in the design file (Default: group)",
                        type=str,
                        default="group")

    parser.add_argument("--query",
                        dest="QUERY",
                        help="Query value in the class column (Default: Query)",
                        type=str,
                        default="Query")

    parser.add_argument("--gene",
                        dest="GENE",
                        help="Activate gene level (Default: False)",
                        action='store_true')

    parser.add_argument("--biotype",
                        dest="BIOTYPE",
                        help="List of filtred biotype (ex: protein_coding,antisense,pseudogene).",
                        type=str,
                        default=None)

    parser.add_argument("--lfc",
                        dest="LOG_FC",
                        help="abs(LOG_FC) filter value (Default: 0.3).",
                        type=float,
                        default=0.3)

    parser.add_argument("--mcc",
                        dest="MCC",
                        help="MCC filter value(Default: 0).",
                        type=float,
                        default=0.0)

    parser.add_argument("--pvalue",
                        dest="PVALUE",
                        help="pValue filter value(Default: 0.05).",
                        type=float,
                        default=0.05)
    #### added output dir (LS)
    parser.add_argument("--outdir",
                        dest="OUTDIR",
                        help="output directory for evaluation metrics tables and log file",
                        type=str,
                        default='OUT'
        )
    #### added argument subgroup for control
    parser.add_argument("--subgroups",
                        dest="SUBGROUPS",
                        help="name of subgroup to evaluate",
                        nargs='+' ,
                        type=str,
                        default=None
        )

    parser.add_argument("--n_datasets",
                        dest="N_DATASETS",
                        help="number of datasets to use, DEBUG ONLY!",
                        type=str,
                        default='all'
        )

    parser.add_argument("--event_log",
                        dest= "EVENT_LOG",
                        help= "puts logs into a file log with date. Else it will print to terminal.",
                        action = "store_true"
        )

    parser.add_argument("--top_range",
                        dest = 'TOP_RANGE',
                        help = '[MIN, MAX, STEP]',
                        nargs= 3,
                        type = int,
                        default= None
        )

    parser.add_argument("--top_values",
                        dest = "TOP_VALUES",
                        type = int,
                        help = 'explicitly set values of top to search',
                        nargs= '+',
                        default = [1, 3, 5, 10, 50, 100, 200]
        )
    parser.add_argument("--shuffle_seeds",
                        dest = "SHUFFLE_SEEDS",
                        help = 'set seeds for dataset rankings shuffles before predictions',
                        nargs = '+',
                        type = int,
                        default = None
        )

    parser.add_argument("--use_randF",
                        dest = "USE_RANDF",
                        help = 'set true: will use a random forest classifier for prediction',
                        action = 'store_true'
        )

    parser.add_argument("--use_logR",
                        dest = "USE_LOGR",
                        help = 'set true: will use a logistic regression (without regulariser) classifier for prediction',
                        action = 'store_true'
        )
    parser.add_argument("--methods",
                        dest = "METHODS",
                        help = 'methods to be tested',
                        type = str, 
                        nargs = '+',
                        default = ['deseq2', 'edger_qlf', 'limma', 'kt']
        )
    parser.set_defaults(LOO=False)
    parser.set_defaults(GENE=False)
