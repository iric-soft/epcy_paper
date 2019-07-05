from .common import *

def get_argparser_heatmap_cv(parser):

    parser.add_argument("-f",
                        dest="FOLD",
                        help="fold value (Default: 3).",
                        type=int,
                        default=3)

    parser.add_argument("-p",
                        dest="PATH",
                        help="path to data folder ",
                        type=lambda x: is_valid_path(parser, x))

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

    parser.add_argument("--loo",
                        dest="LOO",
                        help="Activate leave one out (Default: False)",
                        action='store_true')

    parser.add_argument("--mcc",
                        dest="MCC",
                        help="MCC filter value(Default: 0).",
                        type=float,
                        default=0.0)

    parser.add_argument("--method",
                        dest = "METHOD",
                        help = 'Method analyzed: epcy, deseq2, edger or limma (Default: epcy)',
                        type = str,
                        default = 'epcy')

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

    parser.add_argument("--top",
                        dest = "TOP",
                        type = int,
                        help = 'Number of top x features used. (Default: 10)',
                        default = 10)

    parser.set_defaults(LOO=False)
