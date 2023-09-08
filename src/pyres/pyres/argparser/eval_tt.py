from .common import *

def get_argparser_eval_tt(parser):

    parser.add_argument("-m",
                        dest="MATRIX",
                        help="path to matrix file.",
                        type=str,
                        default=None)

    parser.add_argument("-p",
                        dest="PATH",
                        help="path to data folder ",
                        type=lambda x: is_valid_path(parser, x))

    parser.add_argument("-q",
                        dest="QUANT",
                        help="Software used for quantification (STAR, htseq). (Default: STAR)",
                        type=str,
                        default='STAR')

    parser.add_argument("-t",
                        dest="TYPE_QUANT",
                        help="Type of quantification (readcounts, tpm). (Default: readcounts)",
                        type=str,
                        default='readcounts')

    parser.add_argument("--biotype",
                        dest="BIOTYPE",
                        help="List of filtred biotype (ex: protein_coding,antisense,pseudogene).",
                        type=str,
                        default=None)

    parser.add_argument("--bf",
                        dest="BF",
                        help="path to biotype annotation file.",
                        type=str)

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

    parser.add_argument("--outdir",
                        dest="OUTDIR",
                        help="Output directory for evaluation metrics tables and log file",
                        type=str,
                        default='OUT')

    parser.add_argument("--pvalues",
                        dest="PVALUES",
                        help="PValues cutt-off used for LDE method (Default: [0.05]).",
                        type=float,
                        nargs= '+',
                        default=[0.05])

    parser.add_argument("--query",
                        dest="QUERY",
                        help="Query value in the class column (Default: Query)",
                        type=str,
                        default="Query")

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
                        default = [2, 3, 5, 10, 50, 100, 200])

    parser.set_defaults(USE_LR=False)
