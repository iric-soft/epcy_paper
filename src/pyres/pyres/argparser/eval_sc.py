from .common import *

def get_argparser_eval_ss(parser):

    parser.add_argument("-p",
                        dest="PATH",
                        help="path to data folder ",
                        type=lambda x: is_valid_path(parser, x))

    parser.add_argument("-q",
                        dest="QUANT",
                        help="Software used for quantification (STAR, htseq). (Default: STAR)",
                        type=str,
                        default='STAR_RSEM')

    parser.add_argument("-t",
                        dest="TYPE_QUANT",
                        help="Type of quantification (readcounts, tpm). (Default: readcounts)",
                        type=str,
                        default='readcounts')

    parser.add_argument("--lfc",
                        dest="LOG_FC",
                        help="abs(LOG_FC) filter value (Default: 0 (disable)).",
                        type=float,
                        default=0)

    parser.add_argument("--mcc",
                        dest="MCC",
                        help="MCC filter value(Default: -2 (disable)).",
                        type=float,
                        default=-2)

    parser.add_argument("--methods",
                        dest="METHODS",
                        help='methods to be tested',
                        type=str,
                        nargs='+',
                        default=['deseq2', 'edger', 'limma', 'epcy'])

    parser.add_argument("--reps",
                        dest="REPS",
                        help='Replicates to be tested',
                        type=str,
                        nargs='+',
                        default=['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20'])

    parser.add_argument("--numsamples",
                        dest="NUMSAMPLES",
                        help='Mean number of samples by replicates',
                        type=str,
                        nargs='+',
                        default=['5000', '8000', '10000'])

    parser.add_argument("--outdir",
                        dest="OUTDIR",
                        help="Output directory for evaluation metrics tables and log file",
                        type=str,
                        default='OUT')

    parser.add_argument("--ids",
                        dest="IDS",
                        help="List genes ids to display",
                        type=str,
                        nargs='+',
                        default=[])

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
