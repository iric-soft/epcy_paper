from .common import *

def get_argparser_density(parser):

    parser.add_argument("-m",
                        dest="MATRIX",
                        help="Directoy path to matrix file.",
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
                        help="(Not available for density tools)",
                        type=str,
                        default=None)

    parser.add_argument("--bf",
                        dest="BF",
                        help="path to biotype annotation file.",
                        type=str)

    parser.add_argument("--cpm",
                        dest="CPM",
                        help="To normalize the matrix, as Count Par Million (CPM)",
                        action='store_true')

    parser.add_argument("--design",
                        dest="DESIGN",
                        help="Name of design used",
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

    parser.add_argument("--subgroup",
                        dest="SUBGROUP",
                        help="Name column of group sample in the design file (Default: subgroup)",
                        type=str,
                        default="subgroup")

    parser.add_argument("--top",
                        dest = "TOP",
                        type = int,
                        help = 'Top x for each methods',
                        default = 3)

    parser.set_defaults(BIOTYPE=None)
