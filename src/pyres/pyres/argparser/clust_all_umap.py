from .common import *

def get_argparser_clust_all_umap(parser):

    parser.add_argument("-m",
                        dest="MATRIX",
                        help="path to matrix file.",
                        type=str,
                        default=None)

    parser.add_argument("-p",
                        dest="PATH",
                        help="Folder path of all input files",
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

    parser.add_argument("--cpm",
                        dest="CPM",
                        help="To normalize the matrix, as Count Par Million (CPM)",
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

    parser.add_argument("--ext",
                        dest="EXT",
                        help="Extension to select output format.",
                        type=str,
                        default='pdf')

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
                        help="pValue filter value(Default: 0.05).",
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

    parser.add_argument("--scaled",
                        dest="SCALED",
                        help="Scale the fig, to read each samples and 'features'.",
                        action='store_true')

    parser.set_defaults(SCALED=False)
