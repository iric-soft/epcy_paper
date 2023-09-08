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
                        help="abs(LOG_FC) filter value (Default: 0 (disable)).",
                        type=float,
                        default=0)

    parser.add_argument("--size",
                        dest="SIZE",
                        help="Dot size (Default: 3).",
                        type=float,
                        default=3)

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

    parser.add_argument("--cohorts",
                        dest="COHORTS",
                        help='cohorts used for mean',
                        type=str,
                        nargs='+',
                        default=['3_vs_3', '4_vs_28', '6_vs_65', '17_vs_314', '27_vs_563'])

    parser.add_argument("--outdir",
                        dest="OUTDIR",
                        help="Output directory for evaluation metrics tables and log file",
                        type=str,
                        default='OUT')

    parser.add_argument("--design",
                        dest="DESIGN",
                        help="Name of designs to analyse.",
                        type=str,
                        default="30_t15_17")

    parser.add_argument("--ids",
                        dest="IDS",
                        help="List genes ids to display",
                        type=str,
                        nargs='+',
                        default=[])

    parser.add_argument("--p_ss",
                        dest="P_SS",
                        help="% used for each subssampling",
                        type=float,
                        nargs='+',
                        default=[0, 0.01, 0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])

    parser.add_argument("--reps",
                        dest="REPS",
                        help="Select replicates you want use",
                        type=int,
                        nargs='+',
                        default=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10])

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
