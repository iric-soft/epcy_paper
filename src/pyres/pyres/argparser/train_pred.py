from .common import *

def get_argparser_train_pred(parser):
    parser.add_argument("-p",
                        dest="PATH",
                        help="path to the fold foder ",
                        type=lambda x: is_valid_path(parser, x))

    parser.add_argument("-e",
                        dest="EXP",
                        help="Path folder of expression of dataset which need to be predict",
                        type=lambda x: is_valid_file(parser, x))

    parser.add_argument("-d",
                        dest="DESIGN",
                        help="Design of dataset which need to be predict",
                        type=lambda x: is_valid_path(parser, x))

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

    parser.set_defaults(GENE=False)
