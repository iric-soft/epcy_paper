from .common import *

def get_argparser_compare_mcc(parser):
    parser.add_argument("-p",
                        dest="PATH",
                        help="Folder path of all input files",
                        type=lambda x: is_valid_path(parser, x))

    parser.add_argument("--prefix",
                        dest="PREFIX",
                        help="prefix for output files",
                        type=str,
                        default="")

    parser.add_argument("--dataset",
                        dest="DATASET",
                        help="list of dataset separet by ',' (Default: 28_inv16,28_inv16_vs28,33_MLL,18_t8_21)",
                        type=str,
                        default="28_inv16,28_inv16_vs28,33_MLL,18_t8_21")

    parser.add_argument("--gene",
                        dest="GENE",
                        help="Calculate AUCs on genes (vs Transcripts)",
                        action='store_true')

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

    parser.add_argument("--scaled",
                        dest="SCALED",
                        help="Scale the fig, to read each samples and 'features'.",
                        action='store_true')

    parser.add_argument("--biotype",
                        dest="BIOTYPE",
                        help="List of filtred biotype (ex: protein_coding,antisense,pseudogene).",
                        type=str,
                        default=None)
