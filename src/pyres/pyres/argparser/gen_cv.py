from .common import *

def get_argparser_gen_cv(parser):
    parser.add_argument("-p",
                        dest="PATH",
                        help="path to the foder which have design.tsv file",
                        type=lambda x: is_valid_path(parser, x))

    parser.add_argument("--loo",
                        dest="LOO",
                        help="Activate leave one out (Default: False)",
                        action='store_true')

    parser.add_argument("-f",
                        dest="FOLD",
                        help="fold value (Default: 3).",
                        type=int,
                        default=3)

    parser.add_argument("--subgroup",
                        dest="SUBGROUP",
                        help="Name column of group sample in the design file (Default: group)",
                        type=str,
                        default="subgroup")

    parser.add_argument("--query",
                        dest="QUERY",
                        help="Query value in the class column (Default: Query)",
                        type=str,
                        default="Query")

    parser.set_defaults(LOO=False)
