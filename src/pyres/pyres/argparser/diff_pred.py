from .common import *


def get_argparser_diff_pred(parser):
    parser.add_argument("--outdir",
                        dest="OUTDIR",
                        help="Output directory",
                        type=str,
                        default='OUT')

    parser.add_argument("--size",
                        dest="SIZE",
                        help="Radius of a dot, in points. (Default:5.0)",
                        type=float,
                        default=5.0)

    parser.add_argument("--strip",
                        dest="STRIP",
                        help='To create a strip plot.',
                        action='store_true')

    parser.set_defaults(STRIP=False)
