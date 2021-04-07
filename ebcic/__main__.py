import argparse
import ebcic
# ebcic. in ebcic.function can be removed by using below:
# from ebcic import print_interval, exact, Params, __version__
# from . import print_interval, exact, Params, __version__


def main():

    parser = argparse.ArgumentParser(
        description=(
            "Exact Binomial Interval Calculator, "
            f"ver. {ebcic.__version__}"),
        # For __main__.py, c.f. https://bugs.python.org/issue22240
        prog=(
            None if globals().get('__spec__') is None else
            'python -m {}'.format(__spec__.name.partition('.')[0])),
    )
    parser.add_argument(
        "-k",
        "--errors",
        type=int,
        required=True,
        help="number of errors",
    )
    parser.add_argument(
        "-n",
        "--trials",
        type=int,
        required=True,
        help="number of trials",
    )
    parser.add_argument(
        "-c",
        "--confi-perc",
        type=float,
        help=(
            "Confidence percentage for two-sided [0-100]. "
            "For one-sided, set '--confi-perc <cp_t>' where"
            " cp_t ="
            "   2 * (confidence percentage for one-sided [50-100]) - 100."
            )
    )
    parser.add_argument(
        "-a",
        "--alpha",
        type=float,
        help="ALPHA = 1 - CONFI_PERC/100"
    )
    parser.add_argument(
        "-u",
        "--upper",
        action='store_true',
        help="print upper interval of given confidence percentage "
             "for two-sided"
    )
    parser.add_argument(
        "-l",
        "--lower",
        action='store_true',
        help="print lower interval of given confidence percentage "
             "for two-sided"
    )
    args = parser.parse_args()

    # args check
    if (
            ((args.alpha is not None) and (args.confi_perc is not None)) or
            ((args.alpha is None) and (args.confi_perc is None))):
        parser.error(
            "either --alpha (-a) or --confi-perc (-c) shall be set!!")

    if (args.confi_perc is not None):
        params = ebcic.Params(
            k=args.errors,          # Number of errors
            n=args.trials,          # number of trials
            confi_perc=args.confi_perc
        )
    elif (args.alpha is not None):
        params = ebcic.Params(
            k=args.errors,          # Number of errors
            n=args.trials,          # number of trials
            alpha=args.alpha
        )

    # body
    if args.lower or args.upper:
        interval = ebcic.exact(params)
        if args.lower:
            print(interval[0])
        if args.upper:
            print(interval[1])
    else:
        # more info
        ebcic.print_interval(params)


if __name__ == "__main__":
    main()
