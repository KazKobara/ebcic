# %load ./ebcic/ebcic.py
# %%writefile ./ebcic/ebcic.py
# To be synchronized with ebcic.py copy either of the above magic functions
# depending on whether to overwrite or to load, and then run this cell.
"""Exact Binomial Confidence Interval Calculator

This module provides functions to calculate exact binomial intervals
and some of the approximated intervals, such as normal approximation,
generalized rule of three, Wilson score intervals and so on.

This module also provides functions to depict their intervals.


When you use or publish the confidence interval obtained with the software,
please **refer to the program name, version, platform**, and so on, so that
readers can verify the correctness and reproducibility of the interval
with the parameters, n, k and the confidence level.

Example of the reference is:

**"This confidence interval is obtained by EBCIC x.x.x on Python 3."**

where X.X.X is the version of EBCIC.

The initial software is based on results obtained from a project,
JPNP16007, commissioned by
the New Energy and Industrial Technology Development Organization (NEDO).

THE SOFTWARE COMES WITH ABSOLUTELY NO WARRANTY.

Copyright (C) 2020-2022 National Institute of Advanced Industrial Science
and Technology (AIST).
"""
__version__ = '0.0.4'
import os
import sys
import logging
import numpy as np
import math
import warnings
from scipy import optimize
from scipy.stats import norm, binom, poisson, beta
from matplotlib import pyplot as plt
from matplotlib import cm as mplcm
from matplotlib import colors as colors
from decimal import Decimal, ROUND_HALF_UP
from platform import python_version
from packaging import version
if version.parse(python_version()) >= version.parse("3.8"):
    from statistics import NormalDist
    USE_NORMAL_DIST = True
'''
import scipy
import matplotlib

print(f"scipy     : {scipy.version.full_version}")
print(f"matplotlib: {matplotlib.__version__}")
print(f"python    : {python_version()}")  # os, math, warnings

print(f"np        : {np.__version__}")
print(f"sys       : {sys.version[:6]}")
print(f"logging   : {logging.__version__}")
'''

warnings.filterwarnings(
        'ignore', 'The iteration is not making good progress')

# Log handler's loglevels
# Choose from 'DEBUG', 'INFO', 'WARNING' 'ERROR', 'CRITICAL'.
loglevel_stderr = 'INFO'  # for stderr (stdout in jupyter note)
loglevel_file = 'DEBUG'   # for log file

# logger.setLevel should be DEBUG to handle all the levels at handlers,
# change the above handler's loglevels instead of this logger's.
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

# Check loglevels
numeric_level_stderr = getattr(logging, loglevel_stderr.upper(), None)
if not isinstance(numeric_level_stderr, int):
    raise ValueError('Invalid log level: %s' % loglevel_stderr)
numeric_level_file = getattr(logging, loglevel_file.upper(), None)
if not isinstance(numeric_level_file, int):
    raise ValueError('Invalid log level: %s' % loglevel_file)

# Create handlers with loglevels

# For Jupyter, existing 'sh' and 'fh' shall be removed first,
# or restart Jupyter kernel.
if 'sh' in globals():
    # print(sh)
    logger.removeHandler(sh)  # noqa: F821  # needed for Jupyter
if 'fh' in globals():
    # print(fh)
    logger.removeHandler(fh)  # noqa: F821  # needed for Jupyter

sh = logging.StreamHandler(sys.stdout)  # for Jupyter
# sh = logging.StreamHandler()  # for non Jupyter
sh.setLevel(numeric_level_stderr)
# sh.setLevel(logger.INFO)
fh = logging.FileHandler(
    mode='w',
    # filename=os.path.splitext(os.path.basename(__file__))[0] + ".log")
    filename=os.path.splitext(os.path.basename('ebcic'))[0] + ".log")
fh.setLevel(numeric_level_file)
# fh.setLevel(logger.DEBUG)
# Create formatter
formatter = logging.Formatter(
    # WARNI CRITI
    # "%(asctime)s %(levelname)5s %(lineno)d %(funcName)s %(message)s")
    "%(levelname)5.5s %(lineno)4d %(funcName)12.12s %(message)s")
formatter_simple = logging.Formatter(
    "%(levelname)5s: %(message)s")
# Cdd formatter to channels
sh.setFormatter(formatter_simple)
fh.setFormatter(formatter)
# Cdd channels to logger
logger.addHandler(sh)
logger.addHandler(fh)


# Definitions
def cli(cli_args=None):
    """EBCIC Command Line Interface

    See the command line examples on https://github.com/KazKobara/ebcic
    and below.

    Args:
        cli_args: Command line args.
            If 'cli_args is None', it uses args given by the command interface.

    Examples:
        # Show help
        # python -m ebcic -h
        >>> cli(['-h'])  # doctest: +SKIP

        ### Two-sided confidence interval for 0<k<n ###
        # python -m ebcic -k 1 -n 100 -c 95 -lu
        >>> cli(['-k', '1', '-n', '100', '-c', '95', '-lu'])
        0.00025314603297743815
        0.0544593853920806

        # python -m ebcic -k 1 -n 100 -r 2.5 -s 2.5 -lu
        >>> cli(['-k', '1', '-n', '100', '-r', '2.5', '-s', '2.5', '-lu'])
        0.00025314603297743815
        0.0544593853920806

        ### One-sided upper confidence interval for k=0 ###
        # python -m ebcic -k 0 -n 100 -c 95 -u
        >>> cli(['-k', '0', '-n', '100', '-c', '95', '-u'])
        0.029513049607039932

        # python -m ebcic -k 0 -n 100 --rej-perc-lower 5 -u
        >>> cli(['-k', '0', '-n', '100', '-r', '5', '-u'])
        0.029513049607039932

        ### One-sided lower confidence interval for k=n ###
        # python -m ebcic -k 100 -n 100 --rej-perc-upper 5 -l
        >>> cli(['-k', '100', '-n', '100', '-s', '5', '-l'])
        0.9704869503929601

        ### One-sided upper confidence interval for 0<k<n ###
        # python -m ebcic -k 1 -n 100 -r 5 -u
        >>> cli(['-k', '1', '-n', '100', '-r', '5.', '-u'])
        0.04655981145353891

        # python -m ebcic -k 1 -n 100 -c 90 -u
        >>> cli(['-k', '1', '-n', '100', '-c', '90.', '-u'])
        0.04655981145353891

        # python -m ebcic -k 1 -n 100 -a 0.1 -u
        >>> cli(['-k', '1', '-n', '100', '-a', '0.1', '-u'])
        0.04655981145353891

        ### One-sided lower confidence interval for 0<k<n ###
        # python -m ebcic -k 99 -n 100 -s 5 -l
        >>> cli(['-k', '99', '-n', '100', '-s', '5.', '-l'])
        0.9534401885464611

        # python -m ebcic -k 99 -n 100 -c 90 -l
        >>> cli(['-k', '99', '-n', '100', '-c', '90.', '-l'])
        0.9534401885464611

        # python -m ebcic -k 99 -n 100 -a 0.1 -l
        >>> cli(['-k', '99', '-n', '100', '-a', '0.1', '-l'])
        0.9534401885464611
    """
    import argparse
    parser = argparse.ArgumentParser(
        description=(
            "Exact Binomial Confidence Interval Calculator, "
            f"ver. {__version__}"),
        # f"ver. {ebcic.__version__}"),
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
        "-r",
        "--rej-perc-lower",
        type=float,
        help=(
            "set percentage (0.0<= x <50.0) of lower rejection area in "
            "assuming population"
        )
    )
    parser.add_argument(
        "-s",
        "--rej-perc-upper",
        type=float,
        help=(
            "set percentage (0.0<= x <50.0) of upper rejection area in "
            "assuming population"
        )
    )
    parser.add_argument(
        "-c",
        "--confi-perc",
        type=float,
        help=(
            "set confidence percentage for two-sided of 0<k<n where "
            "0 < CONFI_PERC < 100, or for one-sided of k=0 or k=n; "
            "for one-sided of 0<k<n, set "
            "CONFI_PERC = (2 * confi_perc_for_one_sided - 100) "
            "where 50 < confi_perc_for_one_sided < 100"
        )
    )
    parser.add_argument(
        "-a",
        "--alpha",
        type=float,
        help=(
            "ALPHA = 1 - CONFI_PERC/100, "
            "cannot be used with -c (--confi-perc)"
        )
    )
    parser.add_argument(
        "-u",
        "--upper",
        action='store_true',
        help="print upper interval"
    )
    parser.add_argument(
        "-l",
        "--lower",
        action='store_true',
        help="print lower interval"
    )
    args = parser.parse_args(cli_args)
    # args check should be done in each function.
    params = Params(
        k=args.errors,          # Number of errors
        n=args.trials,          # number of trials
        confi_perc=args.confi_perc,
        alpha=args.alpha,
        rej_perc_lower=args.rej_perc_lower,
        rej_perc_upper=args.rej_perc_upper,
    )

    # body
    if args.lower or args.upper:
        # interval = ebcic.exact(params)
        interval = exact(params)
        if args.lower:
            print(interval[0])
        if args.upper:
            print(interval[1])
    else:
        # more info
        # ebcic.print_interval(params)
        print_interval(params)


def round_h_up(x, sig_digits=-5):
    """Round half up.

    Round half up $x$ to the decimal point corresponding to $10^sig_digits$.

    Args:
        x (int or float): Input value.
        sig_digits (int): Significant digits.

    Returns:
        int or float: The rounded-half-up number of x.

    Examples:
        >>> round_h_up(123.45678, -1)
        123.5

        >>> round_h_up(123.45678,  0)
        123.0

        >>> round_h_up(123.45678,  1)
        120.0
    """
    # Input check
    if isinstance(x, float):
        x_type = 'float'
    elif isinstance(x, int):
        x_type = 'int'
    else:
        raise TypeError("'x' must be either float or integer.")
    if not isinstance(sig_digits, int):
        sig_digits = int(sig_digits)
        # raise TypeError("'sig_digits' must be an integer")

    # dec_input for Decimal()
    if sig_digits == 0:
        dec_input = str(0)
    elif sig_digits < 0:
        dec_input = str(10**(sig_digits))
    else:  # sig_digits > 0
        dec_input = "1E" + str(sig_digits)

    # ROUND_HALF_UP
    tmp = Decimal(str(x)).quantize(Decimal(dec_input), rounding=ROUND_HALF_UP)
    if x_type == 'float':
        rounded = float(tmp)
    else:
        rounded = int(tmp)
    return rounded


def check_n(n):
    if (n < 1):
        logger.error(f"'1 <= n' must hold where n={n}!!")
        sys.exit(1)


def check_k(k):
    if (k < 0):
        logger.error(f"'0 <= k' must hold where k={k}!!")
        sys.exit(1)


def check_between_k_and_n(k, n):
    if (n < k):
        logger.error(f"'k <= n' must hold where k={k} and n={n}!!")
        sys.exit(1)


def confi_perc_to_alpha_wo_check(confi_perc):
    return round_h_up(1 - confi_perc / 100)


def alpha_to_confi_perc_wo_check(alpha):
    return round_h_up(100 * (1 - alpha))


def check_confi_perc(confi_perc):
    check_range_float(0., False, 100., False, confi_perc, 'confi_perc')


def check_alpha(alpha):
    check_range_float(0., False, 1., False, alpha, 'alpha')


def check_range_float(
        min: float, min_inclusive: bool,
        max: float, max_inclusive: bool,
        value: float, val_name: str, do_exit: bool = True) -> bool:
    """Check the range of the float variable.
    Args:
        min: Minimum of `value`.
        min_inclusive: '<=' if True '<' otherwise.
        max: Maximum of `value`.
        max_inclusive: '<=' if True '<' otherwise.
        value: Float value to check.
        val_name: The value's variable name.

    Returns:
        True: If holds.
        False: If not hold and do_exit=False.

    Examples:
        >>> check_range_float(0., True, 1., False, 0, 'test1', do_exit=False)
        True

        >>> check_range_float(0., False, 1., False, 0, 'test2', do_exit=False)
        False

        >>> check_range_float(0., True, 1., True, 1, 'test3', do_exit=False)
        True

        >>> check_range_float(0., True, 1., False, 1, 'test4', do_exit=False)
        False

    """
    if min_inclusive:
        min_inequality = ' <= '
        min_holds = (min <= value)
    else:
        min_inequality = ' < '
        min_holds = (min < value)
    if max_inclusive:
        max_inequality = ' <= '
        max_holds = (value <= max)
    else:
        max_inequality = ' < '
        max_holds = (value < max)
    if not (min_holds and max_holds):
        if do_exit:
            logger.error(
                "'" + str(min) + min_inequality + val_name + max_inequality
                + str(max) + "' must hold where '" + val_name +
                f" = {value}'!")
            sys.exit(1)
        else:
            return False
    return True


class Params:
    """Parameter pool class for Binomial Confidence Interval.

    This stores parameters regarding binomial confidence intervals after
    checking and avoiding duplication.

    Set the confidence with either of the following parameters:
        - `confi_perc`
        - `alpha`
        - both or either of `rej_perc_upper` and `rej_perc_lower`
        - both or either of `rej_upper` and `rej_lower`

    Args:
        k (int):
            Number of errors.
        n (int):
            Number of trials.
        rej_perc_upper (float):
            Upper rejection area of assuming population in percentage.
        rej_perc_lower (float):
            Lower rejection area of assuming population in percentage.
        rej_upper (float):
            Upper rejection area of assuming population in ratio, or
            rej_perc_upper/100.
        rej_lower (float):
            Lower rejection area of assuming population in ratio, or
            rej_perc_lower/100.
        alpha (float):
            Rejection area of assuming population in ratio
            where 0 < alpha < 1.
            If `alpha` is not given, either of the following is set
            depending on the availability:
            - `rej_upper + rej_lower`
            - `rej_perc_upper/100 + rej_perc_lower/100`
            - `1 - confi_perc/100`
        confi_perc (float):
            Confidence percentage (confidence coefficient in percentage
            or confidence level) for two-sided of 0<k<n
            where 0 < confi_perc < 100, or for one-sided of k=0 or k=n.
            For one-sided of 0<k<n, use rej_* variables or set
            `confi_perc=(2 * s - 100)` where `s` is
            the confidence percentage for one-sided of 0<k<n and
            `50 < s < 100.
            If `confi_perc` is not given, either of the following is set
            depending on the availability:
            - `(1 - alpha)*100`
            - `(1 - (rej_upper + rej_lower))*100`
            - `100 - (rej_perc_upper + rej_perc_lower)`
        warn_line_confi_perc (float):
            Warned, if `confi_perc` is smaller than this value.
    """

    def __init__(
            self,
            k=None,
            n=None,
            # confidence parameters
            alpha=None,
            confi_perc=None,
            rej_lower=None,
            rej_upper=None,
            rej_perc_lower=None,
            rej_perc_upper=None,
            warn_line_confi_perc=60,  # for two-sided confidence
            # FIDO Biometrics Requirements v2.0
            # https://fidoalliance.org/specs/biometric/requirements/
            # uses 80% for one-sided
            # (which corresponds with 60% for two-sided) confidence.
            # (confi_perc < warn_line_confi_perc) is warned
            ):
        self.n = None
        self.k = None
        # self.itvls = None  # dict of interval values
        # NOTE:
        #   self.rej*, self.alpha, self.confi_per can be moved to self.itvls.
        self.rej_lower = None
        self.rej_upper = None
        self.alpha = None
        self.confi_perc = None
        # 'warn_line' is set to both confi_perc and alpha
        check_confi_perc(warn_line_confi_perc)  # may exit
        self.warn_line_confi_perc = warn_line_confi_perc
        self.warn_line_alpha = confi_perc_to_alpha_wo_check(
            warn_line_confi_perc)
        self.num_of_warned = 0  # for test
        # Check and set parameters.
        # Set k and n
        if n is not None:
            self.set_n(n)
        if k is not None:
            self.set_k(k)
        # 'rej_perc', 'rej_{lower,upper}', 'confi_perc' or 'alpha'
        itvl_input_method = None
        if (rej_perc_lower is not None) or (rej_perc_upper is not None):
            itvl_input_method = 'rej_perc'
        if (rej_lower is not None) or (rej_upper is not None):
            if (itvl_input_method is not None):
                itvl_input_method = 'duplicated'
            else:
                itvl_input_method = 'rej'
        if (confi_perc is not None) or (alpha is not None):
            if (itvl_input_method is not None):
                itvl_input_method = 'duplicated'
            else:
                if (confi_perc is not None) and (alpha is not None):
                    logger.warning("giving either 'alpha' or "
                                   "'confi_perc' is enough.")
                    if alpha != confi_perc_to_alpha_wo_check(
                            warn_line_confi_perc):
                        logger.error(
                            f"'alpha [{alpha}] != "
                            f"confi_perc_to_alpha_wo_check(confi_perc "
                            f"[confi_perc])"
                            f"[{confi_perc_to_alpha_wo_check(confi_perc)}]'!!")
                        sys.exit(1)
                    # process as 'confi_perc' even for 'alpha'
                    itvl_input_method = 'confi_perc'
                else:  # either (confi_perc is not None) or (alpha is not None)
                    if alpha is not None:
                        itvl_input_method = 'alpha'
                    else:
                        itvl_input_method = 'confi_perc'
        # print("itvl_input_method =", itvl_input_method)
        # Check and set 'confi_perc', 'alpha' and 'rej_*'.
        area_err = False
        if itvl_input_method is None:
            if warn_line_confi_perc is None:
                area_err = True
            elif self.k is not None and self.n is not None:
                # not warn_line set only
                area_err = True
            # else: warn_line set only
        elif itvl_input_method == 'duplicated':
            area_err = True
        else:
            if itvl_input_method == 'confi_perc':
                self.set_confi_perc(confi_perc)
            elif itvl_input_method == 'alpha':
                self.set_alpha(alpha)
            elif itvl_input_method == 'rej_perc':
                self.set_rej_perc(rej_perc_lower, rej_perc_upper)
            elif itvl_input_method == 'rej':
                self.set_rej(rej_lower, rej_upper)
            else:
                logger.error(
                    "unknown itvl_input_method="
                    f"{itvl_input_method}!!")
                sys.exit(1)
        if area_err:
            logger.error(
                "give confidence with 'confi_perc' (-c), 'alpha' (-a), or "
                "'rej_perc_lower' (-r) and/or 'rej_perc_upper' (-s)!")
            sys.exit(1)

    def set_k(self, k):
        check_k(k)  # may exit
        if (self.n is not None):
            check_between_k_and_n(k, self.n)  # may exit
        self.k = k

    def set_n(self, n):
        check_n(n)  # may exit
        if (self.k is not None):
            check_between_k_and_n(self.k, n)  # may exit
        self.n = n

    def check_and_warn_confi_perc(self, confi_perc):
        check_confi_perc(confi_perc)
        if confi_perc < self.warn_line_confi_perc:
            logger.warning(
                f"confi_perc={confi_perc} is too small "
                f"(than {self.warn_line_confi_perc}).")
            warned = 1
            self.num_of_warned += 1
        else:
            warned = 0
        return warned

    def interval_update(self, itvls: dict) -> bool:
        """Update confidence parameters if new to the stored ones.

        Args:
            itvls: Dict of interval values.

        Returns:
            warned: The number of warned.
        """
        warned = 0
        if (
                (self.confi_perc != itvls['confi_perc'])
                or (self.alpha != itvls['alpha'])
                or (self.rej_lower != itvls['rej_lower'])
                or (self.rej_upper != itvls['rej_upper'])
                ):
            # check
            tmp = itvls['confi_perc']
            warned = self.check_and_warn_confi_perc(tmp)
            # set
            # self.alpha_or_confi_perc_checked = True
            self.confi_perc = tmp
            # self.alpha_checked = True
            self.alpha = itvls['alpha']
            # rejection area
            self.rej_lower = itvls['rej_lower']
            self.rej_upper = itvls['rej_upper']
            # warning to rej_*
            if (self.n is not None) and (self.k is not None):
                nh = self.n / 2.
                if self.rej_upper == 0. and self.k > nh:
                    logger.warning(
                        "rej_upper=0 is not appropriate for k>n/2.")
                elif self.rej_lower == 0. and self.k < nh:
                    logger.warning(
                        "rej_lower=0 is not appropriate for k<n/2.")
        return warned

    def set_confi_perc(self, confi_perc: float):
        """Set confidence percentage and more.

        This also sets its corresponding 'alpha'.

        Args:
            confi_perc (0 < x < 100):
                Confidence percentage.

        Returns:
            warned: The number of warned.

        Note:
            v0.0.4 and newer do not set Z_{alpha} and Z_{alpha/2}
            but generate them from 'alpha'.
        """
        alpha = confi_perc_to_alpha_wo_check(confi_perc)
        # rejection area
        rej_lower = None
        rej_upper = None
        itvls = {
            'confi_perc': confi_perc,
            'alpha': alpha,
            'rej_lower': rej_lower,
            'rej_upper': rej_upper,
        }
        return self.interval_update(itvls)

    def check_and_warn_alpha(self, alpha):
        check_alpha(alpha)  # may exit
        if alpha > self.warn_line_alpha:
            logger.warning(
                f"'alpha > {self.warn_line_alpha}'"
                " is unusual.\n"
                f"Is 'alpha={alpha}' correct?")
            warned = 1
            self.num_of_warned += 1
        else:
            warned = 0
        return warned

    def set_alpha(self, alpha):
        """Set alpha, i.e. (1 - (confidence percentage)/100), and more.

        This also sets its corresponding
        'perc_confi' (confidence percentage).

        Args:
            alpha (0 < x < 1): (1 - (confidence percentage)/100).

        Returns:
            warned: The number of warned.

        Note:
            v0.0.4 and newer do not set Z_{alpha} and Z_{alpha/2}
            but generate them from 'alpha'.
        """
        confi_perc = alpha_to_confi_perc_wo_check(alpha)
        # rejection area
        rej_lower = None
        rej_upper = None
        itvls = {
            'confi_perc': confi_perc,
            'alpha': alpha,
            'rej_lower': rej_lower,
            'rej_upper': rej_upper,
        }
        return self.interval_update(itvls)

    def set_rej_perc(
            self,
            rej_perc_lower: float = 0., rej_perc_upper: float = 0.) -> int:
        """Set rejection area in percentage.

        This also sets its corresponding
        'confi_perc' (confidence percentage) and 'alpha'.

        Args:
            rej_perc_lower (0 <= x < 50):
                Percentage of lower rejection area in assuming population.
            rej_perc_upper (0 <= x < 50):
                Percentage of upper rejection area in assuming population.

        Returns:
            warned: The number of warned.
        """
        if (rej_perc_lower is None):
            rej_perc_lower = 0.
            rej_lower = 0.
        else:
            rej_lower = rej_perc_lower / 100.
        if (rej_perc_upper is None):
            rej_perc_upper = 0.
            rej_upper = 0.
        else:
            rej_upper = rej_perc_upper / 100.
        # Check. They may exit.
        check_range_float(0., True, 0.5, False, rej_lower, 'rej_lower')
        check_range_float(0., True, 0.5, False, rej_upper, 'rej_upper')
        # alpha is checked with confi_perc in interval_update()
        alpha = rej_lower + rej_upper
        confi_perc = 100. - (rej_perc_lower + rej_perc_upper)
        # self.check_and_warn_confi_perc(confi_perc)
        itvls = {
            'confi_perc': confi_perc,
            'alpha': alpha,
            'rej_lower': rej_lower,
            'rej_upper': rej_upper,
        }
        # print("itvls =", itvls)
        return self.interval_update(itvls)

    def set_rej(self, rej_lower: float = 0., rej_upper: float = 0.) -> int:
        """Set rejection area in ratio.

        This also sets its corresponding
        'perc_confi' (confidence percentage) and 'alpha'.

        Args:
            rej_lower (0 <= x < 0.5):
                Lower rejection area of assuming population in ratio.
            rej_upper (0 <= x < 0.5):
                Upper rejection area of assuming population in ratio.

        Returns:
            warned: The number of warned.
        """
        if (rej_lower is None):
            rej_lower = 0.
        if (rej_upper is None):
            rej_upper = 0.
        # Check. They may exit.
        check_range_float(0., True, 0.5, False, rej_lower, 'rej_lower')
        check_range_float(0., True, 0.5, False, rej_upper, 'rej_upper')
        # alpha is checked with confi_perc in interval_update()
        alpha = rej_lower + rej_upper
        confi_perc = 100. * (1 - (rej_lower + rej_upper))
        itvls = {
            'confi_perc': confi_perc,
            'alpha': alpha,
            'rej_lower': rej_lower,
            'rej_upper': rej_upper,
        }
        return self.interval_update(itvls)

    def confi_perc_to_alpha(self, confi_perc):
        """Confidence percentage to alpha.

        This converts 'confi_perc' given in percentage between
        0 to 100 exclusive
        to alpha in between 0 to 1 exclusive.
        """
        self.set_confi_perc(confi_perc)
        alpha = confi_perc_to_alpha_wo_check(confi_perc)
        return alpha

    def alpha_to_confi_perc(self, alpha):
        """Alpha to confidence percentage.

        Convert given 'alpha' in between 0 to 1 exclusive
        to 'confi_perc' in percentage between 0 to 100 exclusive.
        """
        self.set_alpha(alpha)
        confi_perc = alpha_to_confi_perc_wo_check(alpha)
        return confi_perc

    def exact(self):
        """Return exact Binomial confidence interval for the parameter.

        Args:
            params (Params): Instance including k, n and
                alpha (or confi_perc).

        Returns:
            tuple: tuple containing:
                lower_p (float): lower interval of Binomial confidence.
                upper_p (float): upper interval of Binomial confidence.

        Examples:
            >>> Params(k=1, n=10000, confi_perc=95.0).exact()
            (2.5317775934704154e-06, 0.0005570369979470232)

            >>> Params(k=1, n=10000, confi_perc=90.0).exact()
            (5.129316283730961e-06, 0.00047429765916542926)

            >>> Params(k=1, n=10000, rej_perc_lower=2.5, rej_perc_upper=2.5).exact()  # noqa: E501
            (2.5317775934704154e-06, 0.0005570369979470232)

            >>> Params(k=1, n=10000, rej_perc_lower=5.).exact()
            (0.0, 0.00047429765916542926)

            >>> Params(k=9999, n=10000, confi_perc=90.0).exact()
            (0.9995257023408346, 0.9999948706837163)

            >>> Params(k=9999, n=10000, rej_perc_upper=5.).exact()
            (0.9995257023408346, 1.0)
            """
        n = self.n
        k = self.k

        if k == 0 or k == n:
            # one-sided
            lower_p, upper_p = self.exact_border()
        else:
            # 0 < k < n
            nh = n / 2
            if (self.rej_lower is None) and (self.rej_upper is None):
                # even two-sided
                rl = self.alpha / 2.
                ru = rl
            else:
                # odd two-sided
                rl = self.rej_lower
                ru = self.rej_upper
            # print("ru =", ru, " rl =", rl)
            reverse_mode = False
            # Cumulative error becomes large for k >> n/2,
            # so make k < n/2.
            if k > nh:
                reverse_mode = True
                k = n - k
                tmp = rl
                rl = ru
                ru = tmp

            def upper(p):
                """
                Upper interval is the p making the following 0.
                """
                return binom.cdf(k, n, p) - rl

            def lower(p):
                """
                Lower interval for k>0 is the p making the following 0.
                """
                return binom.cdf(n - k, n, 1 - p) - ru
            # Alternatives:
            # return 1-sum([binom.pmf(i,n,result_lower) for i in range(k)]) - r
            # return sum([binom.pmf(i,n,result_lower) for i in range(k,n)]) - r

            # set init value
            # 0 < k < n
            tmp = k / n
            u_init = tmp
            l_init = tmp
            if k == 1:
                l_init = 0
            elif k == 2:
                l_init = k / (2 * n)
            # print(f"l_init={l_init}, u_init={u_init}")
            # solve lower
            if ru == 0.:
                lower_p = 0.
            else:
                lower_p = optimize.fsolve(lower, l_init)[0]
            # solve upper
            if rl == 0.:
                upper_p = 1.
            elif (
                    (n % 2) == 0 and k == nh and
                    (self.rej_lower is None) and (self.rej_upper is None)):
                # even two-sided, k/n = 1/2, and n is even
                upper_p = 1 - lower_p
            else:
                upper_p = optimize.fsolve(upper, u_init)[0]

            if reverse_mode:
                tmp = lower_p
                lower_p = 1 - upper_p
                upper_p = 1 - tmp
        return lower_p, upper_p

    def exact_border(self) -> tuple:
        """Exact Binomial Confidence Interval for k==0 or k==n.

        Return the exact binomial confidence interval for given parameters
        for k==0 or k==n.

        Args:
            params (Params):
                Instance including k, n and
                alpha (or confi_perc) where k==0 or k==n.

        Returns:
            tuple: tuple containing:
                lower_p (float):
                    lower interval of Binomial confidence
                    for k==0 or k==n.
                upper_p (float):
                    upper interval of Binomial confidence
                    for k==0 or k==n.

        Examples:
            >>> Params(k=0, n=10000, confi_perc = 95.0).exact_border()[1]
            0.00029952835977664627

            >>> round_h_up(Params(k=0, n=10000, confi_perc = 95.0).exact_border()[1],-7)  # noqa: E501
            0.0002995

            >>> Params(k=0, n=10000, rej_perc_lower = 5.0).exact_border()[1]
            0.00029952835977664627

            >>> Params(k=10000, n=10000, confi_perc = 95.0).exact_border()[0]
            0.9997004716402234

            >>> Params(k=10000, n=10000, rej_perc_upper = 5.).exact_border()[0]
            0.9997004716402234
        """
        # NOTE: Do not replace self.alpha with self.alpha/2.
        #       This function is only for one-sided confidence intervals.
        alpha = self.alpha
        tmp = (alpha)**(1 / self.n)
        if self.k == 0:
            if (self.rej_upper is not None) and (self.rej_upper != 0.):
                logger.error(
                    "'rej_upper' and 'rej_perc_upper' shall be 0.0 for 'k=0'!")
                sys.exit(1)
            lower_p = 0.0
            upper_p = 1 - tmp
        elif self.k == self.n:
            if (self.rej_lower is not None) and (self.rej_lower != 0.):
                logger.error(
                    "'rej_lower' and 'rej_perc_lower' shall be 0.0 for 'k=n'!")
                sys.exit(1)
            lower_p = tmp
            upper_p = 1.0
        else:  # 0 < k < n
            # Use exact()
            logger.error(
                "either k==0 or k==n must hold where "
                f"k={self.k} and n={self.n}!!")
            sys.exit(1)
        return lower_p, upper_p

    # Approximated intervals
    def rule_of_ln_alpha(self):
        """ Generalized rule of three.

        Interval of rule of -ln(alpha), i.e generalized version of
        'rule of three'.
        This can be used only when (k == 0 or k == n) and
        reliable if n is large enough, say n > 50.

        Args:
            params (Params): Instance including k, n, alpha (or confi_perc).

        Returns:
            tuple: tuple containing:
                lower_p (float): lower interval of Binomial confidence.
                upper_p (float): upper interval of Binomial confidence.

        Examples:
            >>> Params(k=0, n=10000, confi_perc = 95.0).rule_of_ln_alpha()[1]
            0.0002995732273553991

            >>> round_h_up(Params(k=0, n=10000, confi_perc = 95.0).rule_of_ln_alpha()[1],-7)  # noqa: E501
            0.0002996

            >>> Params(k=0, n=10000, rej_perc_lower = 5.0).rule_of_ln_alpha()[1]
            0.0002995732273553991

            >>> Params(k=10000, n=10000, rej_perc_upper = 5.0).rule_of_ln_alpha()[0]
            0.9997004267726446
        """
        n = self.n
        k = self.k
        alpha = self.alpha
        # print(f"alpha={alpha}")
        # Get interval
        lower_p = None
        upper_p = None
        if k == 0:
            if (self.rej_upper is not None) and (self.rej_upper != 0.):
                logger.error(
                    "'rej_upper' and 'rej_perc_upper' shall be 0.0 for 'k=0'!")
                sys.exit(1)
            lower_p = 0
            upper_p = -np.log(alpha) / n
        elif k == n:
            if (self.rej_lower is not None) and (self.rej_lower != 0.):
                logger.error(
                    "'rej_lower' and 'rej_perc_lower' shall be 0.0 for 'k=n'!")
                sys.exit(1)
            lower_p = 1 - (-np.log(alpha) / n)
            upper_p = 1
        else:  # 0 < k < n
            # raise ValueError("either k==0 or k==n must hold!!")
            logger.error(
                f"either 'k=0' or 'k=n' must hold where 'k={k}' and 'n={n}'!")
            sys.exit(1)
        return lower_p, upper_p

    def beta_approx(self):
        """Approximated interval using beta distribution.

        Good approximation.

        Args:
            params (Params): Instance including k, n, alpha (or confi_perc).

        Returns:
            tuple: tuple containing:
                lower_p (float): lower interval of Binomial confidence.
                upper_p (float): upper interval of Binomial confidence.

        Examples:
            >>> Params(k=0, n=10000, confi_perc=95.0).beta_approx()[1]
            0.0002995283597770252

            >>> Params(k=10000, n=10000, confi_perc=95.0).beta_approx()[0]
            0.999700471640223

            >>> Params(k=1, n=10000, confi_perc=95.0).beta_approx()
            (2.531777593461957e-06, 0.000557036997946958)

            >>> Params(k=1, n=10000, rej_perc_lower=2.5, rej_perc_upper=2.5).beta_approx()  # noqa E501
            (2.531777593461957e-06, 0.000557036997946958)

            >>> Params(k=1, n=10000, confi_perc=90.0).beta_approx()[1]
            0.00047429765916540134

            >>> Params(k=1, n=10000, rej_perc_lower=5.).beta_approx()[1]
            0.00047429765916540134

            >>> Params(k=9999, n=10000, confi_perc=90.0).beta_approx()[0]
            0.9995257023408346

            >>> Params(k=9999, n=10000, rej_perc_upper=5.).beta_approx()[0]
            0.9995257023408346
        """
        n = self.n
        k = self.k
        alpha = self.alpha

        if k == 0:
            # one-sided
            lower_p = 0.0
            upper_p = beta.ppf(1-alpha, k+1, n-k)
        elif k == n:
            # one-sided
            lower_p = beta.ppf(alpha, k, n-k+1)
            upper_p = 1.0
        else:
            # 0 < k < n
            if (self.rej_lower is None) and (self.rej_upper is None):
                # even two-sided
                rl = alpha / 2.
                ru = rl
            else:
                # odd two-sided or one-sided
                rl = self.rej_lower
                ru = self.rej_upper
            lower_p = beta.ppf(ru, k, n-k+1)
            upper_p = beta.ppf(1-rl, k+1, n-k)
        return lower_p, upper_p

    def normal_approx(self):
        """Approximated interval using normal distribution

        Interval obtained by approximating Binomial distribution
        to normal distribution, which is available only 0<k<n.

        Args:
            params (Params): Instance including k, n and 'alpha'
                (or 'confi_perc').

        Returns:
            tuple: tuple containing:
                lower_p (float): lower interval of Binomial confidence.
                upper_p (float): upper interval of Binomial confidence.

        Note:
            - Normal approximation does not give good approximation for small
              n or k being close to 0 or n.
            - v0.0.4 and newer do not set Z_{alpha} and Z_{alpha/2} to Params
              but generate them from 'alpha'.

        Examples:
            >>> Params(k=1, n=10000, confi_perc=95.0).normal_approx()
            (0.0, 0.00029600580053504625)

            >>> Params(k=1, n=10000, rej_perc_lower=2.5, rej_perc_upper=2.5).normal_approx()  # noqa E501
            (0.0, 0.00029600580053504625)

            >>> Params(k=1, n=10000, confi_perc = 90.0).normal_approx()[1]
            0.00026449322486687016

            >>> Params(k=1, n=10000, rej_perc_lower = 5.0).normal_approx()[1]
            0.00026449322486687016

            >>> Params(k=9999, n=10000, confi_perc = 90.0).normal_approx()[0]
            0.9997355067751331

            >>> Params(k=9999, n=10000, rej_perc_upper = 5.0).normal_approx()[0]
            0.9997355067751331
        """
        n = self.n
        k = self.k
        p = k / n

        if n != 1:
            sigma = np.sqrt(n * p * (1 - p))/(n - 1)
        else:
            # n == 1
            # sigma = np.sqrt(p * (1 - p)/n)
            sigma = np.sqrt(p * (1 - p))

        if (self.rej_lower is None) and (self.rej_upper is None):
            # use alpha (= rej_lower + rej_upper)
            half_width = alpha_to_za(self.alpha / 2.) * sigma
            lower_p = max(0., p - half_width)
            upper_p = min(1., p + half_width)
        else:
            # use rej_lower and rej_upper
            if self.rej_lower == 0.:
                upper_p = 1.
            else:
                half_width_u = alpha_to_za(self.rej_lower) * sigma
                upper_p = min(1., p + half_width_u)
            if self.rej_upper == 0.:
                lower_p = 0.
            else:
                half_width_l = alpha_to_za(self.rej_upper) * sigma
                lower_p = max(0., p - half_width_l)
        return lower_p, upper_p

    def wilson_score(self):
        """Wilson score interval.

        Args:
            params (Params): Instance including k, n and 'alpha'
                (or 'confi_perc').

        Returns:
            tuple: tuple containing:
                lower_p (float): lower interval of Binomial confidence.
                upper_p (float): upper interval of Binomial confidence.

        Note:
            - Upper interval seems close to the exact one, but
              lower one seems a little greater than the exact one.
            - v0.0.4 and newer do not set Z_{alpha} and Z_{alpha/2} to Params
              but generate them from 'alpha'.

        Examples:
            >>> Params(k=0, n=10000, confi_perc = 95.0).wilson_score()[1]
            0.00027047997304067333

            >>> Params(k=10000, n=10000, confi_perc = 95.0).wilson_score()[0]
            0.9997295200269593

            >>> Params(k=1, n=10000, confi_perc = 95.0).wilson_score()
            (1.76527238382147e-05, 0.0005662672867662839)

            >>> Params(k=1, n=10000, rej_perc_lower=2.5, rej_perc_upper=2.5).wilson_score()  # noqa E501
            (1.76527238382147e-05, 0.0005662672867662839)

            >>> Params(k=9999, n=10000, confi_perc = 95.0).wilson_score()
            (0.9994337327132337, 0.9999823472761619)

            >>> Params(k=9999, n=10000, rej_perc_lower=2.5, rej_perc_upper=2.5).wilson_score()  # noqa E501
            (0.9994337327132337, 0.9999823472761619)

            >>> Params(k=1, n=10000, confi_perc = 90.).wilson_score()
            (2.230960071836048e-05, 0.0004481162763277046)

            >>> Params(k=1, n=10000, rej_perc_lower=5.).wilson_score()
            (0.0, 0.0004481162763277046)

            >>> Params(k=9999, n=10000, confi_perc = 90.0).wilson_score()
            (0.9995518837236722, 0.9999776903992816)

            >>> Params(k=9999, n=10000, rej_perc_upper=5.).wilson_score()
            (0.9995518837236722, 1.0)
        """
        n = self.n
        k = self.k
        p = k / n

        def wilson_components(
                z: float, k: int, n: int, p: float) -> (float, float, float):
            mu = 2 * k + z**2
            half_width = z * np.sqrt(z**2 + 4 * k * (1 - p))
            denomi = 2 * (n + z**2)
            return mu, half_width, denomi

        if k == 0 or k == n:
            # one-sided
            z = alpha_to_za(self.alpha)
            mu, half_width, denomi = wilson_components(z, k, n, p)
            if k == 0:
                lower_p = 0.0
                upper_p = min(1, (mu + half_width) / denomi)
            else:  # k == n:
                lower_p = max(0, (mu - half_width) / denomi)
                upper_p = 1.0
        else:
            # two-sided
            if (self.rej_lower is None) and (self.rej_upper is None):
                # use alpha (= rej_lower + rej_upper)
                z = alpha_to_za(self.alpha / 2.)
                mu, half_width, denomi = wilson_components(z, k, n, p)
                lower_p = max(0., (mu - half_width) / denomi)
                upper_p = min(1., (mu + half_width) / denomi)
            else:
                # use rej_lower and rej_upper
                if self.rej_lower == 0.:
                    upper_p = 1.
                else:
                    z = alpha_to_za(self.rej_lower)
                    mu, half_width, denomi = wilson_components(z, k, n, p)
                    upper_p = min(1., (mu + half_width) / denomi)
                if self.rej_upper == 0.:
                    lower_p = 0.
                else:
                    z = alpha_to_za(self.rej_upper)
                    mu, half_width, denomi = wilson_components(z, k, n, p)
                    lower_p = max(0., (mu - half_width) / denomi)
        return lower_p, upper_p

    def wilson_score_cc(self):
        """Wilson score interval with continuity correction.

        Args:
            params (Params): Instance including k, n and 'alpha'
            (or 'confi_perc').

        Returns:
            tuple: tuple containing:
                lower_p (float): lower interval of Binomial confidence.
                upper_p (float): upper interval of Binomial confidence.

        Note:
            - Lower interval is not a good approximation for small k,
              i.e. it becomes higher than the exact lower interval,
              though it is better than Wilson score interval
              (without continuity correction).
            - v0.0.4 and newer do not set Z_{alpha} and Z_{alpha/2} to Params
              but generate them from 'alpha'.

        Examples:
            >>> Params(k=0, n=10000, confi_perc = 95.0).wilson_score_cc()
            (0.0, 0.00025428323087746505)

            >>> Params(k=10000, n=10000, confi_perc = 95.0).wilson_score_cc()
            (0.9996364213055936, 1.0)

            >>> Params(k=1, n=10000, confi_perc = 95.0).wilson_score_cc()
            (5.220053865461516e-06, 0.000578699956739037)

            >>> Params(k=1, n=10000, rej_perc_lower=2.5, rej_perc_upper=2.5).wilson_score_cc()  # noqa E501
            (5.220053865461516e-06, 0.000578699956739037)

            >>> Params(k=9999, n=10000, confi_perc = 95.0).wilson_score_cc()
            (0.9993507610423017, 1.0)

            >>> Params(k=9999, n=10000, rej_perc_lower=2.5, rej_perc_upper=2.5).wilson_score_cc()  # noqa E501
            (0.9993507610423017, 1.0)

            >>> Params(k=1, n=10000, confi_perc = 90.0).wilson_score_cc()
            (6.874230637052348e-06, 0.0004635516464090127)

            >>> Params(k=1, n=10000, rej_perc_lower = 5.).wilson_score_cc()
            (0.0, 0.0004635516464090127)

            >>> Params(k=9999, n=10000, confi_perc = 90.0).wilson_score_cc()
            (0.9994722211307717, 1.0)

            >>> Params(k=9999, n=10000, rej_perc_upper = 5.).wilson_score_cc()
            (0.9994722211307717, 1.0)
        """
        n = self.n
        k = self.k
        p = k / n

        def wilson_cc_components(
                z: float, k: int, n: int, p: float) -> (float, float, float):
            mu = 2 * k + z**2
            half_width = 1 + z * np.sqrt(
                z**2 - 1 / n + 4 * k * (1 - p) + (4 * p - 2))
            denomi = 2 * (n + z**2)
            return mu, half_width, denomi

        if k == 0 or k == n:
            # one-sided
            z = alpha_to_za(self.alpha)
            mu, half_width, denomi = wilson_cc_components(z, k, n, p)
            if k == 0:
                lower_p = 0.0
                upper_p = min(1, (mu + half_width) / denomi)
            else:  # k == n:
                lower_p = max(0, (mu - half_width) / denomi)
                upper_p = 1.0
        else:
            # two-sided
            if (self.rej_lower is None) and (self.rej_upper is None):
                # use alpha (= rej_lower + rej_upper)
                z = alpha_to_za(self.alpha / 2.)
                mu, half_width, denomi = wilson_cc_components(z, k, n, p)
                lower_p = max(0., (mu - half_width) / denomi)
                upper_p = min(1., (mu + half_width) / denomi)
            else:
                # use rej_lower and rej_upper
                if self.rej_lower == 0.:
                    upper_p = 1.
                else:
                    z = alpha_to_za(self.rej_lower)
                    mu, half_width, denomi = wilson_cc_components(z, k, n, p)
                    upper_p = min(1., (mu + half_width) / denomi)
                if self.rej_upper == 0.:
                    lower_p = 0.
                else:
                    z = alpha_to_za(self.rej_upper)
                    mu, half_width, denomi = wilson_cc_components(z, k, n, p)
                    lower_p = max(0., (mu - half_width) / denomi)
        return lower_p, upper_p

    def verify_interval_of_p(
            self, lower_p, upper_p, sig_digits=-5, verbose=1):
        """Check the reliability of the given interval.

        This checks, for given parameteres, if integral of outside of
        the interval is alpha or not, and so on.

        Args:
            params (Params): Instance including k, n, alpha (or confi_perc).
            lower_p: Lower interval of p for the params.
            upper_p: Upper interval of p for the params.
            sig_digits: Significant digits.
            verbose:
                0: no message
                1: show only errors
                2: show always

        Returns:
            bool: 1 if unreliable, 0 reliable,
            i.e. no evidence was found to consider them unreliable.
        """
        MIN_SIG_DIGITS = -10
        n = self.n
        k = self.k
        alpha = self.alpha

        ret = 0
        if sig_digits < MIN_SIG_DIGITS:
            sig_digits = MIN_SIG_DIGITS
        acceptable_range = 10**sig_digits
        # Test that cumulative error = r
        if (k == 0) or (k == n):
            if (k == 0):
                if abs(binom.cdf(k, n, upper_p) - alpha) > acceptable_range:
                    ret = 1
                if (verbose == 1 and ret == 1) or verbose >= 2:
                    print(
                        f"Upper : check if abs(alpha ({alpha}) - "
                        f"cumulative error ({binom.cdf(k, n, upper_p)})) <= "
                        f"{acceptable_range}")
            else:
                # k == n
                if abs(
                        binom.cdf(n - k, n, 1 - lower_p) -
                        alpha) > acceptable_range:
                    ret = 1
                if (verbose == 1 and ret == 1) or verbose >= 2:
                    print(
                        f"Lower : check if abs(alpha ({alpha}) - "
                        "cumulative error "
                        f"({binom.cdf(n - k, n, 1- lower_p)})) "
                        f"<= {acceptable_range}")
        else:
            # 0 < k < n
            # r = alpha / 2
            if (self.rej_lower is None) and (self.rej_upper is None):
                # even two-sided
                rl = self.alpha / 2.
                ru = rl
            else:
                # odd two-sided
                rl = self.rej_lower
                ru = self.rej_upper
            tmp_ret = 0
            if abs(binom.cdf(k, n, upper_p) - rl) > acceptable_range:
                ret = 1
                tmp_ret = 1
            if (verbose == 1 and tmp_ret == 1) or verbose >= 2:
                # alpha/2 or self.rej_lower
                print(
                    f"Upper : check if abs({rl} - "
                    f"cumulative error ({binom.cdf(k, n, upper_p)})) <= "
                    f"{acceptable_range}")
                # print(f"Upper : check if {round_h_up(r, sig_digits)} == "
                #      f"{round_h_up(binom.cdf(k, n, upper_p), sig_digits)}")
            tmp_ret = 0
            if abs(binom.cdf(n - k, n, 1 - lower_p) - ru) > acceptable_range:
                ret = 1
                tmp_ret = 1
            if (verbose == 1 and tmp_ret == 1) or verbose >= 2:
                # alpha/2 or self.rej_upper
                print(
                    f"Lower : check if abs({ru} - "
                    f"cumulative error ({binom.cdf(n - k, n, 1 - lower_p)}))"
                    f" <= {acceptable_range}")

        # Test of position of the intervals
        pp = k / n
        tmp_ret = 0
        if upper_p < pp:
            ret = 1
            tmp_ret = 1
        if (verbose == 1 and tmp_ret == 1) or verbose >= 2:
            print(f"Upper : check if upper interval ({upper_p}) >= k/n ({pp})")
        tmp_ret = 0
        if lower_p > pp:
            ret = 1
            tmp_ret = 1
        if (verbose == 1 and tmp_ret == 1) or verbose >= 2:
            print(f"Lower : check if lower interval ({lower_p}) <= k/n ({pp})")

        return ret

    # ===== Text Output =====
    def print_interval(self):
        """Print exact interval as text.

        Args:
            params (Params): Instance including k, n and
                'confi_perc' (confidence percentage) or 'alpha'.

        Examples:
            >>> print_interval(Params(k=1, n=10000, confi_perc = 95.0))  # noqa: E501  # doctest: +SKIP
        """
        # Num of significant digits.
        # Negative means below decimal point.
        # sig_digits = -10
        sig_digits = -14  # To show lower interval around 10^-8.

        print("\nExact Binomial Confidence Interval Calculator, "
              f"ver. {__version__}\n")

        # Print Parameters
        print("===== Parameters =====")
        print("n (num of trials)  :", self.n)
        print("k (num of failures):", self.k)
        print("k/n (observed p')  :",
              round_h_up(self.k/self.n, sig_digits))

        # Instantiation of parameters where alpha is set as well.
        lower_p, upper_p = exact(self)

        print("\n===== Exact interval of p with ", end="")
        if self.k == 0 or self.k == self.n:
            print(f"{self.confi_perc} [%] one-sided confidence", end="")
        else:
            if (self.rej_lower is not None) and (self.rej_upper is not None):
                # odd two-sided or one-sided
                print(
                    "rejection area of lower " +
                    f"{100 * self.rej_lower} [%] and upper " +
                    f"{100 * self.rej_upper} [%]", end="")
            else:
                # even two-sided
                print(
                    f"{self.confi_perc} [%] two-sided (or " +
                    f"{100 - (100 - self.confi_perc)/2}" +
                    " [%] one-sided) confidence", end="")
        print(" =====")

        print("Upper : ", round_h_up(upper_p, sig_digits))
        print("Lower : ", round_h_up(lower_p, sig_digits))
        print("Width : ", round_h_up(upper_p - lower_p, sig_digits))
        '''
        # NOTE: if the numbers are too small, use below.
        print(f"Upper : {result_upper[0]:.6g}")
        print(f"Lower : {result_lower[0]:.6g}")
        print(f"Width : {result_upper[0] - result_lower[0]:.6g}")
        # NOTE: be careful using round(),
            which rounds 0.*5 to the nearest even,
            say, 1.2345 to 1.234 (not 1.235).
        '''
        print("\n===== Verification =====")
        vmes = ("If the following does not hold, "
                "the results might not be reliable.\n")
        print(vmes)
        ret = self.verify_interval_of_p(
            lower_p, upper_p,
            sig_digits, verbose=2)
        if ret != 0:
            print("\nNG: obtained interval is not reliable, "
                  "try other parameters.")
        else:
            print("\nOK: obtained interval is reliable.")


def verify_interval_of_p(
        params, lower_p, upper_p, sig_digits=-5, verbose=1):
    return params.verify_interval_of_p(lower_p, upper_p, sig_digits, verbose)


def exact(params):
    return params.exact()


# Approximated intervals
def rule_of_ln_alpha(params):
    return params.rule_of_ln_alpha()


def beta_approx(params):
    return params.beta_approx()


def normal_approx(params):
    return params.normal_approx()


def wilson_score(params):
    return params.wilson_score()


def wilson_score_cc(params):
    return params.wilson_score_cc()


def alpha_to_zah(alpha: float, sig_digits: int = -5) -> float:
    """alpha to Z_{alpha/2} for two-sided test in normal distribution.

    Args:
        alpha (float): (1 - (confidence percentage)/100)

    Returns:
        float: zah (Z_{alpha/2})

    Examples:
        >>> alpha_to_zah(0.05)
        1.95996

        >>> alpha_to_zah(0.01)
        2.57583

    Note:
        Not used in v0.0.4 (except in test_for_z()) and will be deprecated.
    """
    return alpha_to_za(alpha / 2., sig_digits=sig_digits)


def alpha_to_za(alpha: float, sig_digits: int = -5) -> float:
    """alpha to Z_alpha for one-sided test in normal distribution.

    Args:
        alpha (float): (1 - (confidence percentage)/100)

    Returns:
        float: za (Z_alpha)

    Examples:
        >>> alpha_to_za(0.05)
        1.64485

        >>> alpha_to_za(0.01)
        2.32635
    """
    if USE_NORMAL_DIST:  # Python >= 3.8
        tmp_zah = NormalDist(mu=0, sigma=1).inv_cdf(1.0 - alpha)
    else:
        tmp_zah = norm.ppf(1.0 - alpha)
    return round_h_up(tmp_zah, sig_digits)


def zah_to_alpha(zah, sig_digits=-5):
    """Z_{alpha/2} to alpha.

    Args:
        zah: Z_{alpha/2} for two-sided test in normal distribution.

    Returns:
        float: alpha, i.e. (1 - (confidence percentage)/100))

    Examples:
        >>> zah_to_alpha(1.95996)
        0.05

        >>> zah_to_alpha(2.57583)
        0.01
    """
    if USE_NORMAL_DIST:  # Python >= 3.8
        cdf_zah = NormalDist(mu=0, sigma=1).cdf(zah)
    else:
        cdf_zah = norm.cdf(zah)
    return round_h_up(2 * (1 - cdf_zah), sig_digits)


def za_to_alpha(za, sig_digits=-5):
    """Z_alpha to alpha.

    Args:
        za: Z_alpha for one-sided test in normal distribution.

    Returns:
        float: alpha, i.e. (1 - (confidence percentage)/100))

    Examples:
        >>> za_to_alpha(1.64485)
        0.05

        >>> za_to_alpha(2.32635)
        0.01
    """
    return round_h_up(1 - norm.cdf(za), sig_digits)


def print_interval(params):
    return params.print_interval()


# ===== Graphs =====
class GraProps:
    """Parameter pool class for graph properties.

    Args:
        k_start: Start range of k.
        k_end: End range of k.
        k_step: Step of k in range(k_start, k_end, k_step).
        log_n_end: Max of x-axis given by (k_end * 10**log_n_end).
        confi_perc_list: List of confidence percentages, e.g. [95, 99].
        savefig: Whether save the figure or not.
        fig_file_name: File name of the saved figure where extension may
            include, eps, jpeg, jpg, pdf, pgf, png, ps,
            raw, rgba, svg, svgz, tif, tiff.
        leg_pos: Position of legend in the figure. Choose from 'auto',
            'upper_right', 'upper_right_nm' (no margin)
            'out_right' (outside right).

    Note:
        It seems that saveifg, i.e. plt.savefig(),
        in jupyter saves blank figure.
        In jupyter, drag and drop the displayed figure to a folder,
        open it with Paint app and so on and then overwrite the figure.
    """

    def __init__(
            self,
            # Defaults
            k_start=0,
            k_end=2,
            k_step=1,
            log_n_end=6,
            confi_perc_list=[99],
            line_list=['with_exact'],
            savefig=False,
            fig_file_name='intervals.png',
            leg_pos='auto',
            dpi=150,
            ):
        self.savefig = savefig
        self.fig_file_name = fig_file_name
        self.leg_pos = leg_pos
        self.confi_perc_list = confi_perc_list
        self.line_list = line_list
        self.k_start = None
        self.k_end = None
        self.k_step = None
        self.k_diff = None
        self.log_n_end = None
        self.set_k(k_start, k_end, k_step)  # may exit
        self.set_n(log_n_end)               # may exit
        # For interval_graph()
        self.dpi = dpi
        self.n_list = None
        self.colorlist = None
        self.col_max_per_k = None
        self.col_offset = None
        self.trans_alpha = None
        self.leg = None

    def set_k(self, k_start, k_end, k_step):
        # Input check
        if (k_start < 0) or (k_end < k_start):
            logger.error(f"'0 <= k_start <= k_end' must hold "
                         f"where k_start={k_start} and k_end={k_end}!!")
            sys.exit(1)
        if k_step <= 0:
            logger.error(f"'k_step > 0' must hold "
                         f"where k_step={k_step}!!")
            sys.exit(1)
        self.k_start = k_start
        self.k_end = k_end
        self.k_step = k_step
        self.k_diff = k_end - k_start

    def set_n(self, log_n_end):
        # Input check
        if (log_n_end < 1):
            logger.error(f"'1 <= log_n_end' must hold "
                         f"where log_n_end={log_n_end}!!")
            sys.exit(1)
        # n <= k_end * 10^log_n_end <= 10**7.
        RELIABLE_LOG_N_MAX = 7
        if (10**RELIABLE_LOG_N_MAX < self.k_end * 10**log_n_end):
            log_n_end_org = log_n_end
            log_n_end = math.floor(
                np.log10((10**RELIABLE_LOG_N_MAX)/(self.k_end)))
            logger.warning(
                f"log_n_end={log_n_end_org} was so large and changed to"
                f" {log_n_end}.")
        self.log_n_end = log_n_end

    def set_graph(
            self, params,
            lower_upper: tuple,
            this_linestyle: str, this_label: str,
            ) -> bool:
        """Set graphs for interval_graph()

        Subroutine to set graph parameters.

        Args:
            params (Params): Instance including k, n, alpha (or confi_perc).
            lower_upper (tuple): tuple containing
                both for 0<k or upper_p only for k=0:

                lower_p (float): lower interval of Binomial confidence

                upper_p (float): upper interval of Binomial confidence

            this_linestyle (str): linestyle,
                such as 'solid', 'dashed', 'dashdot' and 'dotted'.
            this_label (str): label string in the graph legend.

        Returns:
            bool: 0 for no error, 1 otherwise.
        """
        if params.k == 0:
            lower_upper_size = 1  # upper only
        else:
            lower_upper_size = len(lower_upper.shape)

        if (lower_upper_size < 0) or (2 < lower_upper_size):
            logger.error(
                f"lower_upper_size:{lower_upper_size} is out of range")
            return 1
        lower_upper_text = [' lower', ' upper']
        lower_upper_marker = [2, 3]
        for pre_i in range(lower_upper_size):
            i = 1 - pre_i
            # upper, lower for lower_upper_size == 2
            # upper only for lower_upper_size == 1
            plt.plot(
                self.n_list,
                lower_upper[:, i],
                color=self.colorlist[
                    (params.k - self.k_start) * self.col_max_per_k
                    + self.col_offset],
                alpha=self.trans_alpha,
                marker=lower_upper_marker[i],
                linestyle=this_linestyle,
                label=(
                    f"$k$={params.k} {params.confi_perc}% {this_label}"
                    + lower_upper_text[i]))
            self.leg += 1

        if self.k_diff == 0:
            self.col_offset += 1
        return 0


def interval_graph(gra_props):
    """Show interval graphs.

    Coloring rule:

    - Exact lines and the line for k/n for the same k are depicted
        in the same color.

    - If k_start == k_end, approximated lines are depicted in uniq colors.

    - Otherwise, all the lines per the same k are depicted in the same color.

    Args:
        gra_props (GraProps): Instance including k_start, k_end, k_step,
            log_n_end, confi_perc_line, line_list.

    Examples:
        Show intervals for multiple confidence percentages.
        See ebcic.ipynb for more examples.

        >>> interval_graph(GraProps(  # doctest: +SKIP
            k_start=1,    # >= 0
            k_end=1,      # > k_start
            k_step=1,     # > 0
            log_n_end=6,  # max(n) = k_end*10**log_n_end
            confi_perc_list=[90, 95, 99, 99.9, 99.99],
            line_list=[
                'with_exact',
                # The followings are approximation.
                'with_line_kn',  # Line of k/n
                # Available only when k = 0
                # 'with_rule_of_la',
                # Available only when k > 0
                # 'with_normal',
                # 'with_wilson',
                # 'with_wilson_cc',
                # 'with_beta_approx',
            ],
            # savefig = True,
            # fig_file_name='intervals.png',
            # leg_pos='upper_right_nm',
            ))

    Note:
        It seems that saveifg, i.e. plt.savefig(),
        in jupyter saves blank figure.
        In jupyter, drag and drop the displayed figure to a folder,
        open it with Paint app and so on and then overwrite the figure.
    """
    '''
    mul = 2.2
    plt.rcParams['figure.figsize'] = 4 * mul, 3 * mul
    '''
    plt.rcParams['figure.dpi'] = 150
    # Customize of graph
    k_start = gra_props.k_start
    k_end = gra_props.k_end
    k_step = gra_props.k_step
    k_diff = gra_props.k_diff
    log_n_end = gra_props.log_n_end
    line_list = gra_props.line_list
    leg_pos = gra_props.leg_pos
    # confi_perc = params.confi_perc

    title = (
        # f'Interval with {confi_perc}% confidence '
        'Interval of $p$ '
        'after $k$ errors are observed among $n$ trials')
    x_label = 'Number of Trials ($n$)'
    y_label = 'Error Rate ($p$)'
    gra_props.trans_alpha = 0.6  # Transparency

    # Set colorlist
    # With other lines
    line_types_per_k = 1  # Exact interval
    if k_diff == 0:
        # commmon
        if 'with_wilson' in line_list:
            line_types_per_k += 1
        if 'with_wilson_cc' in line_list:
            line_types_per_k += 1
        if 'with_beta_approx' in line_list:
            line_types_per_k += 1
        # k specific
        if k_start == 0:
            if 'with_rule_of_la' in line_list:
                line_types_per_k += 1
        else:  # k_start > 0
            if 'with_normal' in line_list:
                line_types_per_k += 1

    # Max num of colors
    gra_props.col_max_per_k = line_types_per_k * len(gra_props.confi_perc_list)
    col_max = (
        gra_props.k_end - gra_props.k_start + 1) * gra_props.col_max_per_k
    cm = plt.get_cmap('brg')
    # yellow lines are haed to see.
    # cm = plt.get_cmap('gist_rainbow')
    # cm = plt.get_cmap('rainbow')
    cNorm = colors.Normalize(vmin=0, vmax=col_max - 1)
    scalarMap = mplcm.ScalarMappable(norm=cNorm, cmap=cm)
    gra_props.colorlist = [scalarMap.to_rgba(i) for i in range(col_max)]

    n_list_base = (
        sorted(
            [10**i for i in range(0, log_n_end + 1, 1)]
            # intermediates where 3 is floor(sqrt(10))
            + [3 * 10**i for i in range(0, log_n_end, 1)]
            + [2]))
    # print(n_list_base)

    params = Params()
    gra_props.leg = 0  # Num of legend
    for params.k in range(k_start, k_end + 1, k_step):
        gra_props.n_list = np.array(n_list_base)
        if params.k > 1:
            gra_props.n_list *= params.k

        if ('with_line_kn' in line_list) and (params.k > 0):
            # Line of k/n
            k_expect = [params.k / n for n in gra_props.n_list]
            plt.plot(gra_props.n_list,
                     k_expect,
                     color=gra_props.colorlist[(
                        params.k - gra_props.k_start)
                        * gra_props.col_max_per_k],
                     alpha=gra_props.trans_alpha,
                     linestyle='dashed',
                     label=f"$k$={params.k} k/n={params.k}/n")
            gra_props.leg += 1

        gra_props.col_offset = 0
        for confi_perc in gra_props.confi_perc_list:
            params.set_confi_perc(confi_perc)
            # Exact lower and upper
            if 'with_exact' in line_list:
                gra_props.set_graph(
                    params,
                    lower_upper=np.array(
                        [exact(params) for params.n in gra_props.n_list]),
                    this_linestyle='solid',
                    this_label="exact",
                    )
            # Approximated Intervals
            if params.k == 0:
                if 'with_rule_of_la' in line_list:
                    gra_props.set_graph(
                        params,
                        lower_upper=np.array(
                            [rule_of_ln_alpha(params)
                             for params.n in gra_props.n_list]),
                        this_linestyle='dashdot',
                        this_label="rule of -ln(a)",
                        )

            if 0 < params.k:
                if 'with_normal' in line_list:
                    # Normal
                    gra_props.set_graph(
                        params,
                        lower_upper=np.array(
                            [normal_approx(params)
                             for params.n in gra_props.n_list]),
                        this_linestyle='dashdot',
                        this_label="normal",
                        )

            if 'with_wilson' in line_list:
                # Wilson
                gra_props.set_graph(
                    params,
                    lower_upper=np.array(
                        [wilson_score(params)
                            for params.n in gra_props.n_list]),
                    this_linestyle='dashdot',
                    this_label="Wilson",
                    )

            if 'with_wilson_cc' in line_list:
                # Wilson with continuity correction
                gra_props.set_graph(
                    params,
                    lower_upper=np.array(
                        [wilson_score_cc(params)
                         for params.n in gra_props.n_list]),
                    this_linestyle='dashdot',
                    this_label="Wilson cc",
                    )

            if 'with_beta_approx' in line_list:
                # Approximation using beta distribution
                gra_props.set_graph(
                    params,
                    lower_upper=np.array(
                        [beta_approx(params)
                            for params.n in gra_props.n_list]),
                    this_linestyle='dashed',
                    this_label="beta approx",
                    )

    # Show
    # plt.figure(figsize=(400, 300), dpi=300)
    # plt.rcParams["figure.figsize"]=25,20
    plt.title(title)
    plt.xlabel(x_label, fontsize=12)
    plt.ylabel(y_label, fontsize=12)

    leg_fontsize = 8
    # print(leg_pos, (leg_pos == 'auto'))
    if (gra_props.leg <= 5 and (leg_pos == 'auto')) or (
             leg_pos == 'upper_right'):
        lgd = plt.legend(
            bbox_to_anchor=(1, 1),
            loc='upper right',
            borderaxespad=1,
            fontsize=leg_fontsize)
    elif (gra_props.leg <= 6 and (leg_pos == 'auto')) or (
            leg_pos == 'upper_right_nm'):
        lgd = plt.legend(
            # bbox_to_anchor=(1, 1),
            bbox_to_anchor=(0.99, 0.99),
            loc='upper right',
            borderaxespad=0,
            fontsize=leg_fontsize)
    elif (leg_pos == 'auto') or (leg_pos == 'out_right'):
        lgd = plt.legend(
            # bbox_to_anchor=(1.05, 1),
            bbox_to_anchor=(1.01, 1),
            loc='upper left',
            borderaxespad=0,
            fontsize=leg_fontsize)
    plt.grid(which='major', color='black', linestyle='-')
    plt.grid(which='minor', color='black', linestyle='dotted')
    plt.yscale('log')
    plt.xscale('log')
    plt.show()

    if gra_props.savefig:
        plt.savefig(
            gra_props.fig_file_name,
            dpi=gra_props.dpi,
            bbox_extra_artists=(lgd,),
            bbox_inches='tight',)

    plt.cla()
    plt.clf()


def compare_dist(params):
    """Compare major probability distributions in a graph.

    Args:
        params (Params): Instance including k, n and zah (Z_{alpha/2})
            where zah is set by giving either 'alpha' or 'confi_perc'
            to Params().

    Examples:
        >>> compare_dist(Params(k=2, n=20))  # doctest: +SKIP
    """
    # mul = 1
    # plt.rcParams["figure.figsize"] = 4 * mul, 3 * mul
    plt.rcParams['figure.dpi'] = 150
    # plt.rcParams['figure.dpi'] = 100

    # Settings
    k = params.k
    n = params.n

    p = k / n
    x_max = max(8, 3 * k)

    title = (
        f'Probability distribution of $k$ errors in {n} trials'
        f' for $p={k/n:.3g}$')
    x_label = 'Number of Errors ($k$)'
    y_label = 'Probability'

    # Distribution
    x = range(0, x_max + 1, 1)
    # Binomial
    plt.plot(x, binom.pmf(x, n, p), label="Binomial")

    # Poisson
    plt.plot(x, poisson.pmf(x, k), label="Poisson")

    # Normal distribution
    variance = n * p * (1 - p)
    sigma = np.sqrt(variance)
    x = range(-3, x_max + 1, 1)
    plt.plot(x, norm.pdf(x, k, sigma), label="Normal")

    # Plot
    plt.title(title)
    plt.xlabel(x_label, fontsize=12)
    plt.ylabel(y_label, fontsize=12)

    plt.legend(
        bbox_to_anchor=(1, 1),
        loc='upper right',
        borderaxespad=1,
        fontsize=10)

    plt.grid(which='major', color='black', linestyle='dotted')
    # plt.grid(which='minor', color='black', linestyle='dotted')


# ===== Tests =====
def test_for_z():
    """Self test regarding Z_{alpha/2} and Z_{alpha}.

    alpha_to_za(), za_to_alpha(),
    alpha_to_zah() and zah_to_alpha() are tested.

    Returns:
        bool: 0 for no error, 1 otherwise.
    """
    # Right patterns
    right_alpha_zah = [[0.05, 1.95996], [0.01, 2.57583]]
    right_alpha_za = [[0.05, 1.64485], [0.01, 2.32635]]

    print("=== Test for Z_{alpha} and Z_{alpha/2} ===")

    ret = 0
    print("\n== zah_to_alpha() and alpha_to_zah() ==")
    for alpha, zah in right_alpha_zah:
        print(f"Check if zah_to_alpha({zah}) = {zah_to_alpha(zah)} = {alpha}")
        if zah_to_alpha(zah) != alpha:
            logger.error(f"{zah_to_alpha(zah)} = {alpha}")
            ret = 1
        print(
            f"Check if alpha_to_zah({alpha}) = {alpha_to_zah(alpha)} = {zah}")
        if alpha_to_zah(alpha) != zah:
            logger.error(f"{alpha_to_zah(alpha)} = {zah}")
            ret = 1

    print("\n== za_to_alpha() and alpha_to_za() ==")
    for alpha, za in right_alpha_za:
        print(f"Check if za_to_alpha({za}) = {za_to_alpha(za)} = {alpha}")
        if za_to_alpha(za) != alpha:
            logger.error(f"{za_to_alpha(za)} = {alpha}")
            ret = 1
        print(f"Check if alpha_to_za({alpha}) = {alpha_to_za(alpha)} = {za}")
        if alpha_to_za(alpha) != za:
            logger.error(f"{alpha_to_za(alpha)} = {za}")
            ret = 1

    if ret == 0:
        print("\n== OK : Test for Z_{alpha} and Z_{alpha/2} ==\n")
    else:
        print("\n== NG : Test for Z_{alpha} and Z_{alpha/2} ==\n")
    return ret


def test_warning_once():
    """Self test of warning once.

    Test of warning once for the same setting.

    Returns:
        bool: 0 for no error, 1 otherwise.
    """
    ret = 0
    params = Params(warn_line_confi_perc=80)
    print("==== Test of no duplicate warning ====")
    print("Test is success, if one WARNING in each section is displayed.")

    print("\n=== params.set_confi_perc() -> params.set_alpha() ===")
    confi_perc = 70
    cou = 0
    cou += params.set_confi_perc(confi_perc)  # Should be warned
    cou += params.set_confi_perc(confi_perc)
    alpha = 0.3
    cou += params.set_alpha(alpha)
    if cou != 1:
        ret = 1
        logger.error("duplicate warnings!!")

    print("\n=== params.set_alpha() ===")
    alpha = 0.4
    cou = 0
    cou += params.set_alpha(alpha)  # Should be warned
    cou += params.set_alpha(alpha)
    confi_perc = 60
    cou += params.set_confi_perc(confi_perc)
    if cou != 1:
        ret = 1
        logger.error("duplicate warnings!!")

    print("\n=== Input check of alpha ===")
    msg_wo = "\n= Should be warned once. ="
    msg_nw = "\n= Should not be warned. ="

    print("\n== params.alpha_to_confi_perc(alpha) -> "
          "params.confi_perc_to_alpha(alpha_to_confi_perc_wo_check(alpha)) ==")
    alpha = 0.5
    params.num_of_warned = 0
    if alpha != params.alpha:
        print(msg_wo)
        right_num_of_warned = 1
    else:
        print(msg_nw)
        right_num_of_warned = 0
    print(f"alpha={alpha}, "
          f"checked one={params.alpha}.")
    print(f"params.alpha_to_confi_perc({alpha}) = "
          f"{params.alpha_to_confi_perc(alpha)}")
    print(f"params.alpha_to_confi_perc({alpha}) = "
          f"{params.alpha_to_confi_perc(alpha)}")
    tmp = params.confi_perc_to_alpha(
        alpha_to_confi_perc_wo_check(alpha))
    print("params.confi_perc_to_alpha(alpha_to_confi_perc_wo_check"
          f"({alpha})) = {tmp}")
    tmp = params.confi_perc_to_alpha(
        alpha_to_confi_perc_wo_check(alpha))
    print("params.confi_perc_to_alpha(alpha_to_confi_perc_wo_check"
          f"({alpha})) = {tmp}")
    if params.num_of_warned != right_num_of_warned:
        ret = 1
        logger.error("duplicate warnings!!")

    print(
        "\n== params.confi_perc_to_alpha(alpha_to_confi_perc_wo_check(alpha))"
        " -> params.alpha_to_confi_perc(alpha) ==")
    alpha = 0.6
    params.num_of_warned = 0
    if alpha != params.alpha:
        print(msg_wo)
        right_num_of_warned = 1
    else:
        print(msg_nw)
        right_num_of_warned = 0
    print(f"alpha={alpha}, "
          f"checked one={params.alpha}.")
    tmp = params.confi_perc_to_alpha(
        alpha_to_confi_perc_wo_check(alpha))
    print(
        "params.confi_perc_to_alpha(alpha_to_confi_perc_wo_check"
        f"({alpha})) = {tmp}")
    tmp = params.confi_perc_to_alpha(
        alpha_to_confi_perc_wo_check(alpha))
    print(
        "params.confi_perc_to_alpha(alpha_to_confi_perc_wo_check"
        f"({alpha})) = {tmp}")
    print(f"params.alpha_to_confi_perc({alpha}) = "
          f"{params.alpha_to_confi_perc(alpha)}")
    print(f"params.alpha_to_confi_perc({alpha}) = "
          f"{params.alpha_to_confi_perc(alpha)}")
    if params.num_of_warned != right_num_of_warned:
        ret = 1
        logger.error("duplicate warnings!!")

    print("\n=== Input check of confi_perc ===")

    print("\n== params.confi_perc_to_alpha(confi_perc) -> "
          "params.alpha_to_confi_perc("
          "confi_perc_to_alpha_wo_check(confi_perc)) ==")
    confi_perc = 55
    params.num_of_warned = 0
    if confi_perc != params.confi_perc:
        print(msg_wo)
        right_num_of_warned = 1
    else:
        print(msg_nw)
        right_num_of_warned = 0
    print(f"confi_perc={confi_perc}, "
          f"checked one={params.confi_perc}.")
    print(f"params.confi_perc_to_alpha({confi_perc}) = "
          f"{params.confi_perc_to_alpha(confi_perc)}")
    print(f"params.confi_perc_to_alpha({confi_perc}) = "
          f"{params.confi_perc_to_alpha(confi_perc)}")
    tmp = params.alpha_to_confi_perc(
        confi_perc_to_alpha_wo_check(confi_perc))
    print("params.alpha_to_confi_perc(confi_perc_to_alpha_wo_check"
          f"({confi_perc})) = {tmp}")
    tmp = params.alpha_to_confi_perc(
        confi_perc_to_alpha_wo_check(confi_perc))
    print("params.alpha_to_confi_perc(confi_perc_to_alpha_wo_check"
          f"({confi_perc})) = {tmp}")
    if params.num_of_warned != right_num_of_warned:
        ret = 1
        logger.error("duplicate warnings!!")

    print("\n== params.alpha_to_confi_perc("
          "confi_perc_to_alpha_wo_check(confi_perc)) -> "
          "params.confi_perc_to_alpha(confi_perc) ==")
    confi_perc = 65
    params.num_of_warned = 0
    if confi_perc != params.confi_perc:
        print(msg_wo)
    else:
        print(msg_nw)
    print(f"alpha={confi_perc}, "
          f"checked one={params.confi_perc}.")
    tmp = params.alpha_to_confi_perc(
            confi_perc_to_alpha_wo_check(confi_perc))
    print(
        "params.alpha_to_confi_perc(confi_perc_to_alpha_wo_check"
        f"({confi_perc})) = {tmp}")
    tmp = params.alpha_to_confi_perc(
            confi_perc_to_alpha_wo_check(confi_perc))
    print(
        "params.alpha_to_confi_perc(confi_perc_to_alpha_wo_check"
        f"({confi_perc})) = {tmp}")
    print(f"params.confi_perc_to_alpha({confi_perc}) = "
          f"{params.confi_perc_to_alpha(confi_perc)}")
    print(f"params.confi_perc_to_alpha({confi_perc}) = "
          f"{params.confi_perc_to_alpha(confi_perc)}")
    if params.num_of_warned != right_num_of_warned:
        ret = 1
        logger.error("duplicate warnings!!")

    if ret == 0:
        print("\n== OK : No duplicate warning ==\n")
    else:
        print("\n== NG : No duplicate warning ==\n")
    return ret


def test_of_intervals(
        n_start,
        n_end,
        n_step,
        confi_perc,
        sig_digits=-5):
    """Test of reliability of obtained intervals for given parameters.

    Check if 'alpha' is equal to the integral of the outside of the interval
    for given parameters, such as k and n, and then return 0 if it holds.

    Returns:
        int: num_of_wrongs, the number of unreliable intervals.

    Examples:
        See ebcic.ipynb for the parameters where reliable intervals
        might not be able to be obtained.

        >>> test_of_intervals(  # doctest: +SKIP
            n_start=1,
            n_end=1000,
            n_step=1,
            confi_perc=99.0,  # Confidence Percentage (90 <= confi_perc < 100)
            sig_digits=-5)
    """
    warnings.filterwarnings(
        'ignore', 'The iteration is not making good progress')

    check_n(n_start)
    check_n(n_end)
    params = Params(confi_perc=confi_perc)
    '''
    # ALTERNATIVE:
    # Check k selectively.
    K_CHANGE_LINE = 500
    PARTIAL_RANGE = 100  # < K_CHANGE_LINE/3
    '''
    num_of_wrongs = 0
    print(f"Checking from n = {n_start} with {n_step} step(s)")
    for n in range(n_start, n_end + 1, n_step):
        if (n % 100) == 0:
            print(f"Have checked n < {n}")
        '''
        # ALTERNATIVE:
        # Check k selectively.
        if n < K_CHANGE_LINE:
            k_list = range(0, n + 1, 1)
        else:
            k_list = list(range(0, PARTIAL_RANGE, 1))
            k_list.extend(
                range(int((n - PARTIAL_RANGE)/2), int(
                    (n + PARTIAL_RANGE)/2), 1))
            k_list.extend(range(n - (PARTIAL_RANGE + 1), n + 1, 1))
        for k in k_list:
        '''
        # Check all of the k.
        for k in range(0, n + 1, 1):
            if n > 100000:
                if (k % 1000000) == 0:
                    print(f"Have checked k < {k}")
            # NOTE: use params.set_*()
            #       unless range of the parameter is 100% sure.
            # e.g.
            # params.set_n(n)
            # params.set_k(k)
            params.n = n
            params.k = k
            lower_p, upper_p = exact(params)
            ret = verify_interval_of_p(
                # k, n, alpha, lower_p, upper_p, sig_digits, verbose=1)
                params, lower_p, upper_p, sig_digits, verbose=1)
            if ret != 0:
                print(f"Wrong interval at n={n} k={k}")
                num_of_wrongs += 1
    print(f"Have checked up to n = {n}")
    if num_of_wrongs == 0:
        print("=== Test OK ===")
    return num_of_wrongs


def test_all():
    """Test of all or some of the above tests.

    Higher 'ret' is more critical.
    """
    ret = 0
    tmp_ret = test_for_z()
    if tmp_ret > 0:
        ret = tmp_ret

    tmp_ret = test_warning_once()
    if tmp_ret > 0:
        ret = tmp_ret

    num_of_wrongs = test_of_intervals(
            n_start=1,
            n_end=10,
            n_step=1,
            confi_perc=99.0,  # Confidence Percentage (90 <= confi_perc < 100)
            sig_digits=-5)
    if num_of_wrongs > 0:
        ret = num_of_wrongs

    num_of_wrongs = test_of_intervals(
            n_start=11,
            n_end=400,
            n_step=101,
            confi_perc=95.0,  # Confidence Percentage (90 <= confi_perc < 100)
            sig_digits=-5)
    if num_of_wrongs > 0:
        ret = num_of_wrongs

    num_of_wrongs = test_of_intervals(
            n_start=1000,
            n_end=1000,
            n_step=1,
            confi_perc=90.0,  # Confidence Percentage (90 <= confi_perc < 100)
            sig_digits=-5)
    if num_of_wrongs > 0:
        ret = num_of_wrongs

    if ret == 0:
        print("\n== OK : All the tests are succeeded. ==\n")
    else:
        print("\n== NG : Some tests are failed!! ==\n")
    return ret
