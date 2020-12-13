# %load ebcic.py
# %%writefile ebcic.py
# To be syncronized with ebcic.py copy either of the above magic functions
# depending on whether to overwrite or to load, and then run this cell.
"""Exact Binomial Interval Calculator

This module provides functions to calculate exact binomial intervals
and some of the approximated intervals, such as normal approximation,
generalized rule of three, Wilson score intervals and so on.

This module also provides functions to depict their intervals.


When you use or publish the confidence interval obtained with the software,
please **refer to the program name, version, platform,** etc, so that
readers can verify the correctness and reproducibility of the interval
with the parameters, n, k and the confidence level.

Example of the reference is:

**"This confidence interval is obtained by EBCIC 0.0.0 on Python 3."**

The initial software is based on results obtained from a project,
JPNP16007, commissioned by
the New Energy and Industrial Technology Development Organization (NEDO).

THE SOFTWARE COMES WITH ABSOLUTELY NO WARRANTY.

Copyright (C) 2020 National Institute of Advanced Industrial Science
and Technology (AIST).
"""
__version__ = '0.0.0'
import os
import sys
import logging
import numpy as np
import scipy
import math
import matplotlib
import matplotlib.cm as mplcm
import matplotlib.colors as colors
import warnings
from numpy import log
from scipy import optimize
from scipy import stats as st
from scipy.stats import binom
from matplotlib import pyplot as plt
from platform import python_version
from decimal import Decimal, ROUND_HALF_UP

print(f"python    : {python_version()}")  # os, math, warnings
print(f"np        : {np.__version__}")
print(f"sys       : {sys.version[:6]}")
print(f"logging   : {logging.__version__}")
print(f"scipy     : {scipy.version.full_version}")
print(f"matplotlib: {matplotlib.__version__}")

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
        123.4

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
    if (confi_perc <= 0) or (100 <= confi_perc):
        # raise ValueError("'0 < confi_perc={confi_perc} < 100' must hold")
        logger.error(
            f"'0 < confi_perc < 100' must hold where "
            f"confi_perc={confi_perc}!!")
        sys.exit(1)


def check_alpha(alpha):
    if (alpha <= 0) or (1 <= alpha):
        logger.error(f"'0 < alpha < 1' must hold where alpha={alpha}!!")
        sys.exit(1)


class Params:
    """Parameter pool class for Binomial and confidence.

    This stores parameters regarding Binomial and confidence after
    checking them and avoiding duplicate checks and warnings.
    """

    def __init__(
            self,
            k=None,
            n=None,
            alpha=None,
            confi_perc=None,
            # FIDO Biometrics Requirements v1.1
            # https://fidoalliance.org/specs/biometric/requirements/
            # Biometrics-Requirements-v1.1-fd-20190606.pdf
            # uses 80% for one-sided (60% for two-sided) confidence.
            # (confi_perc < warn_line_confi_perc) is warned
            warn_line_confi_perc=60):  # for two-sided confidence
        self.n = None
        self.k = None
        self.alpha = None
        self.confi_perc = None
        self.za = None
        self.zah = None
        # 'warn_line' shall be set by warn_line_confi_perc even for alpha.
        check_confi_perc(warn_line_confi_perc)  # may exit
        self.warn_line_confi_perc = warn_line_confi_perc
        self.warn_line_alpha = confi_perc_to_alpha_wo_check(
            warn_line_confi_perc)
        self.num_of_warned = 0  # For test
        # Check and set parameters.
        # 'confi_perc' or 'alpha'
        which_to_be_used = 'confi_perc'
        if (confi_perc is None) and (alpha is None):
            which_to_be_used = None
        elif (confi_perc is not None) and (alpha is not None):
            logger.warning("giving either 'alpha' or "
                           "'confi_perc' is enough.")
            if alpha != confi_perc_to_alpha_wo_check(warn_line_confi_perc):
                logger.error(
                    f"'alpha [{alpha}] != "
                    f"confi_perc_to_alpha_wo_check(confi_perc "
                    f"[confi_perc])"
                    f"[{confi_perc_to_alpha_wo_check(confi_perc)}]'!!")
                sys.exit(1)
        else:  # (either confi_perc) or (alpha is not None)
            if alpha is not None:
                which_to_be_used = 'alpha'
        # Check and set both 'confi_perc' and 'alpha'.
        if which_to_be_used is not None:
            if which_to_be_used == 'confi_perc':
                self.set_confi_perc(confi_perc)  # alpha is also set
            elif which_to_be_used == 'alpha':
                self.set_alpha(alpha)       # confi_perc is also set
            else:
                logger.error(f"unknown which_to_be_used={which_to_be_used}!!")
        # Set k and n
        if n is not None:
            self.set_n(n)
        if k is not None:
            self.set_k(k)

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
                f"'confi_perc < {self.warn_line_confi_perc}'"
                " is unusual.\n"
                f"Is 'confi_perc={confi_perc}' correct?")
            warned = 1
            self.num_of_warned += 1
        else:
            warned = 0
        return warned

    def set_confi_perc(self, confi_perc):
        """Set confidence percentage and more.

        This also sets its corresponding alpha, Z_{alpha} and Z_{alpha/2}.
        """
        # See if ever checked, otherwise check.
        tmp_alpha = confi_perc_to_alpha_wo_check(confi_perc)
        tmp_zah = alpha_to_zah(tmp_alpha)
        tmp_za = alpha_to_za(tmp_alpha)
        warned = 0
        if (
                (self.confi_perc != confi_perc)
                or (self.alpha != tmp_alpha)
                or (self.zah != tmp_zah)
                or (self.za != tmp_za)):
            # Check
            warned = self.check_and_warn_confi_perc(confi_perc)
            # set
            # self.alpha_or_confi_perc_checked = True
            self.confi_perc = confi_perc
            # self.alpha_checked = True
            self.alpha = tmp_alpha
            self.zah = tmp_zah
            self.za = tmp_za
        return warned

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
        'perc_confi' (confidence percentage),
        Z_{alpha} and Z_{alpha/2}.
        """

        # See if ever checked, otherwise check.
        tmp_confi_perc = alpha_to_confi_perc_wo_check(alpha)
        tmp_zah = alpha_to_zah(alpha)
        tmp_za = alpha_to_za(alpha)
        # print(f"self.zah != tmp_zah : {self.zah} != {tmp_zah}")
        # print(f"self.za != tmp_za : {self.za} != {tmp_za}")
        warned = 0
        if (
                (self.alpha != alpha)
                or (self.confi_perc != tmp_confi_perc)
                or (self.zah != tmp_zah)
                or (self.za != tmp_za)
        ):
            # Check
            warned = self.check_and_warn_alpha(alpha)
            # Set
            self.alpha = alpha
            # print(self.confi_perc)
            self.confi_perc = tmp_confi_perc
            self.zah = tmp_zah
            self.za = tmp_za
        return warned

    def confi_perc_to_alpha(self, confi_perc):
        """Confidence percentage to alpha.

        This converts 'confi_perc' given in percentage between
        0 to 100 exclusive
        to alpha in between 0 to 1 exclusive after input check.
        """
        self.set_confi_perc(confi_perc)
        alpha = confi_perc_to_alpha_wo_check(confi_perc)
        return alpha

    def alpha_to_confi_perc(self, alpha):
        """Alpha to confidence percentage.

        Convert given 'alpha' in between 0 to 1 exclusive
        to 'confi_perc' in percentage between 0 to 100 exclusive
        after input check.
        """
        self.set_alpha(alpha)
        confi_perc = alpha_to_confi_perc_wo_check(alpha)
        return confi_perc


def exact(params):
    """Return exact Binomial confidence interval for the parameter.

    Args:
        params (Params): Instance including k, n and
            alpha (or confi_perc).

    Returns:
        tuple: tuple containing:

            lower_p (float): lower bound of Binomial confidence interval.

            upper_p (float): upper bound of Binomial confidence interval.
    """
    n = params.n
    k = params.k
    alpha = params.alpha

    r = alpha / 2
    reverse_mode = False
    # Cumulative error becomes large for k >> n/2,
    # so make k < n/2.
    if k > n / 2:
        reverse_mode = True
        k = n - k

    def upper(p):
        """
        Upper bound is the p making the following 0.
        """
        return binom.cdf(k, n, p) - r

    def lower(p):
        """
        Lower bound for k>0 is the p making the following 0.
        """
        return binom.cdf(n - k, n, 1 - p) - r
        # Alternatives:
        # return 1-sum([binom.pmf(i,n,result_lower) for i in range(k)]) - r
        # return sum([binom.pmf(i,n,result_lower) for i in range(k,n)]) - r

    if k == 0 or k == n:
        tmp = (alpha)**(1 / n)
        if k == 0:
            lower_p = 0.0
            upper_p = 1 - tmp
        else:
            # k == n:
            lower_p = tmp
            upper_p = 1.0
    else:
        # 0 < k < n
        tmp = k / n
        u_init = tmp
        l_init = tmp
        if k == 1:
            l_init = 0
        elif k == 2:
            l_init = k / (2 * n)
        # print(f"l_init={l_init}, u_init={u_init}")
        lower_p = optimize.fsolve(lower, l_init)[0]
        if n == 2 and k == 1:
            # Exception of k/n = 1/2 and n is too small.
            upper_p = 1 - lower_p
        else:
            upper_p = optimize.fsolve(upper, u_init)[0]

    if reverse_mode:
        tmp = lower_p
        lower_p = 1 - upper_p
        upper_p = 1 - tmp

    return lower_p, upper_p


def verify_interval_of_p(
        params,
        lower_p, upper_p,
        sig_digits=-5, verbose=1):
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
    n = params.n
    k = params.k
    alpha = params.alpha

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
                # print(f"Upper : check if {round_h_up(alpha, sig_digits)} == "
                #      f"{round_h_up(binom.cdf(k, n, upper_p), sig_digits)}")
        else:
            # k == n
            if abs(
                    binom.cdf(n - k, n, 1 - lower_p) -
                    alpha) > acceptable_range:
                ret = 1
            if (verbose == 1 and ret == 1) or verbose >= 2:
                print(
                    f"Lower : check if abs(alpha ({alpha}) - "
                    f"cumulative error ({binom.cdf(n - k, n, 1- lower_p)})) "
                    f"<= {acceptable_range}")
                # print(f"Lower : check if {round_h_up(alpha, sig_digits)} == "
                # f"{round_h_up(binom.cdf(n - k, n, 1 - lower_p),
                # sig_digits)}")
    else:
        # 0 < k < n
        r = alpha / 2
        tmp_ret = 0
        if abs(binom.cdf(k, n, upper_p) - r) > acceptable_range:
            ret = 1
            tmp_ret = 1
        if (verbose == 1 and tmp_ret == 1) or verbose >= 2:
            print(
                f"Upper : check if abs(alpha/2 ({r}) - "
                f"cumulative error ({binom.cdf(k, n, upper_p)})) <= "
                f"{acceptable_range}")
            # print(f"Upper : check if {round_h_up(r, sig_digits)} == "
            #      f"{round_h_up(binom.cdf(k, n, upper_p), sig_digits)}")
        tmp_ret = 0
        if abs(binom.cdf(n - k, n, 1 - lower_p) - r) > acceptable_range:
            ret = 1
            tmp_ret = 1
        if (verbose == 1 and tmp_ret == 1) or verbose >= 2:
            print(
                f"Lower : check if abs(alpha/2 ({r}) - "
                f"cumulative error ({binom.cdf(n - k, n, 1 - lower_p)})) <= "
                f"{acceptable_range}")
            # print(f"Lower : check if {round_h_up(r, sig_digits)} == "
            # f"{round_h_up(binom.cdf(n - k, n, 1 - lower_p), sig_digits)}")

    # Test of position of the intervals
    pp = k / n
    tmp_ret = 0
    if upper_p < pp:
        ret = 1
        tmp_ret = 1
    if (verbose == 1 and tmp_ret == 1) or verbose >= 2:
        print(f"Upper : check if upper bound ({upper_p}) >= k/n ({pp})")
    tmp_ret = 0
    if lower_p > pp:
        ret = 1
        tmp_ret = 1
    if (verbose == 1 and tmp_ret == 1) or verbose >= 2:
        print(f"Lower : check if lower bound ({lower_p}) <= k/n ({pp})")

    return ret


# Approximated intervals
def rule_of_ln_alpha(params):
    """ Generalized rule of three.

    Interval of rule of -ln(alpha), i.e generalized version of
    'rule of three'.
    This can be used only when (k == 0 or k == n) and
    reliable if n is large enough, say n > 50.

    Args:
        params (Params): Instance including k, n, alpha (or confi_perc).

    Returns:
        tuple: tuple containing:

            lower_p (float): lower bound of Binomial confidence interval.

            upper_p (float): upper bound of Binomial confidence interval.
    """
    n = params.n
    k = params.k
    alpha = params.alpha
    # print(f"alpha={alpha}")
    # Get interval
    lower_p = None
    upper_p = None
    if k == 0:
        lower_p = 0
        upper_p = -log(alpha) / n
    elif k == n:
        lower_p = 1 - (-log(alpha) / n)
        upper_p = 1
    else:  # 0 < k < n
        # raise ValueError("either k==0 or k==n must hold!!")
        logger.error(f"either k==0 or k==n must hold where k={k} and n={n}!!")
        sys.exit(1)
    return lower_p, upper_p


def alpha_to_zah(alpha, sig_digits=-5):
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
    """
    USE_NORMAL_DIST = True
    if USE_NORMAL_DIST:  # Python > 3.8
        from statistics import NormalDist
        tmp_zah = NormalDist(mu=0, sigma=1).inv_cdf(1 - alpha / 2.)
    else:
        tmp_zah = st.norm.ppf(1 - alpha / 2.)
    return round_h_up(tmp_zah, sig_digits)


def alpha_to_za(alpha, sig_digits=-5):
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
    return round_h_up(st.norm.ppf(1 - alpha), sig_digits)


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
    USE_NORMAL_DIST = True
    if USE_NORMAL_DIST:
        from statistics import NormalDist
        cdf_zah = NormalDist(mu=0, sigma=1).cdf(zah)
    else:
        cdf_zah = st.norm.cdf(zah)
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
    return round_h_up(1 - st.norm.cdf(za), sig_digits)


def normal_approx(params):
    """Approximated interval using normal distribution.

    Interval obtained by approximating Binomial distribution
    to normal distribution.

    Args:
        params (Params): Instance including k, n and zah (Z_{alpha/2})
            where zah is set by giving either 'alpha' or 'confi_perc'
            to Params().

    Returns:
        tuple: tuple containing:

            lower_p (float): lower bound of Binomial confidence interval.

            upper_p (float): upper bound of Binomial confidence interval.

    Note:
        Normal approximation does not give good approximation for small
        n or k being close to 0 or n.
    """
    n = params.n
    k = params.k
    zah = params.zah

    p = k / n
    # sigma = np.sqrt(p * (1 - p)/n)
    if n != 1:
        sigma = np.sqrt(n * p * (1 - p))/(n - 1)
    else:
        sigma = np.sqrt(p * (1 - p))
    half_width = zah * sigma
    lower_p = max(0, p - half_width)
    upper_p = min(1, p + half_width)
    return lower_p, upper_p


def wilson_score(params):
    """Wilson score interval.

    Args:
        params (Params): Instance including k, n and zah (Z_{alpha/2})
            where zah is set by giving either 'alpha' or 'confi_perc'
            to Params().

    Returns:
        tuple: tuple containing:

            lower_p (float): lower bound of Binomial confidence interval.

            upper_p (float): upper bound of Binomial confidence interval.

    Note:
        Upper interval seems close to the exact one, but
        lower one seems a little greater than the exact one.
    """
    n = params.n
    k = params.k
    z = params.zah
    p = k / n
    mu = 2 * k + z**2
    half_width = z * np.sqrt(z**2 + 4 * k * (1 - p))
    denomi = 2 * (n + z**2)
    lower_p = max(0, (mu - half_width) / denomi)
    upper_p = min(1, (mu + half_width) / denomi)
    return lower_p, upper_p


def wilson_score_cc(params):
    """Wilson score interval with continuity correction.

    Args:
        params (Params): Instance including k, n and zah (Z_{alpha/2})
            where zah is set by giving either 'alpha' or 'confi_perc'
            to Params().

    Returns:
        tuple: tuple containing:

            lower_p (float): lower bound of Binomial confidence interval.

            upper_p (float): upper bound of Binomial confidence interval.

    Note:
        Lower interval is not a good approximation for small k,
        i.e. it becomes higher than the exact lower interval,
        though it is better than Wilson score interval
        (without continuity correction).
    """
    n = params.n
    k = params.k
    z = params.zah
    p = k / n
    mu = 2 * k + z**2
    half_width = 1 + z * np.sqrt(z**2 - 1 / n + 4 * k * (1 - p) + (4 * p - 2))
    denomi = 2 * (n + z**2)
    lower_p = max(0, (mu - half_width) / denomi)
    upper_p = min(1, (mu + half_width) / denomi)
    return lower_p, upper_p


# ===== CUI =====
def print_interval(params):
    """Print exact interval as text.

    Args:
        params (Params): Instance including k, n and
            'confi_perc' (confidence percentage) or 'alpha'.

    Examples:
        >>> print_interval(Params(k=1, n=10000, confi_perc = 95.0))
    """
    # Num of significant digits.
    # Negative means below decimal point.
    # sig_digits = -10
    sig_digits = -14  # To show lower interval around 10^-8.

    # Print Parameters
    print("===== Parameters =====")
    print("n (num of trials)  :", params.n)
    print("k (num of failures):", params.k)
    print("k/n (observed p')  :",
          {round_h_up(params.k/params.n, sig_digits)})

    # Instantiation of parameters where alpha is set as well.
    lower_p, upper_p = exact(params)

    print("\n===== Exact interval of p with",
          params.confi_perc, "[%] two-sided (or",
          100 - (100 - params.confi_perc)/2,
          "[%] one-sided) confidence  =====")
    print(f"Upper : {round_h_up(upper_p, sig_digits)}")
    print(f"Lower : {round_h_up(lower_p, sig_digits)}")
    print(f"Width : {round_h_up(upper_p - lower_p, sig_digits)}")
    '''
    # NOTE: if the numbers are too small, use below.
    print(f"Upper : {result_upper[0]:.6g}")
    print(f"Lower : {result_lower[0]:.6g}")
    print(f"Width : {result_upper[0] - result_lower[0]:.6g}")
    # NOTE: be careful using round(), which rounds 0.*5 to the nearest even,
        say, 1.2345 to 1.234 (not 1.235).
    '''
    print("\n===== Verification =====")
    vmes = ("If the following does not hold, "
            "the results might not be reliable.\n")
    print(vmes)
    ret = verify_interval_of_p(
        params,
        lower_p, upper_p,
        sig_digits, verbose=2)
    if ret != 0:
        print("\nNG: obtained interval is not reliable, try other parameters.")
    else:
        print("\nOK: obtained interval is reliable.")


# ===== Graphs =====
class GraProps:
    """Parameter pool class for graph properties.

    Args:
        k_start: Start range of k.
        k_end: End range of k.
        k_step: Step of k in range(k_start, k_end, k_step).
        log_n_end: Max of x-axis given by (3 k_end * 10^log_n_end).
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
            leg_pos='auto'
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
        self.set_n(log_n_end)  # may exit

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
        # n <= 3 * k_end * 10^log_n_end <= 10^7 where 3 is floor(sqrt(10)).
        RELIABLE_LOG_N_MAX = 7
        if (10**RELIABLE_LOG_N_MAX < 3 * self.k_end * 10**log_n_end):
            log_n_end_org = log_n_end
            log_n_end = math.floor(
                np.log10((10**RELIABLE_LOG_N_MAX)/(3 * self.k_end)))
            logger.warning(
                f"log_n_end={log_n_end_org} was so large and changed to"
                f" {log_n_end}.")
            '''
            logger.error(
                f"log_n_end={log_n_end} is too large!!"
                " '3 k_end * 10^log_n_end <= "
                f"10^{RELIABLE_LOG_N_MAX}' must hold, "
                f"where k_end={self.k_end} "
                "and 3 k_end * 10^log_n_end="
                f"{3 * self.k_end * 10**log_n_end}!!")
            sys.exit(1)
            '''
        self.log_n_end = log_n_end


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

        >>> interval_graph(GraProps(
            k_start=1,  # >= 0
            k_end=1,    # > k_start
            k_step=1,   # > 0
            log_n_end=6,  # max(n) = 3*k_end*10**log_n_end
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
            ],
            # savefig = True,
            # fig_file_name='intervals.jpg'
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
    trans_alpha = 0.6  # Transparency

    # Set colorlist
    # With other lines
    line_types_per_k = 1  # Exact interval
    if k_diff == 0:
        if k_start == 0:
            if 'with_rule_of_la' in line_list:
                line_types_per_k += 1
            if 'with_wilson_cc' in line_list:
                line_types_per_k += 1
        else:
            if 'with_normal' in line_list:
                line_types_per_k += 1
            if 'with_wilson' in line_list:
                line_types_per_k += 1
            if 'with_wilson_cc' in line_list:
                line_types_per_k += 1

    # Max num of colors
    col_max_per_k = line_types_per_k * len(gra_props.confi_perc_list)
    col_max = (k_end - k_start + 1) * col_max_per_k
    cm = plt.get_cmap('brg')
    # yellow lines are haed to see.
    # cm = plt.get_cmap('gist_rainbow')
    # cm = plt.get_cmap('rainbow')
    cNorm = colors.Normalize(vmin=0, vmax=col_max - 1)
    scalarMap = mplcm.ScalarMappable(norm=cNorm, cmap=cm)
    colorlist = [scalarMap.to_rgba(i) for i in range(col_max)]

    n_list_base = (
        sorted([10**i for i in range(0, log_n_end + 1, 1)]
               + [3 * 10**i for i in range(0, log_n_end + 1, 1)]))
    # print(n_list_base)

    params = Params()
    leg = 0  # Num of legend
    if (k_start == 0):
        k = 0  # lower_p = 0
        params.set_k(k)
        n_list = n_list_base
        col_offset = 0
        for confi_perc in gra_props.confi_perc_list:
            params.set_confi_perc(confi_perc)
            if 'with_exact' in line_list:
                k0_upper = [exact(params)[1] for params.n in n_list]
                plt.plot(
                    n_list, k0_upper, color=colorlist[k + col_offset],
                    alpha=trans_alpha, marker=3,
                    # linestyle='dashdot',
                    linestyle='solid',
                    label=f"$k$={k} {confi_perc}% exact upper")
                leg += 1
                if k_diff == 0:
                    col_offset += 1
            # Approximated upper bound
            # col_offset = 0
            if 'with_rule_of_la' in line_list:
                k0_rule_of_ln_a = [rule_of_ln_alpha(
                    params)[1] for params.n in n_list]
                plt.plot(n_list,
                         k0_rule_of_ln_a,
                         color=colorlist[k + col_offset],
                         alpha=trans_alpha,
                         marker=3,
                         # linestyle='dotted',
                         linestyle='dashdot',
                         label=f"$k$={k} {confi_perc}% rule of -ln(a)")
                leg += 1
                if k_diff == 0:
                    col_offset += 1
            if 'with_wilson_cc' in line_list:
                # Wilson with continuity correction
                tmp = np.array(
                    [wilson_score_cc(params) for params.n in n_list])
                wilson_cc_lower = tmp[:, 0]
                wilson_cc_upper = tmp[:, 1]
                plt.plot(n_list,
                         wilson_cc_upper,
                         color=colorlist[
                             (k - k_start) * col_max_per_k + col_offset],
                         alpha=trans_alpha,
                         marker=3,
                         # linestyle='dotted',
                         linestyle='dashdot',
                         label=f"$k$={k} {confi_perc}% Wilson cc upper")
                leg += 1
                plt.plot(n_list,
                         wilson_cc_lower,
                         color=colorlist[
                             (k - k_start) * col_max_per_k + col_offset],
                         alpha=trans_alpha,
                         marker=2,
                         # linestyle='dotted',
                         linestyle='dashdot',
                         label=f"$k$={k} {confi_perc}% Wilson cc lower")
                leg += 1
                if k_diff == 0:
                    col_offset += 1

    for k in range(max(1, k_start), k_end + 1, k_step):
        n_list = np.array(n_list_base) * k
        # print(n_list)
        if 'with_line_kn' in line_list:
            # Line of k/n
            k_expect = [k / n for n in n_list]
            plt.plot(n_list,
                     k_expect,
                     color=colorlist[(k - k_start) * col_max_per_k],
                     alpha=trans_alpha,
                     linestyle='dashed',
                     label=f"$k$={k} k/n={k}/n")
            leg += 1

        # Exact lower and upper
        params.set_k(k)
        col_offset = 0
        for confi_perc in gra_props.confi_perc_list:
            params.set_confi_perc(confi_perc)
            if 'with_exact' in line_list:
                tmp = np.array([exact(params) for params.n in n_list])
                k_lower = tmp[:, 0]
                k_upper = tmp[:, 1]
                plt.plot(n_list,
                         k_upper,
                         color=colorlist[
                             (k - k_start) * col_max_per_k + col_offset],
                         alpha=trans_alpha,
                         marker=3,
                         # linestyle='dashdot',
                         linestyle='solid',
                         label=f"$k$={k} {confi_perc}% exact upper")
                leg += 1
                plt.plot(n_list,
                         k_lower,
                         color=colorlist[
                             (k - k_start) * col_max_per_k + col_offset],
                         alpha=trans_alpha,
                         marker=2,
                         # linestyle='dashdot',
                         linestyle='solid',
                         label=f"$k$={k} {confi_perc}% exact lower")
                leg += 1
                if k_diff == 0:
                    col_offset += 1

            # Approximated Intervals
            # col_offset = 0
            if 'with_normal' in line_list:
                # Normal
                tmp = np.array([normal_approx(params) for params.n in n_list])
                norm_lower = tmp[:, 0]
                norm_upper = tmp[:, 1]
                plt.plot(
                    n_list, norm_upper, color=colorlist[
                        (k - k_start) * col_max_per_k + col_offset],
                    alpha=trans_alpha, marker=3,
                    # linestyle='dotted',
                    linestyle='dashdot',
                    label=f"$k$={k} {confi_perc}% normal upper")
                leg += 1
                plt.plot(
                    n_list, norm_lower, color=colorlist[
                        (k - k_start) * col_max_per_k + col_offset],
                    alpha=trans_alpha, marker=2,
                    # linestyle='dotted',
                    linestyle='dashdot',
                    label=f"$k$={k} {confi_perc}% normal lower")
                leg += 1
                if k_diff == 0:
                    col_offset += 1
            if 'with_wilson' in line_list:
                # Wilson
                tmp = np.array([wilson_score(params) for params.n in n_list])
                wilson_lower = tmp[:, 0]
                wilson_upper = tmp[:, 1]
                plt.plot(n_list,
                         wilson_upper,
                         color=colorlist[
                             (k - k_start) * col_max_per_k + col_offset],
                         alpha=trans_alpha,
                         marker=3,
                         # linestyle='dotted',
                         linestyle='dashdot',
                         label=f"$k$={k} {confi_perc}% Wilson upper")
                leg += 1
                plt.plot(n_list,
                         wilson_lower,
                         color=colorlist[
                             (k - k_start) * col_max_per_k + col_offset],
                         alpha=trans_alpha,
                         marker=2,
                         # linestyle='dotted',
                         linestyle='dashdot',
                         label=f"$k$={k} {confi_perc}% Wilson lower")
                leg += 1
                if k_diff == 0:
                    col_offset += 1
            if 'with_wilson_cc' in line_list:
                # Wilson with continuity correction
                tmp = np.array(
                    [wilson_score_cc(params) for params.n in n_list])
                wilson_cc_lower = tmp[:, 0]
                wilson_cc_upper = tmp[:, 1]
                plt.plot(n_list,
                         wilson_cc_upper,
                         color=colorlist[
                             (k - k_start) * col_max_per_k + col_offset],
                         alpha=trans_alpha,
                         marker=3,
                         # linestyle='dotted',
                         linestyle='dashdot',
                         label=f"$k$={k} {confi_perc}% Wilson cc upper")
                leg += 1
                plt.plot(n_list,
                         wilson_cc_lower,
                         color=colorlist[
                             (k - k_start) * col_max_per_k + col_offset],
                         alpha=trans_alpha,
                         marker=2,
                         # linestyle='dotted',
                         linestyle='dashdot',
                         label=f"$k$={k} {confi_perc}% Wilson cc lower")
                leg += 1
                if k_diff == 0:
                    col_offset += 1

    # Show
    # plt.figure(figsize=(400, 300), dpi=300)
    # plt.rcParams["figure.figsize"]=25,20
    plt.title(title)
    plt.xlabel(x_label, fontsize=12)
    plt.ylabel(y_label, fontsize=12)

    leg_fontsize = 8
    # print(leg_pos, (leg_pos == 'auto'))
    if (leg <= 5 and (leg_pos == 'auto')) or (leg_pos == 'upper_right'):
        plt.legend(
            bbox_to_anchor=(1, 1),
            loc='upper right',
            borderaxespad=1,
            fontsize=leg_fontsize)
    elif (leg <= 6 and (leg_pos == 'auto')) or (leg_pos == 'upper_right_nm'):
        plt.legend(
            bbox_to_anchor=(1, 1),
            loc='upper right',
            borderaxespad=0,
            fontsize=leg_fontsize)
    elif (leg_pos == 'auto') or (leg_pos == 'out_right'):
        plt.legend(
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
        # plt.savefig(gra_props.fig_file_name, dpi = 300)
        plt.savefig(gra_props.fig_file_name)


def compare_dist(params):
    """Compare major probability distributions in a graph.

    Args:
        params (Params): Instance including k, n and zah (Z_{alpha/2})
            where zah is set by giving either 'alpha' or 'confi_perc'
            to Params().

    Examples:
        >>> compare_dist(Params(k=2, n=20))
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
    plt.plot(x, st.binom.pmf(x, n, p), label="Binomial")

    # Poisson
    plt.plot(x, st.poisson.pmf(x, k), label="Poisson")

    # Normal distribution
    variance = n * p * (1 - p)
    sigma = np.sqrt(variance)
    x = range(-3, x_max + 1, 1)
    plt.plot(x, st.norm.pdf(x, k, sigma), label="Normal")

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
    print(
        "params.confi_perc_to_alpha(alpha_to_confi_perc_wo_check"
        "({alpha})) = {tmp}")
    tmp = params.confi_perc_to_alpha(
        alpha_to_confi_perc_wo_check(alpha))
    print(
        "params.confi_perc_to_alpha(alpha_to_confi_perc_wo_check"
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
    print(
        "params.alpha_to_confi_perc(confi_perc_to_alpha_wo_check"
        f"({confi_perc})) = {tmp}")
    tmp = params.alpha_to_confi_perc(
        confi_perc_to_alpha_wo_check(confi_perc))
    print(
        "params.alpha_to_confi_perc(confi_perc_to_alpha_wo_check"
        "({confi_perc})) = {tmp}")
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

    Check if integral of outside of the interval is alpha for given k, n, etc
    and return 0 if it holds.

    Returns:
        int: num_of_wrongs, the number of unreliable intervals.

    Examples:
        See ebcic.ipynb for the parameters where reliable intervas
        might not be able to be obtained.

        >>> test_of_intervals(
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
            n_end=100,
            n_step=1,
            confi_perc=99.0,  # Confidence Percentage (90 <= confi_perc < 100)
            sig_digits=-5)
    if num_of_wrongs > 0:
        ret = num_of_wrongs

    if ret == 0:
        print("\n== OK : All the tests are succeeded. ==\n")
    else:
        print("\n== NG : Some tests are failed!! ==\n")


if __name__ == "__main__":
    '''
    Examples:
    '''
