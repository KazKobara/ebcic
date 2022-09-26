# EBCIC: Exact Binomial Confidence Interval Calculator

[![Downloads](https://pepy.tech/badge/ebcic)](https://pepy.tech/project/ebcic)
[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/KazKobara/ebcic.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/KazKobara/ebcic/context:python)
![GitHub tag (latest by date)](https://img.shields.io/github/v/tag/KazKobara/ebcic)
[![Total alerts](https://img.shields.io/lgtm/alerts/g/KazKobara/ebcic.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/KazKobara/ebcic/alerts/)
![GitHub issues](https://img.shields.io/github/issues/kazkobara/ebcic)
![GitHub](https://img.shields.io/github/license/kazkobara/ebcic)

日本語([github.io](https://kazkobara.github.io/ebcic/README-jp.html), [github](https://github.com/KazKobara/ebcic/blob/master/README-jp.md))
<img src="https://raw.githubusercontent.com/lipis/flag-icons/main/flags/4x3/jp.svg" width="20" alt="Japanese" title="Japanese"/>

These programs are mainly for researchers, developers, and designers who calculate Binomial Confidence Intervals.

<!--for given parameters:

- `n`: the number of Bernoulli or Binomial trials.
- `k`: the number of target events happened.
- Confidence parameter given
  - `confi_perc`: confidence percentage:
    - for two-sided of `0<k<n` where `0 < confi_perc < 100`, or for one-sided of `k=0` or `k=n`.
    - for one-sided of `0<k<n`, set `confi_perc = (2 * confi_perc_for_one_sided - 100)` where `50 < confi_perc_for_one_sided < 100`.
-->

[EBCIC](https://kazkobara.github.io/ebcic/) calculates binomial intervals exactly, i.e. by implementing Clopper-Pearson interval [CP34] without simplifying mathematical equations that may deteriorate intervals for certain combinations of parameters. EBCIC can also shows graphs for comparing exact intervals with approximated ones.

## How to use

### Jupyter notebook

1. Open [ebcic.ipynb](https://github.com/KazKobara/ebcic/blob/master/ebcic.ipynb) with Jupyter-notebook-compatible development environment such as Jupyter Notebook, JupyterLab, or Visual Studio Code.
2. Run the following initial cells:

    ~~~python
    # Run this cell, if `ebcic` package has not been installed yet:
    %pip install ebcic
    ~~~

    ~~~python
    import ebcic
    from ebcic import *
    ~~~

3. Run the cells you want to execute.

### Command line

1. Installation

    - When using [PyPI ebcic package](https://pypi.org/project/ebcic/):

        ~~~console
        pip install ebcic
        ~~~

    - When using github `ebcic` repo:

        ~~~console
        git clone https://github.com/KazKobara/ebcic.git
        cd ebcic
        ~~~

2. Command-line help

    - Check the version and options:

        ~~~console
        python -m ebcic -h
        ~~~

3. Cf. the [examples](#command-line-examples) below.

### MATLAB (with Python and `ebcic` package)

1. Install Python for MATLAB and `ebcic` package according to [this page](https://kazkobara.github.io/ebcic/matlab/ebcic_in_matlab.html).
2. Open a sample MATLAB code file [ebcic_in_matlab.m](https://kazkobara.github.io/ebcic/matlab/ebcic_in_matlab.m) as a 'live script' as shown [this page](https://jp.mathworks.com/help/matlab/matlab_prog/create-live-scripts.html?lang=en).
3. Edit and run the sections you want to execute.

> NOTE: If you manage the edited file with git, save it as a MATLAB code file (*\.m) file to commit (or commit the live code file (\*.mlx) to a git LFS (Large File Storage)) since live code files (\*.mlx) are not git friendly. If necessary, save it as a \*.html file as well to check its look.

## Command Line Examples

To print exact intervals as text.

### One-sided upper 95% confidence interval for no error among 100 trials

~~~console
python -m ebcic -k 0 -n 100 -c 95 -u
~~~

> - For `k=0` or `k=n`, give `-c` option **one-sided** confidence percentage.
> - v0.0.4 or newer returns the same value as the above result by setting `--rej-perc-lower` (`-r`) option the percentage of the lower rejection area in assuming population as follows:

~~~console
python -m ebcic -k 0 -n 100 --rej-perc-lower 5 -u
~~~

### Two-sided 95% confidence interval for one error among 100 trials

~~~console
python -m ebcic -k 1 -n 100 -c 95 -lu
~~~

> - For `0<k<n`, give `-c` option **two-sided** confidence percentage.
> - v0.0.4 or newer returns the same value as the above result by setting both `--rej-perc-lower` (`-r`) and `--rej-perc-upper` (`-s`) options equally divided percentages of assuming population as follows:

~~~console
python -m ebcic -k 1 -n 100 -r 2.5 -s 2.5 -lu
~~~

### One-sided upper 95% confidence interval for one error among 100 trials

~~~console
python -m ebcic -k 1 -n 100 -r 5 -u
~~~

> - For v0.0.4 and newer, set `--rej-perc-lower` (`-r`) option the percentage of the lower rejection area in assuming population.
> - For v0.0.3 and older and `0<k<n`, give `-c` option `2*s-100` as follows where `s` is the one-sided confidence percentage (in this case `s=95` and `2*s-100=2*95-100=90`).

~~~console
python -m ebcic -k 1 -n 100 -c 90 -u
~~~

> Giving `-c 90` is the same as giving `--alpha 0.1` (or `-a 0.1`).

~~~console
python -m ebcic -k 1 -n 100 --alpha 0.1 -u
~~~

### One-sided lower 95% confidence interval for 99 errors among 100 trials

~~~console
python -m ebcic -k 99 -n 100 -s 5 -l
~~~

> - For v0.0.4 and newer, set `--rej-perc-upper` (`-s`) option the percentage of the upper rejection area in assuming population.
> - For v0.0.3 and older and `0<k<n`, the equivalent value is obtained in the same way as the previous example using `-c` or `-a` option as follows:

~~~console
python -m ebcic -k 99 -n 100 -c 90 -l
python -m ebcic -k 99 -n 100 -a 0.1 -l
~~~

## Python Interpreter or Jupyter Cell Examples

### Print exact interval as text

Edit the following parameters, `k`, `n`, and `confi_perc` (or `rej_perc_lower` and `rej_perc_upper`), and then run the cell.

~~~python
print_interval(Params(
    k=1,             # Number of errors
    n=501255,        # Number of trials
    confi_perc=99.0  # Confidence percentage
    ))
~~~

where `confi_perc` is set as follows:

- For `k=0` or `k=n`:
  - give the percentage of **one-sided** confidence interval.
- For `0<k<n`:
  - give the percentage of **two-sided** confidence interval.
- For `0<k<n` and to get one-sided confidence interval:
  - give `-c` option `2*s-100` where `s` is the one-sided confidence percentage.

Result:

~~~python
===== Exact interval of p with 99.0 [%] two-sided (or 99.5 [%] one-sided) confidence  =====
Upper : 1.482295806e-05
Lower : 9.99998e-09
Width : 1.481295808e-05
~~~

For v0.0.4 and newer, instead of `confi_perc` or `alpha`, Params() can set the confidence with
either or both of `rej_perc_lower` and `rej_perc_upper` in percentage of `0 <= x < 50` (or either or both of `rej_lower` and `rej_upper` in ratio of `0 <= x < 0.5`. Params()'s class functions are also available.

~~~python
Params(
    k=1,                # Number of errors
    n=501255,           # Number of trials
    # Rejection area in percentage
    rej_perc_lower=0.5  # Lower rejection area (to get upper interval)
    rej_perc_upper=0.5  # Upper rejection area (to get lower interval)
    ).print_interval()
~~~

> Note that it uses the **lower** rejection area to get the **upper** confidence interval and vice versa.

Result:

~~~text
===== Exact interval of p with rejection area of lower 0.5 [%] and upper 0.5 [%] =====
Upper :  1.482295806e-05
Lower :  9.99998e-09
Width :  1.481295808e-05
~~~

### Depict graphs

#### Exact intervals and the line of k/n for k=1

This program can show not only the typical 95% and 99% confidence lines but also any confidence percentage lines.

Python Interpreter or Jupyter cell to run:

~~~python
interval_graph(GraProps(
    # Set the range of k with k_*
    k_start=1,  # >= 0
    k_end=1,    # >= k_start
    k_step=1,   # >= 1
    # Edit the list of confidence percentages to depict, [confi_perc, ...],
    #   for two-sided of 0<k<n where 0 < confi_perc < 100, or
    #   for one-sided of k=0 or k=n.
    # NOTE For one-sided of 0<k<n, set 
    #   confi_perc=(2 * confi_perc_for_one_sided - 100)
    #   where 50 < confi_perc_for_one_sided < 100
    #   (though both lower and upper intervals are shown).
    confi_perc_list=[90, 95, 99, 99.9, 99.99],
    # Lines to depict
    line_list=[
        'with_exact',
        'with_line_kn',  # Line of k/n
    ],
    # savefig=True,  # uncomment on Python Interpreter 
    # fig_file_name='intervals.png',
    ))
~~~

Result:

If figures or links are not shown appropriately, visit [github.io page](https://kazkobara.github.io/ebcic/) or [github page](https://github.com/KazKobara/ebcic).

![Exact intervals and the line of k/n for k=1](https://github.com/KazKobara/ebcic/raw/master/figs/confidence_percentage.png)

#### Exact intervals for k=0 to 5

Python Interpreter or Jupyter cell to run:

~~~python
interval_graph(GraProps(
    k_start=0,  # >= 0
    k_end=5,    # >= k_start
    line_list=['with_exact'],
    # savefig=True,  # uncomment on Python Interpreter 
    # fig_file_name='intervals.png',
    ))
~~~

Result:

![Exact intervals for k=0 to 5](https://github.com/KazKobara/ebcic/raw/master/figs/num_of_errors.png)

#### Comparison of exact and approximated intervals for k=0

Python Interpreter or Jupyter cell to run:

~~~python
interval_graph(GraProps(
    k_start=0,    # >= 0
    k_end=0,      # >= k_start
    log_n_end=3,  # max(n) = k_end*10**log_n_end
    line_list=[
        'with_exact',
        'with_rule_of_la',  # rule of -ln(alpha)
                            # available only for k=0 and k=n
        #'with_normal',     # not available for k=0 and k=n
        'with_wilson',
        'with_wilson_cc',
        'with_beta_approx',
    ],
    # savefig=True,  # uncomment on Python Interpreter 
    # fig_file_name='intervals.png',
    ))
~~~

where interval names to be added in the `line_list` and their conditions are as follows:

Interval name (after 'with_')    | Explanation | Condition
:---------   |:----  |:--
exact | Implementation of Clopper-Pearson interval [CP34] without approximation. |
rule_of_la | '`Rule of -ln(a)`' or '`Rule of -log_e(alpha)`'; Generalization of the '`Rule of three`' [Lou81,HL83,JL97,Way00,ISO/IEC19795-1] that is for `k=0` and `alpha=0.05` (95% confidence percentage), to other confidence percentages than 95% and `k=n`. | `k=0` or `k=n`
wilson | `Wilson score interval` [Wil27].
wilson_cc | `Wilson score interval with continuity correction` [New98].
beta_approx | Approximated interval using beta distribution.
normal | `Normal approximation interval` or `Wald confidence interval`. | `0<k<n`

Result:

As you can see from the following figure, '`rule of -ln(a)`' for large `n` and '`beta_approx`' are good approximations for `k=0`.

> For `k=0`, interval_graph() of EBCIC v0.0.3 and newer, display only one-sided upper intervals since their lower intervals must be `0` (though some approximations, such as '`Wilson cc`', output wrong values than `0`).

![Comparison of exact and approximated intervals for k=0](https://github.com/KazKobara/ebcic/raw/master/figs/comparison_k0.png)

#### Comparison of exact and approximated intervals for `k=1`

Python Interpreter or Jupyter cell to run:

~~~python
interval_graph(GraProps(
    k_start=1,  # >= 0
    k_end=1,    # >= k_start
    line_list=[
        'with_line_kn'
        # 'with_rule_of_la',  # available only for k=0
        'with_exact',
        'with_normal',
        'with_wilson',
        'with_wilson_cc',
        'with_beta_approx',
    ],
    # savefig=True,  # uncomment on Python Interpreter 
    # fig_file_name='intervals.png',
    ))
~~~

Result:

As you can see from the following figures and as warned in many papers such as [BLC01], normal-approximation intervals are not good approximations for small `k`.

The upper intervals of the other approximations look tight.
The approximation using beta distribution looks tight where the confidence interval for `k=n=1` is one-sided.

![Comparison of exact and approximated intervals for k=1](https://github.com/KazKobara/ebcic/raw/master/figs/comparison_k1.png)

#### Comparison of exact and approximated intervals for `k=10`

Python Interpreter or Jupyter cell to run:

~~~python
interval_graph(GraProps(
    k_start=10,   # >= 0
    k_end=10,     # >= k_start
    log_n_end=2,  # max(n) = k_end*10**log_n_end
    line_list=[
        'with_exact',
        'with_normal',
        'with_wilson',
        'with_wilson_cc',
        'with_beta_approx',
    ],
    # savefig=True,  # uncomment on Python Interpreter 
    # fig_file_name='intervals.png',
    ))
~~~

Result:

For `k=10`, '`normal`' still does not provide a good approximation.

![Comparison of exact and approximated intervals for k=10](https://github.com/KazKobara/ebcic/raw/master/figs/comparison_k10.png)

#### Comparison of exact and approximated intervals for `k=100`

Python Interpreter or Jupyter cell to run:

~~~python
interval_graph(GraProps(
    k_start=100,  # >= 0
    k_end=100,    # >= k_start
    log_n_end=2,  # max(n) = k_end*10**log_n_end
    line_list=[
        'with_exact',
        'with_normal',
        'with_wilson',
        'with_wilson_cc',
        'with_beta_approx',
    ],
    # savefig=True,  # uncomment on Python Interpreter 
    # fig_file_name='intervals.png',
    ))
~~~

Result:

At least for `k=100` and confidence percentage, `confi_perc=99.0`, all these approximations look tight.

![Comparison of exact and approximated intervals for k=100](https://github.com/KazKobara/ebcic/raw/master/figs/comparison_k100.png)

### [API Manual](https://github.com/KazKobara/ebcic/tree/master/docs/_build)

1. Download

    ~~~console
    git clone https://github.com/KazKobara/ebcic.git
    ~~~

2. Open the following file with your browser (after replacing `<path to the downloaded ebcic>` appropriately):

    ~~~text
    file://<path to the downloaded ebcic>/docs/_build/index.html
    ~~~

    For WSL Ubuntu-20.04, replace `<username>` and `<path to the downloaded ebcic>` appropriately:

    ~~~text
    file://wsl%24/Ubuntu-20.04/home/<username>/<path to the downloaded ebcic>/docs/_build/index.html
    ~~~

## Bibliography

[CP34] Clopper, C. and Pearson, E.S. "The use of confidence or fiducial limits illustrated in the case of the binomial," Biometrika. 26 (4): pp.404-413, 1934

[Lou81] Louis, T.A. "Confidence intervals for a binomial parameter after observing no successes," The American Statistician, 35(3), p.154, 1981

[HL83] Hanley, J.A. and Lippman-Hand, A. "If nothing goes wrong, is everything all right? Interpreting zero numerators," Journal of the American Medical Association, 249(13), pp.1743-1745, 1983

[JL97] Jovanovic, B.D. and Levy, P.S. "A look at the rule of three," The American Statistician, 51(2), pp.137-139, 1997

[Way00] Wayman, J.L. "Technical testing and evaluation of biometric identification devices," Biometrics: Personal identification in networked society, edited by A.K. Jain, et al., Kluwer, pp.345-368, 2000

[ISO/IEC19795-1] ISO/IEC 19795-1, "Information technology-Biometric performance testing and reporting-Part 1: Principles and framework" <!-- 2006 2021 -->

[New98] Newcombe, R.G. "Two-sided confidence intervals for the single proportion: comparison of seven methods," Statistics in Medicine. 17 (8): pp.857-872, 1998

[Wil27] Wilson, E.B. "Probable inference, the law of succession, and statistical inference," Journal of the American Statistical Association. 22 (158): pp.209-212, 1927

[BLC01] Brown, L.D., Cai, T.T. and DasGupta, A. "Interval Estimation for a Binomial Proportion," Statistical Science. 16 (2): pp. 101-133, 2001

## [Changelog](./CHANGELOG.md)

## License

[MIT License](./LICENSE)

When you use or publish the confidence interval obtained with the software, please **refer to the software name, version, platform**, and so on, so that readers can verify the correctness and reproducibility of the interval with the input parameters.

An example of the reference is:

~~~text
The confidence interval is obtained by EBCIC X.X.X on Python 3."
~~~

where X.X.X is the version of EBCIC.

The initial software is based on results obtained from a project, JPNP16007, commissioned by the New Energy and Industrial Technology Development Organization (NEDO).

Copyright (c) 2020-2022 National Institute of Advanced Industrial Science and Technology (AIST)

---

- [https://github.com/KazKobara/](https://github.com/KazKobara/)
- [https://kazkobara.github.io/](https://kazkobara.github.io/)
