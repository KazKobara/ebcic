# EBCIC: Exact Binomial Confidence Interval Calculator

[![Downloads](https://pepy.tech/badge/ebcic)](https://pepy.tech/project/ebcic)
[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/KazKobara/ebcic.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/KazKobara/ebcic/context:python)
![GitHub tag (latest by date)](https://img.shields.io/github/v/tag/KazKobara/ebcic)
[![Total alerts](https://img.shields.io/lgtm/alerts/g/KazKobara/ebcic.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/KazKobara/ebcic/alerts/)
![GitHub issues](https://img.shields.io/github/issues/kazkobara/ebcic)
![GitHub](https://img.shields.io/github/license/kazkobara/ebcic)

[English <img src="https://raw.githubusercontent.com/lipis/flag-icons/main/flags/4x3/gb.svg" width="20" alt="English" title="English"/>](./README.md)

ベルヌーイ試行において母比率の信頼区間を厳密に計算したり、他の近似区間との差などを図示したりするための Python スクリプト

## 使い方

### Jupyter notebook で使う場合

1. [ebcic.ipynb](https://kazkobara.github.io/ebcic/ebcic.ipynb) を Jupyter, JupyterLab, Visual Studio Code などの Jupyter notebook を実行可能な環境で開きます。
2. まずは、以下のセルを実行:

    ~~~python
    %pip install ebcic
    ~~~

    ~~~python
    import ebcic
    from ebcic import *
    ~~~

3. 実行したいセルを実行

### コマンドラインで使う場合

1. インストール

    - [PyPI ebcic パッケージ](https://pypi.org/project/ebcic/) を使う場合:

        ~~~console
        pip install ebcic
        ~~~

    - github レポジトリのコードを使う場合:

        ~~~console
        git clone https://github.com/KazKobara/ebcic.git
        cd ebcic
        ~~~

2. コマンドラインヘルプ

    - 引数の使い方、バージョン情報などの表示

        ~~~console
        python -m ebcic -h
        ~~~

3. 以下の[コマンドラインでの実行例](#コマンドラインでの実行例)もご参照下さい。

## MATLABから呼び出す場合

1. [こちらのページ](https://kazkobara.github.io/ebcic/matlab/ebcic_in_matlab.html)に従って Python for MATLAB と `ebcic` package をインストール
2. [こちらのページ](https://jp.mathworks.com/help/matlab/matlab_prog/create-live-scripts.html)に従って
 [ebcic_in_matlab.m](https://kazkobara.github.io/ebcic/matlab/ebcic_in_matlab.m) をライブスクリプトとして開く
3. 実行したいセクションを編集して実行

> ※ ライブスクリプトファイル (\*.mlx) は git との相性がよくないため、MATLAB code file (*\.m) (と確認用の \*.html)をコミットするか、git LFS (Large File Storage)にライブスクリプトファイルを保存するようにした方がよいです。

## コマンドラインでの実行例

厳密な信頼区間を出力するためのコマンド例

### 試行100回中、誤り0回の誤り率の95%片側信頼区間の上限の表示

~~~console
python -m ebcic -k 0 -n 100 -c 95 -u
~~~

> `k=0` または `k=n` の場合は、**片側信**頼区間のパーセンテージを -c オプションで指定します。(前者の下限が0で、後者の上限が1であることは自明なため。)

### 試行100回中、誤り1回の誤り率の95%両側信頼区間の下限と上限の表示

~~~console
python -m ebcic -k 1 -n 100 -c 95 -lu
~~~

> `0<k<n` の場合は、**両側**信頼区間のパーセンテージを -c オプションで指定します。

### 上記の95%片側信頼区間の上限の表示

~~~console
python -m ebcic -k 1 -n 100 -c 90 -u
~~~

> - `0<k<n` の場合に片側信頼区間の値を得るには、片側信頼区間のパーセンテージを `s` として `2*s-100` を -c オプションに指定します。
> - 上記の場合は `2*95-100=90`
> - 理由: 両側信頼区間では `100-s` パーセントの棄却域を反対側にも設けるため`2*(100-s)` を `100` から引いて `100-2*(100-s)=2*s-100`。

## Jupyterセル/Pythonインタプリタ での実行例

### 厳密な信頼区間の数値での出力

以下の k, n, confi_perc を変更して実行

~~~python
print_interval(Params(
    k=1,             # 誤りの数
    n=501255,        # 試行数
    confi_perc=99.0  # 信頼区間のパーセンテージ [(1-信頼係数α)*100]
    ))
~~~

ここで信頼区間のパーセンテージの指定方法:

- `k=0` または `k=n` の場合:
  - **片側**信頼区間のパーセンテージを指定します。
- `0<k<n` の場合:
  - **両側**信頼区間のパーセンテージを入力します。
- `0<k<n` の場合に、片側信頼区間のパーセンテージを求めたい場合:
  - 片側信頼区間のパーセンテージを `s` として `2*s-100` を指定します。

実行結果:

~~~python
===== Exact interval of p with 99.0 [%] two-sided (or 99.5 [%] one-sided) confidence  =====
Upper : 1.482295806e-05
Lower : 9.99998e-09
Width : 1.481295808e-05
~~~

### グラフの描画

#### k=1の場合の標本誤り率 k/n と信頼区間のパーセンテージを変えた場合の厳密な信頼区間

よく使用される `95%` や `99%` 以外との比較も可能です。

~~~python
interval_graph(GraProps(
    # k=1の場合のみを表示するため各1を指定
    k_start=1,  # >= 0
    k_end=1,    # >= k_start
    # 描画する信頼区間のパーセンテージの集合を指定
    confi_perc_list=[90, 95, 99, 99.9, 99.99],
    # 描画する線の種類を指定
    line_list=[
        'with_exact',    # 厳密な信頼区間
        'with_line_kn',  # 標本誤り率 k/n
    ],
    # savefig=True,  # Pythonインタプリタで実行する場合にはコメントを外す
    # fig_file_name='intervals.png',  # 描画ファイル名を指定
    ))
~~~

描画結果:

(もし、図が表示されていない場合には、[こちら](https://kazkobara.github.io/ebcic/)をご参照下さい。)

![Exact intervals and the line of k/n for k=1](./figs/confidence_percentage.png)

#### kを0から5まで変更した場合の厳密な信頼区間

~~~python
interval_graph(GraProps(
    k_start=0,  # >= 0
    k_end=5,    # >= k_start
    line_list=['with_exact'],
    # savefig=True,  # Pythonインタプリタで実行する場合にはコメントを外す
    # fig_file_name='intervals.png',  # 描画ファイル名を指定
    ))
~~~

描画結果:

![Exact intervals for k=0 to 5](./figs/num_of_errors.png)

#### `k=0`の場合の厳密な信頼区間と近似的な信頼区間との比較

~~~python
interval_graph(GraProps(
    k_start=0,    # >= 0
    k_end=0,      # >= k_start
    log_n_end=3,  # max(n) = k_end*10**log_n_end
    line_list=[
        'with_exact',        # 厳密な信頼区間
        'with_rule_of_la',   # rule of -ln(alpha)
                             # k=0またはk=nの場合でのみ使用可能
        #'with_normal',      # 0<k<n の場合でのみ使用可能
        'with_wilson',       # Wilson
        'with_wilson_cc',    # Wilson cc
        'with_beta_approx',  # approximation using beta function
    ],
    # savefig=True,  # Pythonインタプリタで実行する場合にはコメントを外す
    # fig_file_name='intervals.png',  # 描画ファイル名を指定
    ))
~~~

ここで、line_list に追加可能な信頼区間名と条件は以下のとおりです。

信頼区間名 ('with_'の右側)    | 説明  | 条件
:---------   |:----  |:--
exact | Clopper-Pearson [CP34] の考え方を近似を行わずに計算した区間 |
rule_of_la | `k=0` の近似信頼区間である '`Rule of three`' [Lou81,HL83,JL97,Way00,ISO/IEC19795-1]を 95% 以外の信頼区間と `k=n` にも適用できるように一般化した近似区間 ('`Rule of -ln(a)`'または'`Rule of -log_e(alpha)`')  | `k=0` or `k=n`
wilson | `Wilson score interval` [Wil27] の近似区間
wilson_cc | `Wilson score interval with continuity correction` [New98] の近似区間
beta_approx | ベータ関数を使った近似区間
normal | 二項分布を正規分布へ近似して求めた区間(`Normal approximation interval` または `Wald confidence interval`) | `0<k<n`

描画結果:

`k=0` の場合は、'`beta_approx`' と `n`が大きな場合に '`Rule of -ln(a)`' が良い近似になっていることが分かります。

> EBCIC 0.0.3以降の interval_graph() では、 `k=0` の信頼区間は片側の上限のみを描画するようにしてあります。(`k=0`の場合の下限は`0`であることが自明なのですが、'`Wilson cc`'などの近似を用いる方法では `0` とは異なる値が出力されるため。)

![Comparison of exact and approximated intervals for k=0](./figs/comparison_k0.png)

#### `k=1`の場合の厳密な信頼区間と近似的な信頼区間との比較

~~~python
interval_graph(GraProps(
    k_start=1,  # >= 0
    k_end=1,    # >= k_start
    line_list=[
        'with_line_kn'
        # 'with_rule_of_la',  # k=0 の場合でのみ使用可能
        'with_exact',
        'with_normal',        # 0<k<n の場合でのみ使用可能
        'with_wilson',
        'with_wilson_cc',
        'with_beta_approx',
    ],
    # savefig=True,  # Pythonインタプリタで実行する場合にはコメントを外す
    # fig_file_name='intervals.png',  # 描画ファイル名を指定
    ))
~~~

描画結果:

- [BLC01]などの多くの文献で指摘されているとおり、'`normal`' (`Normal approximation interval` または `Wald confidence interval`)は `k` が小さい場合によい近似となっていないことが分かります。
- 上限については '`normal`' **以外**はよい近似になっていることが分かります。
- ベータ関数を使った近似は上下限共によい近似になっていることが分かります。(`k=n=1` の場合の信頼区間は片側になります。)

![Comparison of exact and approximated intervals for k=1](./figs/comparison_k1.png)

#### `k=10`の場合の厳密な信頼区間と近似的な信頼区間との比較

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
    # savefig=True,  # Pythonインタプリタで実行する場合にはコメントを外す
    # fig_file_name='intervals.png',  # 描画ファイル名を指定
    ))
~~~

描画結果:

`k=10`の場合でも、'`normal`'などはまだよい近似とはなっていないことが分かります。

![Comparison of exact and approximated intervals for k=10](./figs/comparison_k10.png)

#### `k=100`の場合の厳密な信頼区間と近似的な信頼区間との比較

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
    # savefig=True,  # Pythonインタプリタで実行する場合にはコメントを外す
    # fig_file_name='intervals.png',  # 描画ファイル名を指定
    ))
~~~

描画結果:

`k=100`で`confi_perc=99.0`の場合は、本図で比較した近似的な信頼区間はいずれも厳密な信頼区間のよい近似になっていることが分かります。

![Comparison of exact and approximated intervals for k=20](./figs/comparison_k100.png)

## [APIマニュアル](https://github.com/KazKobara/ebcic/tree/master/docs/_build)

1. ダウンロード

    ~~~console
    git clone https://github.com/KazKobara/ebcic.git
    ~~~

2. 以下の `<path to the downloaded ebcic>` を上記でダウンロードしたフォルダに変更しWebブラウザで開く:

    ~~~text
    file://<path to the downloaded ebcic>/docs/_build/index.html
    ~~~

    上記ダウンロードフォルダが WSL の Ubuntu-20.04 の下なら、以下の `<username>` と `<path to the downloaded ebcic>` を変更しWebブラウザで開く:

    ~~~text
    file://wsl%24/Ubuntu-20.04/home/<username>/<path to the downloaded ebcic>/docs/_build/index.html
    ~~~

## 参考文献

[CP34] Clopper, C. and Pearson, E.S. "The use of confidence or fiducial limits illustrated in the case of the binomial," Biometrika. 26 (4): pp.404-413, 1934

[Lou81] Louis, T.A. "Confidence intervals for a binomial parameter after observing no successes," The American Statistician, 35(3), p.154, 1981

[HL83] Hanley, J.A. and Lippman-Hand, A. "If nothing goes wrong, is everything all right? Interpreting zero numerators," Journal of the American Medical Association, 249(13), pp.1743-1745, 1983

[JL97] Jovanovic, B.D. and Levy, P.S. "A look at the rule of three," The American Statistician, 51(2), pp.137-139, 1997

[Way00] Wayman, J.L. "Technical testing and evaluation of biometric identification devices," Biometrics: Personal identification in networked society, edited by A.K. Jain, et al., Kluwer, pp.345-368, 2000

[ISO/IEC19795-1] ISO/IEC 19795-1, "Information technology - Biometric performance testing and reporting - Part 1: Principles and framework" <!-- 2006 2021 -->

[New98] Newcombe, R.G. "Two-sided confidence intervals for the single proportion: comparison of seven methods," Statistics in Medicine. 17 (8): pp.857-872, 1998

[Wil27] Wilson, E.B. "Probable inference, the law of succession, and statistical inference," Journal of the American Statistical Association. 22 (158): pp.209-212, 1927

[BLC01] Brown, L.D., Cai, T.T. and DasGupta, A. "Interval Estimation for a Binomial Proportion," Statistical Science. 16 (2): pp. 101-133, 2001

## [Changelog](./CHANGELOG.md)

## [MIT License](./LICENSE)

---

- [https://github.com/KazKobara/](https://github.com/KazKobara/)
- [https://kazkobara.github.io/](https://kazkobara.github.io/README-jp.html)