# How to use [EBCIC](https://kazkobara.github.io/ebcic/) in MATLAB

## The following was verified at

```console
MATLAB R2021a (64-bit (win64))
Python 3.8.10 (Windows installer (64-bit))
ebcic 0.0.1
Widows 10 (64-bit 20H2)
```

## Installation of Python to be used from MATLAB

1. Download and install Python Windows installer (`XX`-bit) from [the official page](https://www.python.org/):
   - `XX` shall be the same as that of MATLAB, which can be checked by any of the followings:
     - "[Help]->[About MATLAB]" on a MATLAB desktop menu.
     - Type `computer` or `mexext` on a MATLAB command prompt.
   - Python (3.8.10 at least) from Microsoft Store could not be used from MATLAB presumably due to the [Known Issues](https://docs.python.org/3.8/using/windows.html#known-issues).

## Installation of `ebcic` package to the installed Python

1. Check the absolute path to the install python.exe on a windows terminal, such as PowerShell and Command Prompt:  

    ```console
    where python.exe
    ```

1. Install `ebcic` package with:

    ```console
    <the_full_path_to_python>\python.exe -m pip install ebcic
    ```

    where `<the_full_path_to_python>` is like `C:\Users\<your_account_name>\AppData\Local\Programs\Python\Python38` and `<your_account_name>` is your login name to Windows.

## Setting in MATLAB

1. Set the python version to use on MATLAB according to [MathWorks' page](https://jp.mathworks.com/help/matlab/ref/pyenv.html?lang=en).

    ```matlab
    pe=pyenv("Version", '3.8')
    ```

    - If Python other than 'Windows installer' is installed, tell MATLAB the full path (absolute path) to the installed python.exe with:

        ```matlab
        pyenv("Version", "<the_full_path_to_python>\python.exe")
        ```

    - Even with this `pyenv()` command, Python (3.8.10 at least) from Microsoft Store could not be used from MATLAB.

1. If the following returns blank or `ans = ""`, the installed Python is not usable from MATLAB, and try another Python installation.

    ```matlab
    pe.Version
    pe.Executable
    ```

## How to use [EBCIC](https://kazkobara.github.io/ebcic/) from MATLAB

1. Open [ebcic_in_matlab.m](https://kazkobara.github.io/ebcic/matlab/ebcic_in_matlab.m) as a 'live script' with MATLAB Live Editor.
1. Edit and run the sections you want to execute.
1. (To commit, save it as a MATLAB code (\*.m) file, and then commit \*.m file to git since Live Code File Format (*.mlx) is not git friendly. :-))

## Examples

### Obtain intervals as a matlab cell and variables

```matlab
confidence_interval=cell(py.ebcic.exact(py.ebcic.Params(pyargs( ...
    'k',1,                  ...% # of errors
    'n',100,                ...% # of trials
    'confi_perc',confi_perc ...% Confidence percentage:
      ...%   for two-sided of 0<k<n where 0 < confi_perc < 100, or
      ...%   for one-sided of k=0 or k=n.
      ...% NOTE
      ...%   For one-sided of 0<k<n,
      ...%     set confi_perc=(2 * confi_perc_for_one_sided - 100)
      ...%     where 50 < confi_perc_for_one_sided < 100.
    ))));

% As variables.
lower_interval=confidence_interval{1,1}
upper_interval=confidence_interval{1,2}
```

Result:

```text
lower_interval = 5.0124e-05
upper_interval = 0.0720
```

### Depict graphs

To depict graphs in [this page](https://kazkobara.github.io/ebcic/), edit the parameters below. The following example corresponds with "Exact intervals and the line of k/n for k=1" where members of the 'confi_perc_list' are reduced (due to the following bug on 'the legend outside of the graph area'.).

```matlab
py.ebcic.interval_graph(py.ebcic.GraProps(pyargs( ...
    ...% Set the range of k to depict with k_*
    'k_start',int32(1), ...
    'k_end',  int32(1), ...
    ...% Edit the list of confidence percentages to depict, [confi_perc, ...],
    ...%   for two-sided of 0<k<n where 0 < confi_perc < 100, or
    ...%   for one-sided of k=0 or k=n.
    ...% NOTE For one-sided of 0<k<n, set 
    ...%   confi_perc=(2 * confi_perc_for_one_sided - 100)
    ...%   where 50 < confi_perc_for_one_sided < 100
    ...%   (though both lower and upper intervals are shown).
    'confi_perc_list',[95.0, 99.0], ...
    ...% Add or remove lines to depict.
    'line_list',[ ...
        'with_line_kn',        ...% Line of k/n
        'with_exact',          ...% Exact interval
        ...% 'with_wilson_cc', ...% Wilson score interval with continuity correction
        ...% 'with_wilson',    ...% Wilson score interval;         not available for k=0
        ...% 'with_normal',    ...% Normal approximation interval; not available for k=0
        ...% 'with_rule_of_la',...% Rule of -ln(alpha);           available only for k=0
        ] ...
    )))
```

### Troubleshooting

- If errors like "Can't find a usable xxx in the following directories" are displayed, fix the path so that missing xxx can be found in the path, e.g. for Python 3.8:
  - 'init.tcl' is in 'Python38\tcl\tcl8.6\', but expected to be in 'Python38\lib\tcl8.6\'
  - 'tk.tcl' is in 'Python38\tcl\tk8.6', but expected to be in 'Python38\lib\tk8.6\' or 'Python38\lib\tcl8.6\tk8.6'.

- To show another graph, the popped-up window with the current graph has to be closed.

### Known bug

- Legend outside of the graph area is not displayed in MATLAB. So it seems better to use Python directly to depict graphs.

---

[EBCIC README](https://kazkobara.github.io/ebcic/)
