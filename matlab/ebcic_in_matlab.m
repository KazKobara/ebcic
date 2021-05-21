%% How to use <https://kazkobara.github.io/ebcic/ EBCIC> (Exact Binomial Confidence Interval Calculator) from MATLAB
% Open this file as a 'live script', but save as a MATLAB code file (*.m) file 
% (and a *.html file to view) since live code files (*.mlx) are not git friendly. 
% :-)
% 
% Commands in each section are executed with 'Ctrl+Enter' on Windows (and '⌘⏎' 
% on MacOS), cf. <https://www.mathworks.com/content/dam/mathworks/fact-sheet/live-editor-quick-start-guide.pdf 
% live-editor-quick-start-guide.pdf>.
%% 1. Check if your MATLAB is 64-bit or 32-bit.

computer
mexext
%% 2. Install Python and ebcic package to use in MATLAB according to <https://kazkobara.github.io/ebcic/README-MATLAB.html this page>.
%% 3. Choose Python version and see if it is usable from MATLAB.
% If "pe.Version" below is empty, the installed Python is not usable from MATLAB.

pe=pyenv("Version",'3.8');
% pe.Executable  % Uncomment to see the absolute path to python.
pe.Version
%% 4. Check if ebcic package is available.

py.help("ebcic")
py.ebcic.print_interval(py.ebcic.Params(pyargs( ...
    'k',1, ...
    'n',100, ...
    'confi_perc',99)))
%% Obtain exact interval

confidence_interval=cell(py.ebcic.exact(py.ebcic.Params(pyargs( ...
    'k',1,                  ...% # of errors
    'n',100,                ...% # of trials
    'confi_perc',confi_perc ...% Confidence percentage for two-sided where 0 < confi_perc < 100.
                            ...% For one-sided, set confi_perc=(2 * confi_perc_for_one_sided - 100)
                            ...% where 50 < confi_perc_for_one_sided < 100.
    ))));

% Show intervals.
fprintf('lower_interval=%f, upper_interval=%f', confidence_interval{1,1}, confidence_interval{1,2});
lower_interval=confidence_interval{1,1}
upper_interval=confidence_interval{1,2}
%% Depict graphs
% Edit the parameters below to depict graphs shown in <https://kazkobara.github.io/ebcic/ 
% this page>.
% 
% *Troubleshooting*
% 
% If errors like "Can't find a usable xxx in the following directories" are 
% displayed, fix the path so that missing xxx can be found in the path, e.g. for 
% Python 3.8:
%% 
% * 'init.tcl' is in 'Python38\tcl\tcl8.6\', but expected to be in 'Python38\lib\tcl8.6\'
% * 'tk.tcl' is in 'Python38\tcl\tk8.6', but expected to be in 'Python38\lib\tk8.6\' 
% or 'Python38\lib\tcl8.6\tk8.6'. 
%% 
% To show another graph, the popped-up window with the current graph has to 
% be closed.
% 
% Legend outside of the graph area is not displayed in MATLAB. So it seems better 
% to use Python directly to depict graphs. :-)

py.ebcic.interval_graph(py.ebcic.GraProps(pyargs( ...
    ...% Set the range of k to depict with k_*
    'k_start',int32(1), ...
    'k_end',  int32(1), ...
    ...% Edit the list of confidence percentages, [confi_perc, ...], to depict
    ...% where 0 < confi_perc < 100 for two-sided.
    ...%
    ...% NOTE For one-sided, set confi_perc=(2 * confi_perc_for_one_sided - 100)
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
%% 
% This file is maintained at <https://github.com/KazKobara/ebcic https://github.com/KazKobara/ebcic> 
% with the <https://github.com/KazKobara/ebcic/blob/master/LICENSE MIT License>.
% 
% When you use or publish the confidence interval obtained with the software, 
% please refer to the software name, version, platform, and so on, so that readers 
% can verify the correctness and reproducibility of the interval with the input 
% parameters.
% 
% An example of the reference is:
% 
% |"The confidence interval is obtained by EBCIC X.X.X on Python 3."|
% 
% where EBCIC is the name of the software, and X.X.X is the version of it.
% 
%