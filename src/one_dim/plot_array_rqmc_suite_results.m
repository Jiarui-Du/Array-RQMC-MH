function [figures, summary, file_index] = ...
        plot_array_rqmc_suite_results(results_dir, varargin)
%PLOT_ARRAY_RQMC_SUITE_RESULTS Plot the complete target--test suite.
%
% The default output is one 1-by-3 figure for each target distribution.
% Each figure contains the three paper test functions
%
%   f_CDF(x),  f(x)=x,  f(x)=x^2.
%
% The function scans individual result MAT files directly and does not
% depend on array_rqmc_suite_manifest.mat.  The default output filenames
% identify the target explicitly:
%
%   variance_normal.pdf
%   variance_double_well_e1.pdf
%   variance_gaussian_mixture3.pdf
%   variance_student_t3.pdf
%
% Example:
% [figures, summary, file_index] = plot_array_rqmc_suite_results( ...
%     results_dir, ...
%     'RequiredR', 200, ...
%     'TailCount', 5);
%
% Optional alternative layouts remain available through
% FigureMode='combined' or FigureMode='paper-pairs'.
%
% Adds paper-layout defaults only when the caller has not supplied the option.
% Returns figures, slope summary and result index from the common plotting core.


    %% Resolve the result folder and fill missing plotting defaults
    if nargin < 1 || isempty(results_dir)
        results_dir = fullfile(pwd, 'array_rqmc_results');
    end

    args = add_default_option(varargin, 'TargetIDs', 'all');
    args = add_default_option(args, 'TestFunctionIDs', ...
        {'cdf', 'x', 'x2'});

    % One target distribution per figure is the paper default.
    args = add_default_option(args, 'FigureMode', 'per-target');
    args = add_default_option(args, 'RequiredR', 200);
    args = add_default_option(args, 'TailCount', 5);
    args = add_default_option(args, 'OutputDir', ...
        array_rqmc_figure_dir('suite'));


    %% Dispatch saved-result plotting
    [figures, summary, file_index] = ...
        plot_array_rqmc_results_core(results_dir, args{:});
end


%% Local helper: add_default_option
% Append a default only if the caller did not already supply that option name.
function args = add_default_option(args, name, value)
    names = args(1:2:end);
    present = any(cellfun(@(x) ...
        (ischar(x) || (isstring(x) && isscalar(x))) && ...
        strcmpi(char(x), name), names));
    if ~present
        args = [args, {name, value}];
    end
end
