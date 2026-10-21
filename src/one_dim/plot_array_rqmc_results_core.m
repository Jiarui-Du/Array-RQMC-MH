function [figures, summary, file_index] = ...
        plot_array_rqmc_results_core(results_dir, varargin)
%PLOT_ARRAY_RQMC_RESULTS_CORE Plot current and legacy result MAT files.
%
% Returns figure handles, fitted-slope metadata and the candidate-file index.
% Nonempty MetadataOutputDir enables CSV exports with SaveFigures. Strict mode raises after
% collecting missing-panel diagnostics; no new trajectories are generated.


    %% Resolve options and create requested export folders
    array_rqmc_setup();
    if nargin < 1 || isempty(results_dir)
        results_dir = fullfile(pwd, 'array_rqmc_results');
    end
    results_dir = char(results_dir);
    opt = parse_options(results_dir, varargin{:});
    if ~isfolder(results_dir)
        error('array_rqmc:ResultsDirNotFound', ...
            'Results directory not found: %s', results_dir);
    end
    output_dir = char(opt.OutputDir);
    metadata_output_dir = char(opt.MetadataOutputDir);
    if opt.SaveFigures && ~isfolder(output_dir)
        [ok, message] = mkdir(output_dir);
        if ~ok
            error('array_rqmc:CreateOutputDir', ...
                'Could not create %s: %s', output_dir, message);
        end
    end
    if opt.SaveFigures && ~isempty(metadata_output_dir) && ~isfolder(metadata_output_dir)
        [ok, message] = mkdir(metadata_output_dir);
        if ~ok
            error('array_rqmc:CreateOutputDir', ...
                'Could not create %s: %s', metadata_output_dir, message);
        end
    end


    %% Index saved files and resolve target/test selections
    file_index = array_rqmc_index_result_files(results_dir, opt.Recursive);
    if isempty(file_index)
        error('array_rqmc:NoResultFiles', ...
            'No Array-RQMC result MAT files found in %s.', results_dir);
    end
    catalog = array_rqmc_case_catalog();
    target_ids = array_rqmc_resolve_plot_selection(opt.TargetIDs, ...
        unique(file_index.TargetID, 'stable'), {catalog.targets.id}, 'target');
    test_ids = array_rqmc_resolve_plot_selection(opt.TestFunctionIDs, ...
        unique(file_index.TestFunctionID, 'stable'), {catalog.tests.id}, ...
        'test function');


    %% Draw panels and export their provenance tables
    [figures, rows, problems] = array_rqmc_draw_figures( ...
        file_index, target_ids, test_ids, opt, output_dir);
    summary = cell2table(rows, 'VariableNames', { ...
        'TargetID','TestFunctionID','Method','Slope','SlopeSE', ...
        'RSquared','TailPointCount','NMinUsed','NMaxUsed', ...
        'ResultFile','SelectionStatus','Replications','DataStatus'});
    if opt.SaveFigures && ~isempty(metadata_output_dir)
        summary_file = array_rqmc_unique_output_path( ...
            fullfile(metadata_output_dir, 'slope_summary.csv'));
        writetable(summary, summary_file);
        index_export = file_index;
        index_export.NValues = cellfun(@mat2str, index_export.NValues, ...
            'UniformOutput', false);
        index_file = array_rqmc_unique_output_path( ...
            fullfile(metadata_output_dir, 'result_file_index.csv'));
        writetable(index_export, index_file);
    end
    if opt.Strict && ~isempty(problems)
        error('array_rqmc:IncompletePlotResults', ...
            'Plotting encountered problems:\n%s', strjoin(problems, '\n'));
    end
end

%% Local helper: parse_options
% Parse and validate the public name/value interface for this operation.
function opt = parse_options(results_dir, varargin)
    function_root = fileparts(mfilename('fullpath'));
    default_figure_dir = fullfile( ...
        fileparts(fileparts(function_root)), 'figures', 'selected');
    ip = inputParser;
    addParameter(ip, 'TargetIDs', 'all', @id_selection);
    addParameter(ip, 'TestFunctionIDs', {'cdf','x','x2'}, @id_selection);
    addParameter(ip, 'RequiredR', [], @optional_positive_integer);
    addParameter(ip, 'OutputDir', default_figure_dir, @text_scalar);
    addParameter(ip, 'MetadataOutputDir', results_dir, @text_scalar);
    addParameter(ip, 'SaveFigures', true, @logical_scalar);
    addParameter(ip, 'TailCount', 5, @(x) positive_integer(x) && x >= 2);
    addParameter(ip, 'CloseAfterSave', false, @logical_scalar);
    addParameter(ip, 'FigureVisible', 'on', ...
        @(x) any(strcmpi(char(x), {'on','off'})));
    addParameter(ip, 'Strict', false, @logical_scalar);
    addParameter(ip, 'Recursive', false, @logical_scalar);
    addParameter(ip, 'FigureMode', 'per-target', ...
        @(x) any(strcmpi(char(x), {'per-target','combined','paper-pairs'})));
    parse(ip, varargin{:});
    opt = ip.Results;
end

%% Local input validators
% Predicates used by the option parser; they do not alter supplied values.
function tf = id_selection(x), tf = ischar(x) || isstring(x) || iscell(x); end
function tf = positive_integer(x), tf = isnumeric(x) && isscalar(x) && ...
    isfinite(x) && x == floor(x) && x >= 1; end
function tf = optional_positive_integer(x), tf = isempty(x) || positive_integer(x); end
function tf = text_scalar(x), tf = ischar(x) || ...
    (isstring(x) && isscalar(x)); end
function tf = logical_scalar(x), tf = islogical(x) && isscalar(x); end
