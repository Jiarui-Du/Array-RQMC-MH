function summary = export_array_rqmc_vrf_tae_slopes_v2( ...
        results_dir, varargin)
%EXPORT_ARRAY_RQMC_VRF_TAE_SLOPES_V2 Export VRF, TaE, and slopes.
%
% TaE = VRF * median(IID standalone time) / median(Array standalone time).
% Current and legacy result schemas are supported.
%
% Returns the combined target/test table and writes CSV plus TeX artifacts.
% Strict mode reports incomplete cases after exporting their diagnostic rows.
% Existing export paths receive a timestamp suffix through the path helper.


    %% Resolve options and index available results
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
    file_index = array_rqmc_index_result_files(results_dir, opt.Recursive);
    if isempty(file_index)
        error('array_rqmc:NoResultFiles', ...
            'No result MAT files found in %s.', results_dir);
    end


    %% Summarize each requested target/test pair
    target_ids = array_rqmc_normalize_id_selection(opt.TargetIDs);
    test_ids = array_rqmc_normalize_id_selection(opt.TestFunctionIDs);
    row_cells = cell(numel(target_ids) * numel(test_ids), 1);
    problems = cell(0, 1);
    row_index = 0;
    for target_index = 1:numel(target_ids)
        for test_index = 1:numel(test_ids)
            [row, problem] = array_rqmc_export_case( ...
                file_index, target_ids{target_index}, ...
                test_ids{test_index}, opt);
            row_index = row_index + 1;
            row_cells{row_index} = row;
            if ~isempty(problem)
                problems{end+1, 1} = sprintf('%s/%s: %s', ...
                    target_ids{target_index}, test_ids{test_index}, ...
                    problem); %#ok<AGROW>
            end
        end
    end

    %% Normalize table text columns
    summary = struct2table(vertcat(row_cells{:}));
    text_variables = {'TargetID','Target','TestFunctionID','TestFunction', ...
        'TaESource','TimingScope','Status','ResultFile'};
    for variable_index = 1:numel(text_variables)
        name = text_variables{variable_index};
        if ~iscell(summary.(name))
            summary.(name) = cellstr(summary.(name));
        end
    end


    %% Export tables and enforce strict completeness if requested
    summary.Status=[];
    if ~isempty(opt.OutputCSV)
        output_csv=array_rqmc_unique_output_path(char(opt.OutputCSV));
        ensure_parent(output_csv);writetable(summary,output_csv);
        fprintf('CSV: %s\n',output_csv);
    end
    if ~isempty(opt.OutputTeX)
        output_tex=array_rqmc_unique_output_path(char(opt.OutputTeX));
        ensure_parent(output_tex);
        write_array_rqmc_latex_table(output_tex,summary,target_ids,test_ids,opt);
        fprintf('LaTeX: %s\n',output_tex);
    end
    fprintf('\nArray-RQMC summary at N=%g:\n', double(opt.TargetN));
    disp(summary(:, {'Target','TestFunction','VRF','TaE', ...
        'ArraySlope','SlopeCILow','SlopeCIHigh','VRFCILow','VRFCIHigh'}));
    if opt.Strict && ~isempty(problems)
        error('array_rqmc:IncompleteExportResults', ...
            'Incomplete results:\n%s', strjoin(problems, '\n'));
    end
end

%% Local helper: parse_options
% Parse and validate the public name/value interface for this operation.
function opt = parse_options(results_dir, varargin)
    ip = inputParser;
    addParameter(ip, 'TargetIDs', {'normal','double_well_e1', ...
        'gaussian_mixture3','student_t3'}, @id_selection);
    addParameter(ip, 'TestFunctionIDs', {'cdf','x','x2'}, @id_selection);
    addParameter(ip, 'TargetN', 2^20, @positive_integer);
    addParameter(ip, 'Bootstrap', 1000, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x>=0 && x==fix(x));
    addParameter(ip, 'TailCount', 5, @(x) positive_integer(x) && x >= 2);
    addParameter(ip, 'RequiredR', 200, @optional_positive_integer);
    addParameter(ip, 'OutputCSV', fullfile(results_dir, ...
        'vrf_tae_slope_summary.csv'), @text_scalar);
    addParameter(ip, 'OutputTeX', fullfile(results_dir, ...
        'table_all_vrf_tae_slopes.tex'), @text_scalar);
    addParameter(ip, 'Strict', false, @logical_scalar);
    addParameter(ip, 'Recursive', false, @logical_scalar);
    parse(ip, varargin{:});
    opt = ip.Results;
end

%% Local helper: ensure_parent
% Create the export parent directory when needed.
function ensure_parent(file_path)
    parent = fileparts(file_path);
    if ~isempty(parent) && ~isfolder(parent)
        [ok, message] = mkdir(parent);
        if ~ok
            error('array_rqmc:CreateDirectory', ...
                'Could not create %s: %s', parent, message);
        end
    end
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
