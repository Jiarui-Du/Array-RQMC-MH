function file_index = array_rqmc_index_result_files(results_dir, recursive)
%ARRAY_RQMC_INDEX_RESULT_FILES Index individual target--test MAT files.
%
% The function scans files following
%
%   array_rqmc__target-<target_id>__f-<test_id>__*.mat
%
% and reads only the metadata required for candidate selection.
%
% Returns candidate metadata and load status, including incomplete or invalid
% files. Filename IDs may be overridden by embedded settings. It loads the
% results variable to inspect metadata and available variance pairs.


    %% Discover matching result files
    if nargin < 2
        recursive = false;
    end

    results_dir = char(results_dir);
    pattern = 'array_rqmc__target-*__f-*__*.mat';
    if recursive
        files = dir(fullfile(results_dir, '**', pattern));
    else
        files = dir(fullfile(results_dir, pattern));
    end


    %% Allocate the candidate index columns
    TargetID = cell(0, 1);
    TestFunctionID = cell(0, 1);
    Path = cell(0, 1);
    Replications = zeros(0, 1);
    NMaxCompleted = zeros(0, 1);
    ValidVariancePairs = zeros(0, 1);
    NValues = cell(0, 1);
    ModifiedDatenum = zeros(0, 1);
    Status = cell(0, 1);

    row = 0;

    %% Inspect saved settings, completed counts and variance availability
    for file_number = 1:numel(files)
        file_path = fullfile(files(file_number).folder, files(file_number).name);
        parsed = regexp(files(file_number).name, ...
            '^array_rqmc__target-(.+?)__f-(.+?)__', 'tokens', 'once');
        if isempty(parsed)
            continue;
        end

        target_id = parsed{1};
        test_id = parsed{2};
        replications = NaN;
        n_max = NaN;
        valid_pairs = 0;
        n_values = [];
        status = 'indexed';

        try
            loaded = load(file_path, 'results');
            if ~isfield(loaded, 'results') || ~isstruct(loaded.results)
                status = 'missing results struct';
            else
                results = loaded.results;

                if isfield(results, 'settings') && isstruct(results.settings)
                    if isfield(results.settings, 'target_id') && ...
                            ~isempty(results.settings.target_id)
                        target_id = char(results.settings.target_id);
                    end
                    if isfield(results.settings, 'test_function_id') && ...
                            ~isempty(results.settings.test_function_id)
                        test_id = char(results.settings.test_function_id);
                    end
                    if isfield(results.settings, 'R') && ...
                            isnumeric(results.settings.R) && ...
                            isscalar(results.settings.R) && ...
                            isfinite(results.settings.R)
                        replications = double(results.settings.R);
                    end
                end

                if isfield(results, 'N_vec')
                    N = double(results.N_vec(:));
                    completed = numel(N);
                    if isfield(results, 'completed_idx') && ...
                            isnumeric(results.completed_idx) && ...
                            isscalar(results.completed_idx) && ...
                            isfinite(results.completed_idx)
                        completed = min(completed, ...
                            max(0, floor(double(results.completed_idx))));
                    end
                    N = N(1:completed);
                    n_values = N(isfinite(N) & N > 0);
                    if ~isempty(n_values)
                        n_max = max(n_values);
                    end
                end

                [variance_iid, variance_array] = extract_variance_fields(results);
                if ~isempty(variance_iid) && ~isempty(variance_array)
                    length_common = min([numel(n_values), ...
                        numel(variance_iid), numel(variance_array)]);
                    if length_common > 0
                        N_pair = n_values(1:length_common);
                        variance_iid = variance_iid(1:length_common);
                        variance_array = variance_array(1:length_common);
                        valid = isfinite(N_pair) & N_pair > 0 & ...
                            isfinite(variance_iid) & variance_iid > 0 & ...
                            isfinite(variance_array) & variance_array > 0;
                        valid_pairs = nnz(valid);
                    end
                end
            end
        catch ME
            status = ['load failed: ', ME.message];
        end

        row = row + 1;
        TargetID{row, 1} = target_id;
        TestFunctionID{row, 1} = test_id;
        Path{row, 1} = file_path;
        Replications(row, 1) = replications;
        NMaxCompleted(row, 1) = n_max;
        ValidVariancePairs(row, 1) = valid_pairs;
        NValues{row, 1} = n_values;
        ModifiedDatenum(row, 1) = files(file_number).datenum;
        Status{row, 1} = status;
    end


    %% Return the candidate table for deterministic selection
    file_index = table(TargetID, TestFunctionID, Path, Replications, ...
        NMaxCompleted, ValidVariancePairs, NValues, ModifiedDatenum, Status);
end


%% Local helper: extract_variance_fields
% Read the current variance fields or their supported historical aliases.
function [variance_iid, variance_array] = extract_variance_fields(results)
    variance_iid = [];
    variance_array = [];

    if isfield(results, 'var_mcmc')
        variance_iid = double(results.var_mcmc(:));
    end

    if isfield(results, 'var_sobol')
        variance_array = double(results.var_sobol(:));
    elseif isfield(results, 'var_owennus')
        variance_array = double(results.var_owennus(:));
    elseif isfield(results, 'var_nus')
        variance_array = double(results.var_nus(:));
    elseif isfield(results, 'var_original')
        variance_array = double(results.var_original(:));
    end
end
