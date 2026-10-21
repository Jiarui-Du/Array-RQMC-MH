function single = select_array_rqmc_test_result(combined, test_selector)
%SELECT_ARRAY_RQMC_TEST_RESULT Extract one integrand from a joint result.
%
% SINGLE = SELECT_ARRAY_RQMC_TEST_RESULT(COMBINED, TEST_SELECTOR) accepts a
% numeric test index or a test-function ID. The output uses the historical
% single-integrand schema, so existing plotting and export workflows can
% read one MAT file per target--test pair.
%
% Raw estimates are sliced on their test axis; summaries are sliced by column.
% Shared acceptance data are retained, and settings preserve the original joint
% test collection as provenance alongside the selected single-test metadata.


    %% Validate the joint result and resolve the requested test
    if ~isstruct(combined) || ~isfield(combined, 'test_ids')
        error('select_array_rqmc_test_result:InvalidInput', ...
            'COMBINED must be a joint result containing test_ids.');
    end

    test_ids = normalize_text_cell(combined.test_ids);
    number_of_tests = numel(test_ids);

    if isnumeric(test_selector)
        if ~isscalar(test_selector) || ~isfinite(test_selector) || ...
                test_selector ~= floor(test_selector) || ...
                test_selector < 1 || test_selector > number_of_tests
            error('select_array_rqmc_test_result:TestIndex', ...
                'The numeric test selector is out of range.');
        end
        test_index = double(test_selector);
    else
        requested_id = char(test_selector);
        test_index = find(strcmpi(requested_id, test_ids), 1);
        if isempty(test_index)
            error('select_array_rqmc_test_result:UnknownTestID', ...
                'Unknown test-function ID ''%s''.', requested_id);
        end
    end


    %% Extract raw estimates, summaries and slopes for one test
    single = combined;
    % Full paired bootstrap draws retain their original output order.
    if isfield(single,'bootstrap'),single.bootstrap_selected_output=test_index;end

    raw_fields = {'est_mcmc', 'est_sobol'};
    for field_index = 1:numel(raw_fields)
        name = raw_fields{field_index};
        if isfield(single, name)
            values = single.(name);
            if ndims(values) >= 3
                single.(name) = values(:, :, test_index);
            end
        end
    end

    column_fields = { ...
        'mean_mcmc', 'mean_sobol', 'var_mcmc', 'var_sobol', ...
        'vrf_sobol', 'vrf_ci_low', 'vrf_ci_high', ...
        'bias_mcmc', 'bias_sobol', ...
        'empirical_mse_mcmc', 'empirical_mse_sobol'};
    for field_index = 1:numel(column_fields)
        name = column_fields{field_index};
        if isfield(single, name)
            values = single.(name);
            if ismatrix(values) && size(values, 2) >= test_index
                single.(name) = values(:, test_index);
            end
        end
    end

    if isfield(single, 'slopes') && isstruct(single.slopes)
        if isfield(single.slopes, 'mcmc') && ...
                numel(single.slopes.mcmc) >= test_index
            single.slopes.mcmc = single.slopes.mcmc(test_index);
        end
        if isfield(single.slopes, 'sobol') && ...
                numel(single.slopes.sobol) >= test_index
            single.slopes.sobol = single.slopes.sobol(test_index);
        end
    end


    %% Update single-test metadata while retaining joint provenance
    test_labels = normalize_text_cell(combined.test_labels);
    reference_values = reshape(double(combined.reference_values), 1, []);

    single.test_ids = test_ids(test_index);
    single.test_labels = test_labels(test_index);
    single.reference_values = reference_values(test_index);

    if isfield(single, 'settings') && isstruct(single.settings)
        single.settings.test_function_id = test_ids{test_index};
        single.settings.test_function_label = test_labels{test_index};
        single.settings.reference_value = reference_values(test_index);

        if isfield(combined.settings, 'test_functions') && ...
                numel(combined.settings.test_functions) >= test_index
            all_functions = normalize_text_cell( ...
                combined.settings.test_functions);
            single.settings.test_function = ...
                all_functions{test_index};
            single.settings.joint_test_functions = all_functions;
            single.settings.test_functions = all_functions(test_index);
        end

        single.settings.joint_test_ids = test_ids;
        single.settings.joint_test_labels = test_labels;
        single.settings.joint_number_of_test_functions = number_of_tests;
        single.settings.number_of_test_functions = 1;
        single.settings.test_ids = test_ids(test_index);
        single.settings.test_labels = test_labels(test_index);
        single.settings.reference_values = reference_values(test_index);
    end
end


%% Local helper: normalize_text_cell
% Convert text containers to a row cell array for consistent metadata indexing.
function values = normalize_text_cell(values)
    if ischar(values)
        values = {values};
    elseif isstring(values)
        values = cellstr(values(:)).';
    else
        values = values(:).';
        values = cellfun(@char, values, 'UniformOutput', false);
    end
end
