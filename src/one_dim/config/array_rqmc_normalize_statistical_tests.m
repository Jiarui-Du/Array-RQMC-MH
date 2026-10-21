function [test_functions, test_ids, test_labels, reference_values] = ...
        array_rqmc_normalize_statistical_tests(opt)
%ARRAY_RQMC_NORMALIZE_STATISTICAL_TESTS Normalize comparison test options.
%
% Returns aligned handles, IDs, labels and reference values. With no explicit
% handles, uses the three standard-normal tests; unknown references remain NaN.


    %% Resolve mutually exclusive single-test and multi-test inputs
    has_multi = ~isempty(opt.TestFunctions);
    has_single = ~isempty(opt.TestFunction);
    if has_multi && has_single
        error('array_rqmc:TestFunctionOptions', ...
            'Use TestFunctions or TestFunction, not both.');
    end

    used_defaults = false;
    if has_multi
        if isa(opt.TestFunctions, 'function_handle')
            test_functions = {opt.TestFunctions};
        else
            test_functions = reshape(opt.TestFunctions, 1, []);
        end
    elseif has_single
        test_functions = {opt.TestFunction};
    else
        threshold = -0.6744897501960817;
        test_functions = {@(x) double(x <= threshold), @(x) x, @(x) x.^2};
        used_defaults = true;
    end


    %% Align labels and reference values with the selected handles
    count = numel(test_functions);
    if used_defaults && isempty(opt.TestFunctionIDs)
        test_ids = {'cdf', 'x', 'x2'};
    else
        test_ids = normalize_text(opt.TestFunctionIDs, count, 'f');
    end
    if used_defaults && isempty(opt.TestFunctionLabels)
        test_labels = {'f_CDF(x)=1{x<=a_pi}', 'f(x)=x', 'f(x)=x^2'};
    else
        test_labels = normalize_text(opt.TestFunctionLabels, count, 'f');
    end

    if isempty(opt.ReferenceValues)
        reference_values = repmat(double(opt.ReferenceValue), 1, count);
    else
        reference_values = reshape(double(opt.ReferenceValues), 1, []);
        if isscalar(reference_values) && count > 1
            reference_values = repmat(reference_values, 1, count);
        end
    end
    if numel(reference_values) ~= count
        error('array_rqmc:ReferenceValues', ...
            'ReferenceValues must contain one value per test function.');
    end
end

%% Local helper: normalize_text
% Generate default test names when omitted and enforce one label per test.
function values = normalize_text(value, count, prefix)
    if isempty(value)
        values = arrayfun(@(k) sprintf('%s%d', prefix, k), 1:count, ...
            'UniformOutput', false);
    elseif ischar(value)
        values = {value};
    elseif isstring(value)
        values = cellstr(value(:)).';
    else
        values = value(:).';
        values = cellfun(@char, values, 'UniformOutput', false);
    end
    if numel(values) ~= count
        error('array_rqmc:TestMetadata', ...
            'Test IDs and labels must match the test-function count.');
    end
end
