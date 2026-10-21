function [test_functions, test_ids, test_labels, reference_values] = ...
        array_rqmc_resolve_tests(tests, target, x0)
%ARRAY_RQMC_RESOLVE_TESTS Resolve target-dependent test definitions.
%
% Builds vectorized handles and aligned metadata for one catalog target.
% A zero reference is supplied only for an odd test, symmetric target and zero
% start; other finite-horizon reference values remain unspecified (NaN).


    %% Allocate handles, references and labels
    test_functions = cell(1, numel(tests));
    reference_values = nan(1, numel(tests));
    test_ids = {tests.id};
    test_labels = {tests.label};


    %% Resolve target-dependent CDF thresholds and symmetry references
    for test_index = 1:numel(tests)
        test = tests(test_index);
        switch test.kind
            case 'fixed_handle'
                test_functions{test_index} = test.handle;
            case 'cdf_indicator'
                threshold = double(target.cdf_threshold);
                test_functions{test_index} = @(x) double(x <= threshold);
            otherwise
                error('array_rqmc:UnknownTestKind', ...
                    'Unknown test kind ''%s''.', test.kind);
        end

        if target.is_symmetric && test.is_odd && x0 == 0
            reference_values(test_index) = 0;
        end
    end
end
