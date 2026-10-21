function logp0 = array_rqmc_validate_handles( ...
        log_pi, test_functions, x0)
%ARRAY_RQMC_VALIDATE_HANDLES Check vectorized target and test handles.
%
% Checks two-element vector calls and returns log_pi(x0). Negative infinity
% is permitted away from x0; the initial density and test probe values must
% be finite. These probes do not certify behavior at every possible state.


    %% Check vectorized log-density evaluation
    probe = [x0; x0];
    try
        log_values = log_pi(probe);
    catch ME
        error('array_rqmc:LogPiEvaluation', ...
            'LogPi failed on a two-element vector: %s', ME.message);
    end
    if ~isnumeric(log_values) || numel(log_values) ~= 2 || ...
            ~isreal(log_values) || any(isnan(log_values(:))) || ...
            any(log_values(:) == Inf)
        error('array_rqmc:LogPiVectorized', ...
            'LogPi must return one real non-NaN value per input.');
    end


    %% Require a finite starting log density
    logp0 = double(log_pi(x0));
    if ~isscalar(logp0) || ~isreal(logp0) || ~isfinite(logp0)
        error('array_rqmc:InitialStateDensity', ...
            'LogPi(InitialState) must be a finite real scalar.');
    end


    %% Check each vectorized test function
    for test_index = 1:numel(test_functions)
        try
            values = test_functions{test_index}(probe);
        catch ME
            error('array_rqmc:TestFunctionEvaluation', ...
                'Test function %d failed: %s', test_index, ME.message);
        end
        if ~isnumeric(values) || numel(values) ~= 2 || ...
                ~isreal(values) || any(~isfinite(values(:)))
            error('array_rqmc:TestFunctionVectorized', ...
                'Test function %d must return finite real values.', ...
                test_index);
        end
    end
end
