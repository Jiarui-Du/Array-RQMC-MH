function [N, variance_iid, variance_array, array_label] = ...
        array_rqmc_extract_variance_data(results, positive_only)
%ARRAY_RQMC_EXTRACT_VARIANCE_DATA Read current or legacy variance fields.
%
% Returns aligned column vectors, sorted by N, from supported saved schemas.
% Truncates to completed_idx when present; positive_only additionally excludes
% pairs unsuitable for logarithmic plots.


    %% Resolve variance fields across result schema versions
    if nargin < 2
        positive_only = false;
    end
    if ~isfield(results, 'N_vec') || ~isfield(results, 'var_mcmc')
        error('array_rqmc:ResultSchema', 'Missing N_vec or var_mcmc.');
    end
    N = double(results.N_vec(:));
    variance_iid = double(results.var_mcmc(:));
    if isfield(results, 'var_sobol')
        variance_array = double(results.var_sobol(:));
        array_label = 'Array-RQMC';
    elseif isfield(results, 'var_owennus')
        variance_array = double(results.var_owennus(:));
        array_label = 'Array-Owen-NUS';
    elseif isfield(results, 'var_nus')
        variance_array = double(results.var_nus(:));
        array_label = 'Array-Owen-NUS';
    elseif isfield(results, 'var_original')
        variance_array = double(results.var_original(:));
        array_label = 'Array-RQMC';
    else
        error('array_rqmc:ResultSchema', ...
            'No supported Array-RQMC variance field.');
    end


    %% Restrict arrays to the recorded completed range
    completed = min([numel(N), numel(variance_iid), numel(variance_array)]);
    if isfield(results, 'completed_idx') && ...
            isnumeric(results.completed_idx) && ...
            isscalar(results.completed_idx) && isfinite(results.completed_idx)
        completed = min(completed, max(0, floor(results.completed_idx)));
    end
    N = N(1:completed);
    variance_iid = variance_iid(1:completed);
    variance_array = variance_array(1:completed);


    %% Optionally filter log-plot data and sort paired values
    if positive_only
        valid = isfinite(N) & N > 0 & isfinite(variance_iid) & ...
            variance_iid > 0 & isfinite(variance_array) & variance_array > 0;
        N = N(valid);
        variance_iid = variance_iid(valid);
        variance_array = variance_array(valid);
    end
    [N, order] = sort(N, 'ascend');
    variance_iid = variance_iid(order);
    variance_array = variance_array(order);
end
