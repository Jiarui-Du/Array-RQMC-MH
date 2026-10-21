function [estimate_iid, estimate_array, acceptance_iid, acceptance_array] = ...
        array_rqmc_run_macro_replication( ...
            N, substream_key, T0, T, sigma, x0, logp0, log_pi, ...
            test_functions, sobol_base, master_seed, replication)
%ARRAY_RQMC_RUN_MACRO_REPLICATION Run paired method replications.
%
% The methods use independent mrg32k3a streams. All test functions are
% accumulated jointly along the same trajectory for each method.
%
% Returns one estimate vector and acceptance rate for each matching kernel.
% The same starting state, target, horizon and test collection are used by both;
% method offsets keep their random streams distinct.


    %% Assign deterministic streams for this macro replication
    % Reserve a particle-number namespace before the experiment-specific key.
    m = round(log2(N));
    substream = double(m) * 1e9 + ...
        double(substream_key) * 1e6 + double(replication);
    stream_iid = array_rqmc_make_method_stream(master_seed, 11, substream);
    stream_array = array_rqmc_make_method_stream(master_seed, 37, substream);


    %% Run both kernels on the same experimental configuration
    [estimate_iid, acceptance_iid] = run_iid_gaussian_rwm( ...
        N, T0, T, sigma, x0, logp0, log_pi, test_functions, stream_iid);
    [estimate_array, acceptance_array] = run_array_rqmc_gaussian_rwm( ...
        N, T0, T, sigma, x0, logp0, log_pi, test_functions, ...
        sobol_base, stream_array);
end
