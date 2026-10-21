function results = array_rqmc_attach_single_test_timing( ...
        results, joint_timing, test_index)
%ARRAY_RQMC_ATTACH_SINGLE_TEST_TIMING Attach one standalone benchmark.
%
% Extracts one timing column and attaches it to a single-test result.
% TAE is formed only when the benchmark N occurs in that result grid and
% both its saved VRF and elapsed-time ratio are finite.


    %% Extract the selected standalone benchmark and metadata
    timing = struct;
    scalar_fields = {'N','m','R','T0','T','sigma','initial_state', ...
        'master_seed','warmup_N','warmup_m'};
    for field_index = 1:numel(scalar_fields)
        name = scalar_fields{field_index};
        timing.(name) = double(joint_timing.(name));
    end
    timing.test_id = joint_timing.test_ids{test_index};
    timing.test_label = joint_timing.test_labels{test_index};
    timing.time_mcmc = joint_timing.time_mcmc(:, test_index);
    timing.time_sobol = joint_timing.time_sobol(:, test_index);
    timing.median_time_mcmc = joint_timing.median_time_mcmc(test_index);
    timing.median_time_sobol = joint_timing.median_time_sobol(test_index);
    timing.time_ratio_iid_to_array = ...
        joint_timing.time_ratio_iid_to_array(test_index);
    timing.execution = joint_timing.execution;
    timing.method_order = joint_timing.method_order;
    timing.scope = joint_timing.scope;
    if isfield(joint_timing,'protocol')
        timing.protocol=joint_timing.protocol;
        timing.computationalThreads=joint_timing.computationalThreads;
    end
    timing.vrf_sobol = NaN;
    timing.time_adjusted_efficiency_sobol = NaN;


    %% Match the timing particle count and compute its TAE
    position = find(abs(double(results.N_vec(:)) - timing.N) <= ...
        8 * eps(timing.N), 1, 'last');
    if ~isempty(position) && numel(results.vrf_sobol) >= position
        timing.vrf_sobol = double(results.vrf_sobol(position));
        if isfinite(timing.vrf_sobol) && ...
                isfinite(timing.time_ratio_iid_to_array)
            timing.time_adjusted_efficiency_sobol = ...
                timing.vrf_sobol * timing.time_ratio_iid_to_array;
        end
    end

    %% Attach timing provenance to the result
    results.timing = timing;
    results.settings.timing_status = 'completed';
    results.settings.timing_scope = timing.scope;
    results.settings.timing_N = timing.N;
    results.settings.timing_R = timing.R;
end
