function timing = benchmark_array_rqmc_single_test_times(varargin)
%BENCHMARK_ARRAY_RQMC_SINGLE_TEST_TIMES Compatibility timing entry point.
%
% Forwards name/value options to the serial standalone timing implementation.
% The returned timing struct contains per-test samples and median elapsed times.


    %% Register modules and dispatch to the implementation
    array_rqmc_setup();
    timing = run_array_rqmc_timing_experiment(varargin{:});
end
