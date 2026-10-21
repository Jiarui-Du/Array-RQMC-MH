function results = compare_array_rqmc_gaussian_rwm(varargin)
%COMPARE_ARRAY_RQMC_GAUSSIAN_RWM Compatibility statistical entry point.
%
% The implementation lives in experiments/. Existing callers can continue
% using this public function and its historical name-value interface.
%
% Forwards name/value options to the joint statistical experiment and returns
% its raw estimates and summaries; separate timing is not performed here.


    %% Register modules and dispatch to the implementation
    array_rqmc_setup();
    results = run_array_rqmc_statistical_experiment(varargin{:});
end
