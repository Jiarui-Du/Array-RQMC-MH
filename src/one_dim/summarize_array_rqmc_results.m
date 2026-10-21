function results = summarize_array_rqmc_results(varargin)
%SUMMARIZE_ARRAY_RQMC_RESULTS Compatibility statistics entry point.
%
% Forwards raw arrays and metadata to the statistics implementation.
% The returned result struct is assembled without generating new trajectories.


    %% Register modules and dispatch to the implementation
    array_rqmc_setup();
    results = summarize_array_rqmc_statistics(varargin{:});
end
