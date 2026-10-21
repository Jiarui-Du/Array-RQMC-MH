function manifest = run_array_rqmc_suite(varargin)
%RUN_ARRAY_RQMC_SUITE Compatibility paper-suite entry point.
%
% Statistical simulations accumulate all selected tests jointly. Standalone
% timing remains a separate serial experiment per test function.
%
% Forwards name/value options and returns the target/test output manifest.
% This wrapper can start simulations; plotting and export wrappers do not call it.


    %% Register modules and dispatch to the implementation
    array_rqmc_setup();
    manifest = run_array_rqmc_paper_suite(varargin{:});
end
