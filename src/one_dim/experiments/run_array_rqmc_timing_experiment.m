function timing = run_array_rqmc_timing_experiment(varargin)
%RUN_ARRAY_RQMC_TIMING_EXPERIMENT Time complete single-test runs serially.
%
% Every recorded interval contains a complete algorithm run. No time is
% allocated from the joint statistical simulation.
%
% Returns R-by-test elapsed-time arrays and their per-test medians.
% Timing streams use separate method offsets from statistical replications;
% alternating method order reduces systematic first-run ordering effects.


    %% Resolve tests, benchmark size and initial state
    opt = array_rqmc_parse_timing_options(varargin{:});
    if isa(opt.TestFunctions, 'function_handle')
        test_functions = {opt.TestFunctions};
    else
        test_functions = reshape(opt.TestFunctions, 1, []);
    end
    test_count = numel(test_functions);
    test_ids = array_rqmc_normalize_text_metadata( ...
        opt.TestFunctionIDs, test_count, 'f');
    test_labels = array_rqmc_normalize_text_metadata( ...
        opt.TestFunctionLabels, test_count, 'f');

    m = double(opt.m);
    N = 2^m;
    T0 = double(opt.T0);
    T = double(opt.T);
    R = double(opt.R);
    sigma = double(opt.sigma);
    x0 = double(opt.InitialState);
    log_pi = opt.LogPi;
    logp0 = array_rqmc_validate_handles(log_pi, test_functions, x0);
    sobol_base = array_rqmc_sobol_base();
    if isempty(opt.WarmupM)
        warmup_m = min(m, 12);
    else
        warmup_m = min(m, double(opt.WarmupM));
    end
    warmup_N = 2^warmup_m;


    %% Apply the shared serial protocol separately to each test
    time_iid=nan(R,test_count);time_array=nan(R,test_count);
    benchmarks=cell(1,test_count);
    for test_index=1:test_count
        fprintf('\nTiming: test=%s\n',test_ids{test_index});
        current_test=test_functions(test_index);
        callback=@(n,mi,r,warmup)run_timed_kernel(n,mi,r,warmup,m, ...
            opt.MasterSeed,T0,T,sigma,x0,logp0,log_pi,current_test,sobol_base);
        b=benchmark_single_test(callback,N,{'iid','array'}, ...
            'Replications',R,'WarmupExponent',warmup_m);
        time_iid(:,test_index)=b.seconds(:,1,1);
        time_array(:,test_index)=b.seconds(:,1,2);
        benchmarks{test_index}=b;
    end

    %% Return timing samples, medians and measurement provenance
    timing = struct( ...
        'N', N, 'm', m, 'R', R, 'T0', T0, 'T', T, ...
        'sigma', sigma, 'initial_state', x0, ...
        'master_seed', double(opt.MasterSeed), ...
        'test_ids', {test_ids}, 'test_labels', {test_labels}, ...
        'time_mcmc', time_iid, 'time_sobol', time_array, ...
        'median_time_mcmc', median(time_iid, 1), ...
        'median_time_sobol', median(time_array, 1), ...
        'time_ratio_iid_to_array', ...
            median(time_iid, 1) ./ median(time_array, 1), ...
        'warmup_N', warmup_N, 'warmup_m', warmup_m, ...
        'execution', 'serial', ...
        'method_order', ...
            'IID first on odd replications; Array-RQMC first on even replications.', ...
        'scope', b.scope, 'protocol',b.protocol, ...
        'computationalThreads',b.computationalThreads,'benchmarks',{benchmarks});
end

%% Local helper: run_timed_kernel
% Stream construction is inside the same timed boundary as multidimensional runs.
function out=run_timed_kernel(N,mi,r,warmup,m,seed,T0,T,sigma,x0,logp0,logpi,test,base)
    substream=m*1e6+r;if warmup,substream=m*1e6+900001;end
    offsets=[211 237];stream=array_rqmc_make_method_stream(seed,offsets(mi),substream);
    if mi==1
        out=run_iid_gaussian_rwm(N,T0,T,sigma,x0,logp0,logpi,test,stream);
    else
        out=run_array_rqmc_gaussian_rwm(N,T0,T,sigma,x0,logp0,logpi,test,base,stream);
    end
end
