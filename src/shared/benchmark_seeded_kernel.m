function timing=benchmark_seeded_kernel(runOne,N,methods,master,caseID,statisticalSeeds,options)
%BENCHMARK_SEEDED_KERNEL Apply the common protocol with disjoint recorded seeds.
% runOne takes (N,methodName,seed); target/model parameters are already resolved.
    arguments
        runOne (1,1) function_handle
        N (1,:) double
        methods (1,2) cell
        master (1,1) double
        caseID (1,1) double
        statisticalSeeds double
        options.Replications (1,1) double {mustBeInteger,mustBePositive} = 20
        options.WarmupExponent (1,1) double {mustBeInteger,mustBePositive} = 12
    end

    %% Resolve and check the exact seed sets before starting any timing
    R=options.Replications;assert(R<=499);
    seeds=zeros(R,numel(N),2);warmupSeeds=zeros(numel(N),2);
    for ni=1:numel(N)
        for mi=1:2
            warmupSeeds(ni,mi)=standalone_timing_seed(master,caseID, ...
                min(N(ni),2^options.WarmupExponent),mi,1,true);
            for r=1:R
                seeds(r,ni,mi)=standalone_timing_seed(master,caseID,N(ni),mi,r,false);
            end
        end
    end
    assert(isempty(intersect([seeds(:);warmupSeeds(:)],statisticalSeeds(:))), ...
        'Timing and statistical seeds overlap; change the master seed');
    assert(numel(unique(seeds(:)))==numel(seeds));

    %% Time serial complete calls with separate warmup randomization
    callback=@(n,mi,r,warmup)runOne(n,methods{mi}, ...
        standalone_timing_seed(master,caseID,n,mi,r,warmup));
    timing=benchmark_single_test(callback,N,methods, ...
        'Replications',R,'WarmupExponent',options.WarmupExponent);
    timing.seeds=seeds;timing.warmupSeeds=warmupSeeds;
    timing.seedFormula='mod(master+2^31+caseID*1e6+log2(N)*1e4+methodID*1e3+replication+500*isWarmup,2^32)';
end
