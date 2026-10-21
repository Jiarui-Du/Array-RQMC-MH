function timing=benchmark_single_test(runOne,N,methods,options)
%BENCHMARK_SINGLE_TEST Independent serial timings of complete single-test runs.
% runOne(N,methodIndex,replication,isWarmup) invokes one complete kernel.
% It must use randomization separate from statistics and restore global streams.
% Model/data preparation and output are outside the recorded interval; per-run
% random-stream setup, initialization, transitions and test evaluation are inside.
    arguments
        runOne (1,1) function_handle
        N (1,:) double {mustBeInteger,mustBePositive}
        methods (1,2) cell
        options.Replications (1,1) double {mustBeInteger,mustBePositive} = 20
        options.WarmupExponent (1,1) double {mustBeInteger,mustBePositive} = 12
    end

    %% Warm up each kernel and alternate method order within every N
    assert(all(log2(N)==fix(log2(N))) && all(diff(N)>0));
    R=options.Replications;seconds=nan(R,numel(N),2);
    warmupN=min(N,2^options.WarmupExponent);
    for ni=1:numel(N)
        fprintf('Timing N=2^%d: serial, %d repetitions per method; warming up.\n',log2(N(ni)),R);
        for mi=1:2,runOne(warmupN(ni),mi,1,true);end
        for r=1:R
            order=1:2;if mod(r,2)==0,order=2:-1:1;end
            for mi=order
                timer=tic;
                runOne(N(ni),mi,r,false);
                seconds(r,ni,mi)=toc(timer);
            end
            % Console I/O is strictly outside each measured kernel interval.
            if mod(r,5)==0 || r==R
                fprintf('Timing N=2^%d: %d/%d repetitions completed per method.\n',log2(N(ni)),r,R);
            end
        end
        med=reshape(median(seconds(:,ni,:),1),1,2);
        fprintf('Timing N=2^%d median: %s %.6f s; %s %.6f s.\n', ...
            log2(N(ni)),methods{1},med(1),methods{2},med(2));
    end

    %% Store raw intervals, aggregation and execution provenance
    timing=struct('protocol','serial-single-test-v1','N',N,'methods',{methods}, ...
        'replications',R,'seconds',seconds, ...
        'medianSeconds',reshape(median(seconds,1),numel(N),2), ...
        'warmupN',warmupN,'execution','serial','aggregation','median', ...
        'methodOrder','IID first on odd replications; Array first on even replications', ...
        'scope','Complete single-test kernel including per-run stream setup; excludes data/model preparation and file I/O', ...
        'matlabVersion',version,'computer',computer,'computationalThreads',maxNumCompThreads);
end
