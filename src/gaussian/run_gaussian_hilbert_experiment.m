function results = run_gaussian_hilbert_experiment(options)
%RUN_GAUSSIAN_HILBERT_EXPERIMENT Fixed-T variance vs N, independent macro runs.
% Example: run_gaussian_hilbert_experiment('Dimensions',[1 2 4 8], ...
%   'NExponents',8:16,'Replications',200,'T',128,'BurnIn',8)
%
% Returns results.cases and derived tables. Raw estimate axes are replication,
% N, test, method (time averages only). Checkpoints are
% written after each N; this runner does not load checkpoints to resume.

    %% Declare experiment and analysis options
    arguments
        options.Dimensions (1,:) double {mustBeInteger,mustBePositive} = [1 2 4 8]
        options.NExponents (1,:) double {mustBeInteger,mustBePositive} = 8:14
        options.Replications (1,1) double {mustBeInteger,mustBePositive} = 128
        options.T (1,1) double {mustBeInteger,mustBePositive} = 64
        options.BurnIn (1,1) double {mustBeInteger,mustBeNonnegative} = 8
        options.Seed (1,1) double {mustBeInteger,mustBeNonnegative} = 2026
        options.Target (1,:) char {mustBeMember(options.Target,{'standard','general'})} = 'general'
        options.Rho (1,1) double {mustBeFinite,mustBeReal} = 0.5
        options.Mean double = []
        options.Covariance double = []
        options.StepScale double = []
        options.Initialization (1,:) char = 'stationary'
        options.InitialState double = []
        options.SortCoordinates (1,:) char = 'whitened'
        options.HilbertBits double = []
        options.TailCount (1,1) double {mustBeInteger,mustBePositive} = 5
        options.Bootstrap (1,1) double {mustBeInteger,mustBeNonnegative} = 1000
        options.OutputDir (1,:) char = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'results','gaussian')
        options.RunTiming (1,1) logical = true
        options.TimingNExponents (1,:) double {mustBeInteger,mustBePositive} = []
        options.TimingReplications (1,1) double {mustBeInteger,mustBePositive} = 20
        options.TimingWarmupExponent (1,1) double {mustBeInteger,mustBePositive} = 12
        options.UseParallel (1,1) logical = true
        options.NumWorkers (1,1) double {mustBeInteger,mustBePositive} = 4
    end

    runTimer=tic;
    if options.RunTiming
        assert(options.TimingReplications<=499,'Use at most 499 independent timing samples');
        select_timing_N(2.^options.NExponents,options.TimingNExponents);
    end

    %% Validate initialization and grid compatibility
    assert(ismember(options.Initialization,{'stationary','zero','fixed'}),'Unknown initialization');
    if strcmp(options.Initialization,'zero')
        assert(isempty(options.InitialState),'InitialState must be empty when Initialization is zero');
    end
    if strcmp(options.Target,'standard')
        assert(isempty(options.Mean) && isempty(options.Covariance),...
            'For Target=standard, leave Mean and Covariance empty');
    end
    assert(options.Replications>=3 && options.T>options.BurnIn);
    assert(numel(options.NExponents)>=3 && all(diff(options.NExponents)>0));
    assert(numel(unique(options.Dimensions))==numel(options.Dimensions));
    assert(abs(options.Rho)<1 && options.TailCount>=3);
    assert(isempty(options.Covariance) || isscalar(options.Dimensions), ...
        'For a custom covariance, specify one dimension.');

    %% Configure methods, execution and output provenance
    R = options.Replications; Nvec = 2.^options.NExponents;
    methods = {'iid','hilbert'};
    methodIDs = [1 3]; % Preserve the existing IID/Hilbert random-stream namespaces.
    nMethods = numel(methods); nFunctions = 1;
    [useParallel,workerCount] = initialize_parallel(options);
    results = struct('schema','gaussian-joint-v2','settings',options,'matlabVersion',version, ...
        'methods',{methods},'cases',[]);
    results.execution = struct('useParallel',useParallel,'workerCount',workerCount, ...
        'implementationProfile','cached-hybrid-key-v3', ...
        'experimentProfile','projection-standalone-timing-v2','scheduler','single-parfor-per-N');
    if ~isfolder(options.OutputDir), mkdir(options.OutputDir); end
    tag = char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'));
    stem = fullfile(options.OutputDir,['gaussian_hilbert_' tag]);
    assert(~isfile([stem '.mat']) && ~isfile([stem '_checkpoint.mat']), ...
        'Output collision; rerun with another directory.');
    results.outputStem = stem;
    rows = {}; slopeRows = {};

    %% Construct each Gaussian target and run its particle grid
    for di = 1:numel(options.Dimensions)
        d = options.Dimensions(di);caseTimer=tic;
        fprintf('\n[%d/%d] dimension=%d; tests=projection\n',di,numel(options.Dimensions),d);
        print_experiment_header('Gaussian joint Hilbert statistical experiment',{'projection'}, ...
            options.NExponents,options.BurnIn,options.T,options.Replications,useParallel);
        bits=resolve_hilbert_bits(d,options.HilbertBits);
        mu = options.Mean;
        if strcmp(options.Target,'standard')
            mu=zeros(1,d);
        elseif isempty(mu)
            mu=linspace(-0.5,0.5,d);
        end
        assert(isequal(size(mu),[1,d]),'Mean must be 1-by-d.');
        Sigma = options.Covariance;
        if strcmp(options.Target,'standard')
            Sigma=eye(d);
        elseif isempty(Sigma)
            scales = linspace(0.75,1.5,d);
            Sigma = (scales'*scales).*toeplitz(options.Rho.^(0:d-1));
        end
        step = options.StepScale;
        if isempty(step), step = 2.38/sqrt(d); end
        validateattributes(step,{'double'},{'scalar','positive','finite'});
        estimates = nan(R,numel(Nvec),nFunctions,nMethods);
        acceptance = zeros(R,numel(Nvec),nMethods);
        collision = zeros(R,numel(Nvec),nMethods);
        seeds = zeros(R,numel(Nvec),nMethods);
        for ni = 1:numel(Nvec)
            nTimer=tic;
            % Preserve original method IDs (1,3) and dimension/grid seed namespaces.
            [repIDs,methodIndex]=ndgrid(1:R,1:nMethods);
            repIDs=repIDs(:);methodIndex=methodIndex(:);
            taskIDs=reshape(methodIDs(methodIndex),[],1);
            taskSeeds=options.Seed+(((di-1)*numel(Nvec)+ni-1)*3+taskIDs-1)*R+repIDs;
            assert(all(taskSeeds<2^32),'Seed range exceeds mt19937ar limit.');
            batch=cell(numel(repIDs),1);N=Nvec(ni);
            parfor (k=1:numel(repIDs),workerCount)
                batch{k}=run_one(N,mu,Sigma,step,methods{methodIndex(k)},taskSeeds(k),options);
            end
            for k=1:numel(repIDs)
                r=repIDs(k);mi=methodIndex(k);out=batch{k};
                seeds(r,ni,mi)=taskSeeds(k);
                estimates(r,ni,:,mi,1)=reshape(out.average,1,1,nFunctions);
                acceptance(r,ni,mi)=out.acceptance;
                collision(r,ni,mi)=out.collisionMax;
            end
            % This runner has no resume API; retain the completed N values for recovery.
            results.progress=struct('dimension',d,'N',Nvec,'completedN',Nvec(1:ni), ...
                'estimates',estimates,'acceptance',acceptance,'seeds',seeds, ...
                'maxCellCollisions',collision,'testNames',{{'projection'}});
            temporary=[stem '_checkpoint_pending.mat'];save(temporary,'results','-v7.3');
            movefile(temporary,[stem '_checkpoint.mat'],'f');
            fprintf('N=2^%d completed in %.2f s\n',log2(Nvec(ni)),toc(nTimer));
        end
        fprintf('Statistical block total: %.2f s\n',toc(caseTimer));
        names = out.testNames; reference = out.reference;
        caseResult = struct('dimension',d,'mu',mu,'Sigma',Sigma,'stepScale',step, ...
            'N',Nvec,'testNames',{names},'reference',reference, ...
            'estimates',estimates,'acceptance',acceptance,'maxCellCollisions',collision, ...
            'seeds',seeds,'hilbertBits',bits,'timing',[]);
        if options.RunTiming
            fprintf('\nTiming: dimension=%d; test=projection\n',d);
            timingN=select_timing_N(Nvec,options.TimingNExponents);
            callback=@(n,method,seed)run_one(n,mu,Sigma,step,method,seed,options);
            caseResult.timing=benchmark_seeded_kernel(callback,timingN,methods, ...
                options.Seed,di,seeds,'Replications',options.TimingReplications, ...
                'WarmupExponent',options.TimingWarmupExponent);
            caseResult.timing.testName='projection';
            caseResult.timing.T=options.T;caseResult.timing.T0=options.BurnIn;
        end
        serialSeconds=standalone_timing_seconds(caseResult,Nvec,methods);
        paired=reshape(permute(estimates,[1 2 3 5 4]),R,numel(Nvec),[],2);
        bs=bootstrap_macro_statistics(paired,Nvec,'Bootstrap',options.Bootstrap, ...
            'Windows',options.TailCount,'Seed',19002026);
        caseResult.bootstrap=bs;
        results.cases = [results.cases; caseResult];
        results=rmfield(results,'progress');

        %% Compute empirical statistics and bootstrap slope intervals
        estimatorNames = {'time_average'};
        for ei = 1:numel(estimatorNames)
            for fi = 1:nFunctions
                iidVariance = var(estimates(:,:,fi,1,ei),0,1);
                for mi = 1:nMethods
                    A = reshape(estimates(:,:,fi,mi,ei),R,[]);
                    variance = var(A,0,1); mse = mean((A-reference(fi)).^2,1);
                    fitCount = min(options.TailCount,numel(Nvec));
                    tail = numel(Nvec)-fitCount+1:numel(Nvec);
                    g=(ei-1)*nFunctions+fi;
                    slope=bs.slope(g,mi);ci=[bs.slopeLow(g,mi),bs.slopeHigh(g,mi)];
                    slopeRows(end+1,:) = {d,d+2,estimatorNames{ei},names{fi},methods{mi}, ...
                        slope,ci(1),ci(2),-1-1/(d+2),Nvec(tail(1)),Nvec(tail(end)),bs.slopeValid(g,mi),bs.slopeStatus(g,mi)}; %#ok<AGROW>
                    for ni = 1:numel(Nvec)
                        low=NaN;high=NaN;valid=0;status="not_max_N";
                        if ni==numel(Nvec) && mi==2
                            low=bs.vrfLow(g);high=bs.vrfHigh(g);valid=bs.vrfValid(g);status=bs.vrfStatus(g);
                        elseif mi==1,status="baseline";end
                        rows(end+1,:) = {d,Nvec(ni),estimatorNames{ei},names{fi},methods{mi}, ...
                            mean(A(:,ni)),reference(fi),mean(A(:,ni))-reference(fi),variance(ni), ...
                            mse(ni),sqrt(mse(ni)),iidVariance(ni)/variance(ni), ...
                            mean(acceptance(:,ni,mi)),serialSeconds(ni,mi), ...
                            max(collision(:,ni,mi)), ...
                            iidVariance(ni)/variance(ni)*serialSeconds(ni,1)/serialSeconds(ni,mi),low,high,valid,status}; %#ok<AGROW>
                    end
                end
            end
        end
    end

    %% Export raw results and derived tables
    results.summary = cell2table(rows,'VariableNames',{'d','N','estimator','test','method', ...
        'estimate_mean','reference','bias','variance','mse','rmse','vrf','acceptance','seconds','max_cell_collisions','tae','vrf_ci_low','vrf_ci_high','vrf_bootstrap_valid','vrf_bootstrap_status'});
    results.slopes = cell2table(slopeRows,'VariableNames',{'d','point_dimension','estimator','test','method', ...
        'variance_slope','bootstrap_low','bootstrap_high', ...
        'point_dimension_reference','fit_N_min','fit_N_max','bootstrap_valid','bootstrap_status'});
    save([stem '.mat'],'results','-v7.3');
    writetable(results.summary,[stem '_summary.csv']);
    writetable(results.slopes,[stem '_slopes.csv']);
    fprintf('Saved %s.mat\n',stem);
    fprintf('Total elapsed time (this invocation): %.2f s\n',toc(runTimer));
end

%% Local helper: initialize_parallel
% Resolve serial/process execution; reject thread pools that share global streams.
function [enabled,workers] = initialize_parallel(options)
    enabled = false; workers = 0;
    if ~options.UseParallel
        return
    end
    if ~license('test','Distrib_Computing_Toolbox') || isempty(ver('parallel'))
        warning('gaussian_hilbert:ParallelUnavailable', ...
            'Parallel Computing Toolbox unavailable; using serial execution.');
        return
    end
    pool = gcp('nocreate');
    if ~isempty(pool) && isa(pool,'parallel.ThreadPool')
        error('gaussian_hilbert:ThreadPoolUnsupported', ...
            ['A thread pool is active. Use a process pool for independent global ', ...
             'random streams, or set UseParallel=false.']);
    end
    if isempty(pool), pool = parpool('Processes',options.NumWorkers); end
    enabled = true; workers = min(pool.NumWorkers,options.NumWorkers);

end

%% Local helper: run_one
% Forward resolved experiment settings to one replication kernel.
function out = run_one(N,mu,Sigma,step,method,seed,options)
    % Resolve zero initialization before dispatching to the replication kernel.
    % This also works with kernels supporting only stationary/fixed.
    initialization=options.Initialization;initialState=options.InitialState;
    if strcmp(initialization,'zero')
        initialization='fixed';initialState=zeros(size(mu));
    end
    out = gaussian_rwm_replication(N,options.BurnIn,options.T,mu,Sigma, ...
        step,method,seed,'Initialization',initialization, ...
        'InitialState',initialState,'SortCoordinates',options.SortCoordinates, ...
        'HilbertBits',options.HilbertBits);
end
