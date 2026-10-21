function result = run_bayesian_single_parfor(options)
% Bayesian runner: one combined IID/Hilbert parfor per N.
% Checkpoints are saved only after all pending tasks at an N complete.
% An interruption loses unsaved work at the current N; serial timing is unchanged.
%
% Returns result and writes raw MAT, checkpoints, CSV and JSON summaries.
% Fresh runs reconstruct and verify input data before simulating. On resume,
% result stores the simulation settings while runtime controls govern workers.

    %% Declare run and runtime options
    arguments
        options.Dataset (1,:) char {mustBeMember(options.Dataset,{'ripley','pima'})} = 'ripley'
        options.NExponents (1,:) double {mustBeInteger,mustBePositive} = 8:16
        options.Replications (1,1) double {mustBeInteger,mustBePositive} = 200
        options.T (1,1) double {mustBeInteger,mustBePositive} = 128
        options.BurnIn (1,1) double {mustBeInteger,mustBeNonnegative} = 8
        options.HilbertBits double = []
        options.PrecisionCheck (1,1) logical = false
        options.CollisionCheck (1,1) logical = false
        options.StatisticalTestFunction (1,:) char = 'projection' % Paper statistical test
        options.TimingTestFunction (1,:) char = 'projection'
        options.RunTiming (1,1) logical = true
        options.TimingNExponents (1,:) double {mustBeInteger,mustBePositive} = []
        options.TimingReplications (1,1) double {mustBeInteger,mustBePositive} = 20
        options.TimingWarmupExponent (1,1) double {mustBeInteger,mustBePositive} = 12
        options.UseParallel (1,1) logical = true
        options.NumWorkers (1,1) double {mustBeInteger,mustBePositive} = 4
        options.Seed (1,1) double {mustBeInteger,mustBeNonnegative} = 2026
        options.TailCount (1,1) double {mustBeInteger,mustBePositive} = 5
        options.Bootstrap (1,1) double {mustBeInteger,mustBeNonnegative} = 1000
        options.OutputDir (1,:) char = ''
        options.ResumePath (1,:) char = ''
    end

    runTimer=tic;
    if options.RunTiming
        assert(options.TimingReplications<=499,'Use at most 499 independent timing samples');
        select_timing_N(2.^options.NExponents,options.TimingNExponents);
    end

    %% Validate the grid and seed range
    assert(options.Replications>=3 && options.Replications<=999 && options.T>options.BurnIn);
    assert(numel(options.NExponents)>=3 && all(diff(options.NExponents)>0));
    assert(options.TailCount>=3);
    assert(options.Seed+2000000+max(options.NExponents)*10000+2000+options.Replications<2^32);

    %% Resolve output paths and initialize or resume a run
    root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
    addpath(genpath(fullfile(root,'src')));
    output=options.OutputDir; if isempty(output),output=fullfile(root,'results','bayesian');end
    if ~isfolder(output),mkdir(output);end
    if isempty(options.ResumePath)
        dataset=options.Dataset;
        model=load_bayesian_model(dataset); preprocess_bayesian_data(dataset);
        [~,testNames]=bayesian_functions(zeros(0,model.p),model,options.StatisticalTestFunction);
        bayesian_functions(zeros(0,model.p),model,options.TimingTestFunction);
        assert(~strcmp(options.TimingTestFunction,'all'),'Timing requires exactly one test');
        bits=resolve_hilbert_bits(model.p,options.HilbertBits);
        R=options.Replications; Nvec=2.^options.NExponents; methods={'iid','hilbert'};
        tag=char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'));
        stem=fullfile(output,[dataset '_single_parfor_v3_' tag]);
        assert(~isfile([stem '.mat']) && ~isfile([stem '_checkpoint.mat']));
        result=struct('schema','bayesian-tests-v3','testNames',{testNames},'timing',[],'phase','formal','model',model,'N',Nvec,'R',R,'T',options.T,'T0',options.BurnIn,...
            'methods',{methods},'estimates',nan(R,numel(Nvec),numel(testNames),2),...
            'acceptance',nan(R,numel(Nvec),2),...
            'seeds',nan(R,numel(Nvec),2),'collisions',nan(R,numel(Nvec),2),...
            'precisionChangedFraction',nan(R,numel(Nvec),2,2),...
            'completed',false(R,numel(Nvec),2),'scheduler','single-parfor-per-N',...
            'seedFormula','settings.Seed + datasetID*1000000 + log2(N)*10000 + methodID*1000 + replication',...
            'HilbertBits',bits,'matlabVersion',version,'outputStem',stem,'settings',options);
        result.dataManifest=jsondecode(fileread(fullfile(root,'data','source_manifest.json')));
        [~,hash]=system(['git -C "' root '" rev-parse HEAD']); result.gitCommit=strtrim(hash);
        result.sobol=struct('dimension',model.p+2,'Skip',0,'Leap',0,'PointOrder','graycode',...
            'scramble','MatousekAffineOwen','freshEachStep',true);
    else
        loaded=load(options.ResumePath,'result'); result=loaded.result;
        assert(isfield(result,'scheduler') && strcmp(result.scheduler,'single-parfor-per-N'), ...
            'Resume only this runner checkpoints; existing batched results are not modified');
        assert(strcmp(result.phase,'formal') && strcmp(result.model.dataset,options.Dataset),...
            'Resume checkpoint must match the requested dataset');
        assert(isfield(result,'schema') && strcmp(result.schema,'bayesian-tests-v3'), ...
            'A v3 multi-test checkpoint is required; start a new run');
        assert(size(result.estimates,5)==1,'Checkpoint must contain only time averages; start a new run');
        assert(isfield(result.settings,'StatisticalTestFunction'), ...
            'Bayesian:MissingStatisticalTestFunction', ...
            'Checkpoint lacks StatisticalTestFunction; start a new run');
        assert(strcmp(result.settings.StatisticalTestFunction,options.StatisticalTestFunction), ...
            'Bayesian:StatisticalTestFunctionMismatch', ...
            'StatisticalTestFunction differs from checkpoint; start a new run');
        keys={'NExponents','Replications','T','BurnIn','HilbertBits','Seed', ...
            'StatisticalTestFunction','TimingTestFunction','RunTiming','TimingNExponents','TimingReplications','TimingWarmupExponent'};
        for k=1:numel(keys)
            key=keys{k};
            assert(isequaln(result.settings.(key),options.(key)), ...
                'Resume setting %s differs; restore saved settings or start a new run',key);
        end
        current=load_bayesian_model(options.Dataset);
        assert(isequal(current.X,result.model.X) && isequal(current.y,result.model.y),...
            'Dataset changed; do not resume');
        [~,resumeName]=fileparts(options.ResumePath);
        resumeName=regexprep(resumeName,'_checkpoint$','');
        stem=fullfile(output,resumeName);result.outputStem=stem;
    end

    %% Configure process workers and statistical diagnostic metadata
    workerLimit=0;
    if options.UseParallel
        pool=gcp('nocreate'); if isempty(pool),pool=parpool('Processes',options.NumWorkers);end
        assert(~isa(pool,'parallel.ThreadPool'),'Use a process pool');
        workerLimit=min(pool.NumWorkers,options.NumWorkers);
    end
    result.workerCount=workerLimit;
    if ~isfield(result,'implementationProfile')
        result.implementationProfile=repmat({''},size(result.completed));
        result.implementationProfile(result.completed)={'legacy-bitstring'};
    end
    if ~isfield(result,'precisionComparisons')
        result.precisionComparisons=zeros(size(result.completed));
        if ~isfield(result.settings,'PrecisionCheck') && result.HilbertBits==20
            for legacyN=find(result.N==65536)
                result.precisionComparisons(:,legacyN,2)=double(result.completed(:,legacyN,2))*result.N(legacyN)*result.T;
            end
        end
    end
    diagnostics=struct('PrecisionCheck',options.PrecisionCheck,'CollisionCheck',options.CollisionCheck);
    result.runtimeDiagnostics=diagnostics;
    result.precisionReferenceBits=[20 28]; % Columns of precisionChangedFraction.
    statisticalTestFunction=result.settings.StatisticalTestFunction;
    seedBase=result.settings.Seed; % Fresh or resumed v3 results retain their own seed.
    logs=fullfile(output,'logs');if ~isfolder(logs),mkdir(logs);end
    [~,name]=fileparts(stem);diary(fullfile(logs,[name '.log'])); cleanup=onCleanup(@()diary('off'));
    fprintf('\nDataset=%s\n',result.model.dataset);
    print_experiment_header('Bayesian statistical experiment',result.testNames, ...
        log2(result.N),result.T0,result.T,result.R,workerLimit>0);
    statsTimer=tic;

    %% Combine both methods in one task list, preserving original seed identities
    for ni=1:numel(result.N)
        if all(result.completed(:,ni,:),'all')
            fprintf('N=2^%d already completed; skipped.\n',log2(result.N(ni)));continue;
        end
        nTimer=tic;
        % Rows are replication IDs; columns are method IDs, independent of scheduling.
        pending=~reshape(result.completed(:,ni,:),result.R,2);
        [replicationIDs,methodIDs]=find(pending);
        datasetID=1+strcmp(result.model.dataset,'pima');
        seeds=seedBase+datasetID*1000000+log2(result.N(ni))*10000+methodIDs*1000+replicationIDs;
        outputs=cell(numel(replicationIDs),1);
        model=result.model;N=result.N(ni);T=result.T;T0=result.T0;
        methods=result.methods;bits=result.HilbertBits;
        parfor (k=1:numel(replicationIDs),workerLimit)
            outputs{k}=bayesian_rwm_replication(N,T0,T,model, ...
                methods{methodIDs(k)},seeds(k),bits,false,diagnostics,statisticalTestFunction);
        end
        %% Scatter by original replication/method IDs, never by completion order
        for k=1:numel(replicationIDs)
            r=replicationIDs(k);mi=methodIDs(k);out=outputs{k};
            result.estimates(r,ni,:,mi,1)=reshape(out.average,1,1,[]);
            result.acceptance(r,ni,mi)=out.acceptance;
            result.seeds(r,ni,mi)=seeds(k);
            result.collisions(r,ni,mi)=out.maxCollisions;
            result.precisionChangedFraction(r,ni,mi,:)=reshape(out.precisionChangedFraction,1,1,1,2);
            result.precisionComparisons(r,ni,mi)=out.precisionComparisons;
            result.implementationProfile{r,ni,mi}='cached-hybrid-key-v3';
            result.completed(r,ni,mi)=true;
        end
        % One checkpoint per completed N; use the run-specific temporary path.
        temp=[stem '_checkpoint_pending.mat'];save(temp,'result','-v7.3');
        movefile(temp,[stem '_checkpoint.mat'],'f');
        fprintf('N=2^%d completed in %.2f s\n',log2(result.N(ni)),toc(nTimer));
    end
    fprintf('Statistical block total: %.2f s\n',toc(statsTimer));

    %% Independently time only the selected test, with diagnostics disabled
    if options.RunTiming && isempty(result.timing)
        fprintf('\nTiming: dataset=%s; test=%s\n',result.model.dataset,options.TimingTestFunction);
        timingN=select_timing_N(result.N,options.TimingNExponents);
        model=result.model;bits=result.HilbertBits;
        callback=@(n,method,seed)bayesian_rwm_replication(n,result.T0,result.T, ...
            model,method,seed,bits,false,struct('PrecisionCheck',false,'CollisionCheck',false),options.TimingTestFunction);
        result.timing=benchmark_seeded_kernel(callback,timingN,result.methods, ...
            result.settings.Seed,1+strcmp(model.dataset,'pima'),result.seeds, ...
            'Replications',options.TimingReplications,'WarmupExponent',options.TimingWarmupExponent);
        result.timing.testName=options.TimingTestFunction;result.timing.T=result.T;result.timing.T0=result.T0;
    end

    %% Verify completion and export derived statistics
    assert(all(result.completed,'all') && numel(unique(result.seeds(:)))==numel(result.seeds));
    save([stem '.mat'],'result','-v7.3');
    save([stem '_checkpoint.mat'],'result','-v7.3');
    if isfield(result,'settings')
        export_bayesian_statistics([stem '.mat'],result.settings.TailCount,result.settings.Bootstrap);
    else
        export_bayesian_statistics([stem '.mat']);
    end
    small=rmfield(result,'estimates');fid=fopen([stem '.json'],'w');
    assert(fid>=0);fprintf(fid,'%s',jsonencode(small,'PrettyPrint',true));fclose(fid);
    fprintf('Finished. Raw results: %s.mat\n',stem);
    fprintf('Total elapsed time (this invocation): %.2f s\n',toc(runTimer));
end
