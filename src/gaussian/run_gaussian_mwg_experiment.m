function results=run_gaussian_mwg_experiment(options)
%RUN_GAUSSIAN_MWG_EXPERIMENT Product Gaussian variance study with one combined parfor per N.
% T counts complete sweeps. Saves time averages only, never particle paths.
% Resume checks statistical settings; worker count and analysis can change.
%
% Raw case estimates have axes replication, N, test and method. Completion
% masks identify pending replications; random substreams depend on case indices,
% not on worker scheduling. Full defaults can be expensive.

    %% Declare experiment, resume and analysis options
    arguments
        options.Dimensions (1,:) double {mustBeInteger,mustBePositive} = [2 4 8]
        options.NExponents (1,:) double {mustBeInteger,mustBePositive} = 10:18
        options.Replications (1,1) double {mustBeInteger,mustBePositive} = 200
        options.T (1,1) double {mustBeInteger,mustBePositive} = 128
        options.BurnIn (1,1) double {mustBeInteger,mustBeNonnegative} = 8
        options.StepScale (1,1) double {mustBePositive,mustBeFinite} = 2.38
        options.Initialization (1,:) char {mustBeMember(options.Initialization,{'zero','fixed','stationary'})} = 'zero'
        options.InitialState double = []
        options.Seed (1,1) double {mustBeInteger,mustBeNonnegative} = 2026
        options.TailCount (1,1) double {mustBeInteger,mustBePositive} = 5
        options.Bootstrap (1,1) double {mustBeInteger,mustBeNonnegative} = 1000
        options.RunTiming (1,1) logical = true
        options.TimingNExponents (1,:) double {mustBeInteger,mustBePositive} = []
        options.TimingReplications (1,1) double {mustBeInteger,mustBePositive} = 20
        options.TimingWarmupExponent (1,1) double {mustBeInteger,mustBePositive} = 12
        options.UseParallel (1,1) logical = true
        options.NumWorkers (1,1) double {mustBeInteger,mustBePositive} = 12
        options.ValidateStratification (1,1) logical = false
        options.OutputDir (1,:) char = ''
        options.ResumePath (1,:) char = ''
    end

    if options.RunTiming
        assert(options.TimingReplications<=499,'Use at most 499 independent timing samples');
        select_timing_N(2.^options.NExponents,options.TimingNExponents);
    end

    %% Validate dimensions, particle grid and initial state
    assert(~isempty(options.Dimensions) && all(options.Dimensions<=8) && ...
        numel(unique(options.Dimensions))==numel(options.Dimensions),'Use distinct dimensions from 1:8');
    assert(numel(options.NExponents)>=2 && all(diff(options.NExponents)>0) && ...
        max(options.NExponents)<=26,'Use increasing exponents between 1 and 26');
    assert(options.T>options.BurnIn && options.Replications>=3 && options.Replications<=100000);
    runTimer=tic;
    assert(options.Seed<2^32 && options.TailCount>=2);
    assert(isreal(options.InitialState) && all(isfinite(options.InitialState),'all'));
    if ~strcmp(options.Initialization,'fixed')
        assert(isempty(options.InitialState),'InitialState is only used with fixed initialization');
    elseif ~isempty(options.InitialState) && ~isscalar(options.InitialState)
        assert(isscalar(options.Dimensions) && isequal(size(options.InitialState),[1 options.Dimensions]), ...
            'For a vector InitialState, specify one matching dimension');
    end

    %% Resolve output folder and initialize or resume compatible cases
    root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
    folder=options.OutputDir;
    if isempty(folder),folder=fullfile(root,'results','gaussian_mwg');end
    if ~isfolder(folder),mkdir(folder);end
    [ok,attr]=fileattrib(folder);assert(ok);folder=attr.Name;
    resume=options.ResumePath;
    if strcmp(resume,'auto')
        candidates=dir(fullfile(folder,'gaussian_mwg_v2_*_checkpoint.mat'));
        resume='';
        if ~isempty(candidates)
            [~,last]=max([candidates.datenum]);
            resume=fullfile(candidates(last).folder,candidates(last).name);
        end
    end
    if isempty(resume)
        tag=char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'));
        stem=fullfile(folder,['gaussian_mwg_v2_' tag]);
        assert(~isfile([stem '.mat']) && ~isfile([stem '_checkpoint.mat']));
        results=struct('schema','gaussian-mwg-v2','settings',options, ...
            'methods',{{'iid_mwg','array_mwg'}},'cases',[], ...
            'target','N_p(0,I_p)','outputStem',stem,'matlabVersion',version);
        results.algorithm=struct('scanOrder','1:p','timeUnit','complete sweep', ...
            'pointDimension',3,'scramble','MatousekAffineOwen', ...
            'freshEachCoordinateUpdate',true,'generator','mrg32k3a', ...
            'substreamFormula','1+8*(100000*(2*(32*(p-1)+log2(N))+methodID-1)+replication-1)+coordinate-1');
        results.testDescription={ ...
            'projection: sum(x_r)/sqrt(p); unbounded, polynomial extension has log(N) factor'};
        R=options.Replications;N=2.^options.NExponents;
        for p=options.Dimensions
            c=struct('dimension',p,'mu',zeros(1,p),'Sigma',eye(p), ...
                'stepScale',options.StepScale,'N',N, ...
                'testNames',{{'projection'}}, ...
                'stationaryReference',0, ...
                'estimates',nan(R,numel(N),1,2),'timing',[], ...
                'acceptance',nan(R,numel(N),2), ...
                'acceptanceByCoordinate',nan(R,numel(N),p,2), ...
                'substreamStarts',nan(R,numel(N),2), ...
                'completed',false(R,numel(N),2));
            results.cases=[results.cases;c];
        end
    else
        saved=load(resume,'results');
        assert(isfield(saved,'results') && isfield(saved.results,'schema') && ...
            strcmp(saved.results.schema,'gaussian-mwg-v2'),'Expected a v2 single-test checkpoint; start a new run for legacy results');
        results=saved.results;
        % Drop obsolete runtime diagnostics when explicitly resuming old checkpoints.
        obsolete=intersect(fieldnames(results.cases),{'statisticalSeconds','batchWallSeconds'});
        if ~isempty(obsolete),results.cases=rmfield(results.cases,obsolete);end
        keys={'Dimensions','NExponents','Replications','T','BurnIn','StepScale','Initialization','InitialState','Seed', ...
            'RunTiming','TimingNExponents','TimingReplications','TimingWarmupExponent'};
        for j=1:numel(keys)
            key=keys{j};
            assert(isequaln(results.settings.(key),options.(key)), ...
                'Resume setting %s differs. Restore the saved setting, or set ResumePath to empty for a new run.',key);
        end
        [~,name]=fileparts(resume);name=regexprep(name,'_checkpoint$','');
        stem=fullfile(folder,name);results.outputStem=stem;
        fprintf('Resuming %s\n',resume);
    end

    %% Configure workers only when simulation work remains
    pending=any(arrayfun(@(c)any(~c.completed,'all'),results.cases));
    workers=0;
    if options.UseParallel && pending
        if ~license('test','Distrib_Computing_Toolbox') || isempty(ver('parallel'))
            warning('gaussian_mwg:NoParallel','Parallel toolbox unavailable; executing serially');
        else
            pool=gcp('nocreate');
            if isempty(pool),pool=parpool('Processes',options.NumWorkers);end
            assert(~isa(pool,'parallel.ThreadPool'),'Use a process pool, or set UseParallel=false');
            workers=min(pool.NumWorkers,options.NumWorkers);
        end
    end
    results.execution=struct('workerLimit',workers,'scheduler','single-parfor-per-N', ...
        'validateStratification',options.ValidateStratification);
    %% Run pending replications and persist completed replications
    for ci=1:numel(results.cases)
        c=results.cases(ci);p=c.dimension;caseTimer=tic;
        fprintf('\n[%d/%d] dimension=%d; tests=projection\n',ci,numel(results.cases),p);
        print_experiment_header('Gaussian MWG statistical experiment',c.testNames, ...
            log2(c.N),options.BurnIn,options.T,options.Replications,workers>0);
        for ni=1:numel(c.N)
            N=c.N(ni);e=log2(N);
            if all(c.completed(:,ni,:),'all')
                fprintf('N=2^%d already completed; skipped.\n',e);continue;
            end
            nTimer=tic;
            [reps,methodIDs]=find(~reshape(c.completed(:,ni,:),size(c.completed,1),2));
            starts=1+8*(100000*(2*(32*(p-1)+e)+methodIDs-1)+reps-1);
            methods=results.methods;batch=cell(numel(reps),1);
            parfor (j=1:numel(reps),workers)
                batch{j}=run_one(N,p,methods{methodIDs(j)},starts(j),options);
            end
            for j=1:numel(reps)
                r=reps(j);mi=methodIDs(j);out=batch{j};
                assert(all(isfinite(out.average)),'Nonfinite MWG estimates');
                c.estimates(r,ni,:,mi)=reshape(out.average,1,1,1);
                c.acceptance(r,ni,mi)=out.acceptance;
                c.acceptanceByCoordinate(r,ni,:,mi)=reshape(out.acceptanceByCoordinate,1,1,p);
                c.substreamStarts(r,ni,mi)=starts(j);c.completed(r,ni,mi)=true;
            end
            results.cases(ci)=c;save_checkpoint(stem,results);
            fprintf('N=2^%d completed in %.2f s\n',e,toc(nTimer));
        end
        fprintf('Statistical block total: %.2f s\n',toc(caseTimer));
    end

    %% Independently time only the projection after all statistical work finishes
    for ci=1:numel(results.cases)
        c=results.cases(ci);p=c.dimension;
        if options.RunTiming && isempty(c.timing)
            fprintf('\nTiming: dimension=%d; test=projection\n',p);
            timingN=select_timing_N(c.N,options.TimingNExponents);
            callback=@(n,method,seed)gaussian_mwg_replication(n,options.BurnIn, ...
                options.T,p,options.StepScale,method,seed,1, ...
                'Initialization',options.Initialization,'InitialState',options.InitialState, ...
                'ValidateStratification',options.ValidateStratification);
            c.timing=benchmark_seeded_kernel(callback,timingN,results.methods, ...
                options.Seed,p,options.Seed,'Replications',options.TimingReplications, ...
                'WarmupExponent',options.TimingWarmupExponent);
            c.timing.testName='projection';c.timing.coordinateSubstreams=1:p;
            c.timing.T=options.T;c.timing.T0=options.BurnIn;
            results.cases(ci)=c;save_checkpoint(stem,results);
        end
    end

    %% Summarize matching kernels and export final run artifacts
    [results.summary,results.slopes]=summarize_gaussian_mwg(results, ...
        'TailCount',options.TailCount,'Bootstrap',options.Bootstrap);
    results.analysisSettings=struct('TailCount',options.TailCount,'Bootstrap',options.Bootstrap, ...
        'protocol','macro-percentile-v1','Seed',19002026,'confidence',0.95, ...
        'pairedN',false,'pairedMethods',false,'pairedOutputs',true,'minimumValidFraction',0.9);
    save([stem '.mat'],'results','-v7.3');
    writetable(results.summary,[stem '_summary.csv']);
    writetable(results.slopes,[stem '_slopes.csv']);
    fprintf('Finished. Raw results: %s.mat\n',stem);
    fprintf('Total elapsed time (this invocation): %.2f s\n',toc(runTimer));
end

%% Local helper: save_checkpoint
% Save to a temporary MAT before replacing this run's checkpoint.
function save_checkpoint(stem,results)
    temporary=[stem '_checkpoint_pending.mat'];
    save(temporary,'results','-v7.3');
    movefile(temporary,[stem '_checkpoint.mat'],'f');
end

%% Local helper: run_one
% Forward resolved experiment settings to one replication kernel.
function out=run_one(N,p,method,start,options)
    out=gaussian_mwg_replication(N,options.BurnIn,options.T,p,options.StepScale, ...
        method,options.Seed,start,'Initialization',options.Initialization, ...
        'InitialState',options.InitialState,'ValidateStratification',options.ValidateStratification);
end
