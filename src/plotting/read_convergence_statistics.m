function output=read_convergence_statistics(matPath,options)
% Read complete selected time-average estimates; compute variance and VRF.
% No figures or output files are created.
%
% Input supports joint Gaussian results and Bayesian result schemas.
% Only estimator index 1 (time average) is read. Output.summary has one row
% per selected dimension, particle count, test and method.

    %% Declare saved-file selection
    arguments
        matPath (1,:) char
        options.NExponents (1,:) double {mustBeInteger,mustBePositive} = 10:18
        options.Bootstrap (1,1) double {mustBeInteger,mustBeNonnegative} = 0
        options.Dimensions (1,:) double = []
    end

    %% Load and normalize supported MAT schemas
    assert(isfile(matPath),'Result file not found: %s',matPath);
    assert(~isempty(options.NExponents) && all(diff(options.NExponents)>0));
    [ok,info]=fileattrib(matPath);assert(ok);matPath=info.Name;
    saved=load(matPath);requestedN=2.^options.NExponents;
    if isfield(saved,'results')
        r=saved.results;cases=r.cases;methods=r.methods;dataset='gaussian';
    elseif isfield(saved,'result') && isfield(saved.result,'model')
        r=saved.result;methods=r.methods;dataset=r.model.dataset;
        names=[{'smooth','indicator'},arrayfun(@(j)sprintf('beta_%d',j),1:r.model.p,'UniformOutput',false),{'predictive'}];
        if isfield(r,'testNames'),names=r.testNames;end
        cases=struct('dimension',r.model.p,'N',r.N,'testNames',{names},'estimates',r.estimates);
        if isfield(r,'timing'),cases.timing=r.timing;end
    else
        error('Expected a Gaussian results or Bayesian result MAT, not a sensitivity-study file');
    end

    %% Select dimensions and verify completeness
    if ~isempty(options.Dimensions)
        assert(all(ismember(options.Dimensions,[cases.dimension])),'Requested dimension missing');
        cases=cases(ismember([cases.dimension],options.Dimensions));
    end
    ii=find(strcmp(methods,'iid'));hi=find(strcmp(methods,'hilbert'));
    assert(isscalar(ii) && isscalar(hi),'Expected exactly one IID and one Hilbert method');
    % Validate every selected case before writing anything.
    for ci=1:numel(cases)
        c=cases(ci);[present,idx]=ismember(requestedN,c.N);
        assert(all(present),'Requested N missing from d=%d: %s',c.dimension,mat2str(requestedN(~present)));
        if isfield(r,'completed')
            for ni=idx
                for mi=[ii hi]
                    assert(all(r.completed(:,ni,mi)),...
                        'N=%d %s has only %d/%d saved repetitions; requested range is incomplete',...
                        c.N(ni),methods{mi},nnz(r.completed(:,ni,mi)),size(r.completed,1));
                end
            end
        end
        A=c.estimates(:,idx,:,[ii hi],1);
        assert(all(isfinite(A),'all'),'Selected d=%d estimates contain missing/nonfinite values',c.dimension);
        assert(size(A,1)>=2,'At least two repetitions are needed');
    end

    %% Recompute variances, matching IID ratios and standalone median times
    rows={};estimatorNames={'time_average'};
    for ci=1:numel(cases)
        c=cases(ci);[~,idx]=ismember(requestedN,c.N);R=size(c.estimates,1);
        medianSeconds=standalone_timing_seconds(c,requestedN,methods([ii hi]));
        tests=c.testNames;components=num2cell(1:numel(tests));
        if ~strcmp(dataset,'gaussian') && ~isfield(r,'testNames')
            tests{end+1}='coefficient_total';components{end+1}=3:c.dimension+2;
        end
        fullA=c.estimates(:,:,:,[ii hi],1);
        bs=bootstrap_macro_statistics(fullA,c.N,'Groups',components,'Bootstrap',options.Bootstrap);
        V=zeros(numel(requestedN),numel(tests),2);
        for ei=1
            for fi=1:numel(tests)
                testSeconds=medianSeconds;
                if ~strcmp(dataset,'gaussian') && isfield(c,'timing') && ~isempty(c.timing) && ~strcmp(c.timing.testName,tests{fi})
                    testSeconds(:)=NaN;
                end
                for mi=1:2
                    sourceMethod=[ii hi];
                    A=reshape(c.estimates(:,idx,components{fi},sourceMethod(mi),ei),R,numel(idx),[]);
                    % For coefficient_total, this is a covariance trace.
                    V(:,fi,mi,ei)=reshape(sum(var(A,0,1),3),[],1);
                end
                for mi=1:2
                    for ni=1:numel(requestedN)
                        low=NaN;high=NaN;valid=0;status="not_max_N";
                        if requestedN(ni)==max(c.N) && mi==2
                            low=bs.vrfLow(fi);high=bs.vrfHigh(fi);valid=bs.vrfValid(fi);status=bs.vrfStatus(fi);
                        elseif mi==1,status="baseline";end
                        rows(end+1,:)={dataset,c.dimension,requestedN(ni),R,estimatorNames{ei},...
                            tests{fi},methods{sourceMethod(mi)},V(ni,fi,mi,ei),...
                            V(ni,fi,1,ei)/V(ni,fi,mi,ei),testSeconds(ni,mi),low,high,valid,status}; %#ok<AGROW>
                    end
                end
            end
        end
    end

    %% Return a common table with source provenance
    summary=cell2table(rows,'VariableNames',{'dataset','d','N','R','estimator','test','method','variance','vrf','seconds','vrf_ci_low','vrf_ci_high','vrf_bootstrap_valid','vrf_bootstrap_status'});
    output=struct('source',matPath,'N',requestedN,'summary',summary);
end
