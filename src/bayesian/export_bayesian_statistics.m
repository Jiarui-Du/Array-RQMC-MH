function [summary,slopes] = export_bayesian_statistics(path,tailCount,bootstrapCount,outputStem)
% Regenerate variance outputs strictly from raw per-replication estimates.
%
% Input MAT must contain result.estimates with axes replication, N, test,
% method; legacy files may have an additional estimator axis. No chains are rerun.
% Legacy coefficient_total is the sum of coefficient variances. Paper drivers use projection.

    %% Load raw estimates and resolve analysis defaults
    if nargin<2,tailCount=5;end
    if nargin<3,bootstrapCount=1000;end
    loaded=load(path,'result'); result=loaded.result; R=result.R; N=result.N;
    names=[{'smooth','indicator'},arrayfun(@(k)sprintf('beta_%d',k),1:result.model.p,'UniformOutput',false),{'predictive'}];
    legacy=~isfield(result,'testNames');
    if ~legacy,names=result.testNames;end
    serialSeconds=standalone_timing_seconds(result,N,result.methods);
    rows={}; fits={}; estimator={'time_average'};

    groups=num2cell(1:numel(names));
    if legacy,groups{end+1}=3:result.model.p+2;end
    baseGroups=groups;
    % Read only retained averages, including when exporting historical files.
    paired=reshape(permute(result.estimates(:,:,:,:,1),[1 2 3 5 4]),R,numel(N),[],2);
    bs=bootstrap_macro_statistics(paired,N,'Groups',groups, ...
        'Windows',tailCount,'Bootstrap',bootstrapCount);

    %% Summarize retained time averages only
    for ei=1:numel(estimator)
        for fi=1:numel(names)+legacy
            total=fi>numel(names); name='coefficient_total';
            if ~total,name=names{fi};end
            component=fi; if total,component=3:result.model.p+2;end
            iid=reshape(result.estimates(:,:,component,1,ei),R,numel(N),[]);
            iidV=sum(var(iid,0,1),3);
            testSeconds=serialSeconds;
            if isfield(result,'timing') && ~isempty(result.timing) && ~strcmp(result.timing.testName,name)
                testSeconds(:)=NaN;
            end
            g=(ei-1)*numel(baseGroups)+fi;
            for mi=1:2
                A=reshape(result.estimates(:,:,component,mi,ei),R,numel(N),[]);
                V=sum(var(A,0,1),3);
                for ni=1:numel(N)
                    vrf=iidV(ni)/V(ni); ti=testSeconds(ni,1); tm=testSeconds(ni,mi);
                    low=NaN;high=NaN;valid=0;status="not_max_N";
                    if ni==numel(N) && mi==2
                        low=bs.vrfLow(g);high=bs.vrfHigh(g);valid=bs.vrfValid(g);status=bs.vrfStatus(g);
                    elseif mi==1,status="baseline";end
                    rows(end+1,:)={result.model.dataset,N(ni),estimator{ei},name,result.methods{mi},V(ni),vrf,vrf*ti/tm,mean(result.acceptance(:,ni,mi)),tm,low,high,valid,status}; %#ok<AGROW>
                end

                %% Fit log-log slopes and bootstrap across replications
                for count=min(tailCount,numel(N))
                    w=find(bs.windows==count);slope=bs.slope(g,mi,w);
                    ci=[bs.slopeLow(g,mi,w),bs.slopeHigh(g,mi,w)];
                    fits(end+1,:)={result.model.dataset,estimator{ei},name,result.methods{mi},count,slope,ci(1),ci(2),-1-1/(result.model.p+2),bs.slopeValid(g,mi,w),bs.slopeStatus(g,mi,w)}; %#ok<AGROW>
                end
            end
        end
    end

    %% Write summary, slope and provenance tables
    summary=cell2table(rows,'VariableNames',{'dataset','N','estimator','test','method','variance','vrf','tae','acceptance','seconds','vrf_ci_low','vrf_ci_high','vrf_bootstrap_valid','vrf_bootstrap_status'});
    slopes=cell2table(fits,'VariableNames',{'dataset','estimator','test','method','window','slope','ci_low','ci_high','reference_slope','bootstrap_valid','bootstrap_status'});
    if nargin<4,outputStem=path(1:end-4);end
    stem=outputStem; writetable(summary,[stem '_summary.csv']); writetable(slopes,[stem '_slopes.csv']);
    derived=struct('source',path,'summary',table2struct(summary),'slopes',table2struct(slopes),'bootstrap',bs);
    fid=fopen([stem '_statistics.json'],'w'); fprintf(fid,'%s',jsonencode(derived,'PrettyPrint',true));fclose(fid);
end
