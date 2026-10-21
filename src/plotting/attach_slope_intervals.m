function V=attach_slope_intervals(V,source,tailCount,B)
% Attach both methods' slopes using the source grid, not displayed N rows.
% Thus a max-N VRF view still fits the largest TailCount available N values.
    saved=load(source);
    if isfield(saved,'results')
        cases=saved.results.cases;
    else
        r=saved.result;
        if isfield(r,'testNames'),names=r.testNames;else
            names=[{'smooth','indicator'},arrayfun(@(j)sprintf('beta_%d',j),1:r.model.p,'UniformOutput',false),{'predictive'}];
        end
        cases=struct('dimension',r.model.p,'N',r.N,'testNames',{names},'estimates',r.estimates);
    end
    %% Allocate slope values, confidence limits and diagnostics
    for method={'iid','array'}
        for field={'slope','slope_ci_low','slope_ci_high'}
            V.([field{1} '_' method{1}])=nan(height(V),1);
        end
        V.(['slope_valid_' method{1}])=zeros(height(V),1);
        V.(['slope_status_' method{1}])=repmat("unavailable",height(V),1);
    end
    V.slope_N_min=nan(height(V),1);V.slope_N_max=V.slope_N_min;
    %% Resample all paired test outputs once per selected case
    for c=cases(:)'
        if ~any(V.d==c.dimension),continue;end
        A=c.estimates(:,:,:,:,1);
        groups=num2cell(1:numel(c.testNames));
        if any(strcmp(V.test,'coefficient_total'))
            c.testNames{end+1}='coefficient_total';groups{end+1}=3:c.dimension+2;
        end
        bs=bootstrap_macro_statistics(A,c.N,'Groups',groups,'Windows',tailCount,'Bootstrap',B);
        for fi=1:numel(c.testNames)
            rows=V.d==c.dimension & strcmp(V.test,c.testNames{fi});
            for mi=1:2
                methods={'iid','array'};m=methods{mi};
                V.(['slope_' m])(rows)=bs.slope(fi,mi);
                V.(['slope_ci_low_' m])(rows)=bs.slopeLow(fi,mi);
                V.(['slope_ci_high_' m])(rows)=bs.slopeHigh(fi,mi);
                V.(['slope_valid_' m])(rows)=bs.slopeValid(fi,mi);
                V.(['slope_status_' m])(rows)=bs.slopeStatus(fi,mi);
            end
            V.slope_N_min(rows)=c.N(max(1,numel(c.N)-tailCount+1));
            V.slope_N_max(rows)=c.N(end);
        end
    end
end
