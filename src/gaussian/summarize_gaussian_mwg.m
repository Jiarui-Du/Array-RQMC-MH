function [summary,slopes]=summarize_gaussian_mwg(results,options)
% Summarize selected complete time averages. No plots or simulations.
% VRF and TAE compare Array-MWG with IID-MWG, never with joint RWM.
%
% Selects complete saved cases and returns summary/slopes tables without file writes.
% Variances are across macro replications; seconds are independent serial medians.
% Legacy files remain readable, with missing independent times reported as NaN.
% Bootstrap=0 skips confidence intervals; empirical slopes do not establish a rate.

    %% Declare case selection and analysis options
    arguments
        results (1,1) struct
        options.Dimensions (1,:) double = []
        options.NExponents (1,:) double = []
        options.TailCount (1,1) double {mustBeInteger,mustBePositive} = 5
        options.Bootstrap (1,1) double {mustBeInteger,mustBeNonnegative} = 1000
    end

    %% Validate the result schema and requested grid
    assert(ismember(results.schema,{'gaussian-mwg-v1','gaussian-mwg-v2'}),'Expected Gaussian MWG results');
    assert(isequal(results.methods,{'iid_mwg','array_mwg'}));
    dims=options.Dimensions;if isempty(dims),dims=[results.cases.dimension];end
    assert(~isempty(dims) && all(ismember(dims,[results.cases.dimension])));
    assert(numel(unique(dims))==numel(dims),'Dimensions must not repeat');
    exponents=options.NExponents;if isempty(exponents),exponents=results.settings.NExponents;end
    assert(~isempty(exponents) && all(diff(exponents)>0) && all(exponents==fix(exponents)));
    assert(options.TailCount>=2);N=2.^exponents;rows={};fits={};

    %% Summarize selected complete replications
    for p=dims
        c=results.cases([results.cases.dimension]==p);[present,idx]=ismember(N,c.N);
        assert(all(present),'Requested particle counts missing for p=%d',p);
        for ni=idx
            for mi=1:2
                assert(all(c.completed(:,ni,mi)), ...
                    'p=%d N=%d %s incomplete: %d/%d replications saved', ...
                    p,c.N(ni),results.methods{mi},nnz(c.completed(:,ni,mi)),size(c.completed,1));
            end
        end
        R=size(c.estimates,1);assert(R>=3);
        seconds=standalone_timing_seconds(c,N,results.methods);
        bs=bootstrap_macro_statistics(c.estimates(:,idx,:,:),N, ...
            'Bootstrap',options.Bootstrap,'Windows',options.TailCount);
        bsFull=bootstrap_macro_statistics(c.estimates,c.N,'Bootstrap',options.Bootstrap);
        for fi=1:numel(c.testNames)
            A=reshape(c.estimates(:,idx,fi,:),R,numel(N),2);
            assert(all(isfinite(A),'all'),'Missing or nonfinite selected estimates');
            V=reshape(var(A,0,1),numel(N),2);
            means=reshape(mean(A,1),numel(N),2);
            for mi=1:2
                v=V(:,mi);tail=max(1,numel(N)-options.TailCount+1):numel(N);
                slope=bs.slope(fi,mi);ci=[bs.slopeLow(fi,mi),bs.slopeHigh(fi,mi)];
                fits(end+1,:)={p,c.testNames{fi},results.methods{mi},slope,ci(1),ci(2), ...
                    -4/3,N(tail(1)),N(tail(end)),bs.slopeValid(fi,mi),bs.slopeStatus(fi,mi)}; %#ok<AGROW>
                for ni=1:numel(N)
                    vrf=V(ni,1)/v(ni);tae=vrf*seconds(ni,1)/seconds(ni,mi);
                    low=NaN;high=NaN;valid=0;status="not_max_N";
                    if N(ni)==max(c.N) && mi==2
                        low=bsFull.vrfLow(fi);high=bsFull.vrfHigh(fi);valid=bsFull.vrfValid(fi);status=bsFull.vrfStatus(fi);
                    elseif mi==1,status="baseline";end
                    rows(end+1,:)={p,N(ni),c.testNames{fi},results.methods{mi}, ...
                        means(ni,mi),v(ni),vrf,seconds(ni,mi),tae,low,high,valid,status}; %#ok<AGROW>
                end
            end
        end
    end

    %% Return variance, efficiency and slope tables
    summary=cell2table(rows,'VariableNames', ...
        {'p','N','test','method','estimate_mean','variance','VRF','seconds','TAE','vrf_ci_low','vrf_ci_high','vrf_bootstrap_valid','vrf_bootstrap_status'});
    slopes=cell2table(fits,'VariableNames', ...
        {'p','test','method','variance_slope','bootstrap_low','bootstrap_high', ...
        'reference_slope','fit_N_min','fit_N_max','bootstrap_valid','bootstrap_status'});
end
