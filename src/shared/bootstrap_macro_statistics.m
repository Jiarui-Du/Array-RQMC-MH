function out=bootstrap_macro_statistics(A,N,options)
%BOOTSTRAP_MACRO_STATISTICS Percentile intervals from whole macro estimates.
% A axes: replication, N, trajectory output, method (IID then array).
% Outputs on one trajectory share indices; N and methods are independent
% unless their explicit pairing flags are true. Groups sum component variances.
% No particles, time points, or timing observations are resampled here.
    arguments
        A double
        N double
        options.Bootstrap (1,1) double {mustBeInteger,mustBeNonnegative} = 1000
        options.Seed (1,1) double {mustBeInteger,mustBeNonnegative} = 19002026
        options.Windows double = 5
        options.Groups cell = {}
        options.PairedN (1,1) logical = false
        options.PairedMethods (1,1) logical = false
    end
    %% Validate dimensions and fix the analysis windows before resampling
    N=N(:)';[R,K,C,M]=size(A);
    assert(R>=2 && K==numel(N) && M==2 && all(isfinite(A),'all'));
    assert(all(N>0) && all(diff(N)>0));
    counts=unique(min(options.Windows,K));
    assert(all(counts>=1 & counts==fix(counts)));
    groups=options.Groups;if isempty(groups),groups=num2cell(1:C);end
    for g=1:numel(groups)
        assert(~isempty(groups{g}) && all(ismember(groups{g},1:C)));
    end
    B=options.Bootstrap;G=numel(groups);W=numel(counts);
    out=struct('protocol','macro-percentile-v1','count',B,'seed',options.Seed, ...
        'confidence',0.95,'minimumValidFraction',0.9,'pairedN',options.PairedN, ...
        'pairedMethods',options.PairedMethods,'pairedOutputs',true, ...
        'N',N,'windows',counts,'groups',{groups},'maxN',N(end));
    V=group_variance(A,groups,K,G);
    out.variance=V;
    out.slope=slopes(V,N,counts,G,W);
    out.vrf=ratio(reshape(V(end,:,1),G,1),reshape(V(end,:,2),G,1));
    out.slopeSamples=nan(B,G,2,W);out.vrfSamples=nan(B,G);
    stream=RandStream('mt19937ar','Seed',mod(options.Seed,2^32));

    %% Draw independent blocks, preserving every within-trajectory output
    for b=1:B
        sampled=zeros(size(A));indices=zeros(R,K,2);
        for m=1:2
            for k=1:K
                if options.PairedMethods && m==2
                    index=indices(:,k,1);
                elseif options.PairedN && k>1
                    index=indices(:,1,m);
                else
                    index=randi(stream,R,R,1);
                end
                indices(:,k,m)=index;sampled(:,k,:,m)=A(index,k,:,m);
            end
        end
        bv=group_variance(sampled,groups,K,G);
        out.slopeSamples(b,:,:,:)=reshape(slopes(bv,N,counts,G,W),1,G,2,W);
        out.vrfSamples(b,:)=ratio(bv(end,:,1),bv(end,:,2));
    end

    %% Apply one degeneracy rule to all slope and VRF intervals
    [out.slopeLow,out.slopeHigh,out.slopeValid,out.slopeStatus]= ...
        intervals(out.slopeSamples,out.slope,B);
    [out.vrfLow,out.vrfHigh,out.vrfValid,out.vrfStatus]= ...
        intervals(out.vrfSamples,out.vrf,B);
end

function V=group_variance(A,groups,K,G)
% Trace groups retain component pairing without replacing trace by var(sum).
    V=zeros(K,G,2);
    for g=1:G
        V(:,g,:)=reshape(sum(var(A(:,:,groups{g},:),0,1),3),K,1,2);
    end
end

function s=slopes(V,N,counts,G,W)
% Fixed windows: never silently remove a zero-variance N from a fit.
    s=nan(G,2,W);
    for w=1:W
        tail=numel(N)-counts(w)+1:numel(N);
        if numel(tail)<2,continue;end
        for m=1:2
            for g=1:G
                v=V(tail,g,m);
                if all(isfinite(v) & v>0)
                    f=polyfit(log2(N(tail)),log2(v(:)'),1);s(g,m,w)=f(1);
                end
            end
        end
    end
end

function r=ratio(a,b)
% Both variances must be finite and positive for a nondegenerate VRF.
    r=nan(size(a));ok=isfinite(a)&isfinite(b)&a>0&b>0;r(ok)=a(ok)./b(ok);
end

function [lo,hi,valid,status]=intervals(samples,point,B)
% Invalid draws are counted, not repaired with eps or redrawn conditionally.
    lo=nan(size(point));hi=lo;valid=zeros(size(point));
    status=repmat("disabled",size(point));
    if B==0,return;end
    samples=reshape(samples,B,[]);
    for j=1:numel(point)
        v=samples(:,j);v=v(isfinite(v));valid(j)=numel(v);
        if ~isfinite(point(j))
            status(j)="degenerate_original";
        elseif numel(v)<max(2,ceil(0.9*B))
            status(j)="insufficient_valid";
        else
            q=quantile(v,[.025 .975]);lo(j)=q(1);hi(j)=q(2);status(j)="ok";
        end
    end
end
