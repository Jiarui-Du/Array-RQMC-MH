function out = bayesian_rwm_replication(N,T0,T,model,method,seed,bits,recordHistory,diagnostics,testFunction)
% Identical fixed initial state/proposal; only sorting and driving inputs differ.
%
% One replication evolves N parameter rows; averages use transitions T0+1:T.
% The Sobol coordinates encode matching rank, p proposal inputs and acceptance.
% The output contains time-average tests and optional diagnostic history.

    %% Validate options and set diagnostic defaults
    if nargin<10,testFunction='all';end
    [~,testNames]=bayesian_functions(zeros(0,model.p),model,testFunction);
    if nargin<8,recordHistory=false;end
    if nargin<7,bits=[];end
    bits=resolve_hilbert_bits(model.p,bits);
    if nargin<9,diagnostics=struct;end
    if ~isfield(diagnostics,'PrecisionCheck'),diagnostics.PrecisionCheck=false;end
    if ~isfield(diagnostics,'CollisionCheck'),diagnostics.CollisionCheck=false;end
    assert(T>T0 && log2(N)==fix(log2(N)) && ismember(method,{'iid','hilbert'}));

    %% Initialize the random stream, state and cached quantities
    stream=RandStream('mt19937ar','Seed',seed);
    previous=RandStream.getGlobalStream;
    cleanup=onCleanup(@() RandStream.setGlobalStream(previous));
    RandStream.setGlobalStream(stream);
    p=model.p; X=repmat(model.center,N,1);
    marginalSD=sqrt(diag(model.Sigma))';
    lp=bayesian_logposterior(X,model.X,model.y);
    base=sobolset(p+2,'Skip',0,'Leap',0,'PointOrder','graycode');
    accumulated=zeros(1,numel(testNames)); accepted=0; collisionMax=NaN; delta=2^-53;
    if diagnostics.CollisionCheck,collisionMax=0;end
    precisionChanges=zeros(1,2); precisionComparisons=0;
    if recordHistory,history=zeros(T,numel(testNames));end
    sortKeys=[]; values=[]; changed=true(N,1);

    %% Advance particles with IID inputs or Hilbert-matched nets
    for t=1:T
        if strcmp(method,'hilbert')
            if isempty(sortKeys)
                psi=0.5*erfc(-(X-model.center)./marginalSD/sqrt(2));
                sortKeys=hilbert_keys(psi,bits);
            elseif any(changed)
                psi=0.5*erfc(-(X(changed,:)-model.center)./marginalSD/sqrt(2));
                sortKeys(changed,:)=hilbert_keys(psi,bits);
            end
            if diagnostics.CollisionCheck
                [order,c]=hilbert_order_from_keys(sortKeys,bits);
                collisionMax=max(collisionMax,c);
            else
                order=hilbert_order_from_keys(sortKeys,bits);
            end
            if diagnostics.PrecisionCheck
                psi=0.5*erfc(-(X-model.center)./marginalSD/sqrt(2));
                for k=1:2
                    referenceBits=[20 28];
                    alternate=hilbert_sort_order(psi,referenceBits(k));
                    precisionChanges(k)=precisionChanges(k)+nnz(order~=alternate);
                end
                precisionComparisons=precisionComparisons+N;
            end
            X=X(order,:); lp=lp(order);
            sortKeys=sortKeys(order,:);
            if ~isempty(values),values=values(order,:);end
            W=net(scramble(base,'MatousekAffineOwen'),N);
            rank=floor(N*W(:,1))+1;
            assert(all(rank>=1 & rank<=N) && all(accumarray(rank,1,[N,1])==1),...
                'First coordinate must contain exactly one point per N-stratum');
            U=zeros(N,p); ua=zeros(N,1); U(rank,:)=W(:,2:p+1); ua(rank)=W(:,p+2);
            Z=-sqrt(2)*erfcinv(2*min(max(U,delta),1-delta));
        else
            Z=randn(stream,N,p); ua=rand(stream,N,1);
        end

        %% Apply the symmetric-proposal Metropolis update
        Y=X+model.scale*Z*model.L'; lq=bayesian_logposterior(Y,model.X,model.y);
        logalpha=min(0,lq-lp); accept=isfinite(logalpha) & log(ua)<=logalpha;
        X(accept,:)=Y(accept,:); lp(accept)=lq(accept); accepted=accepted+nnz(accept);
        changed=accept;

        %% Accumulate tests over the requested observation window
        if t>T0 || recordHistory
            % Rejected particles retain every test-function value. Initialize
            % at the first observed step, then evaluate only accepted states.
            if N<=2048 || isempty(values)
                values=bayesian_functions(X,model,testFunction);
            elseif any(accept)
                % Preserve the original 2048-column BLAS block shape. A
                % short/single-column product can round differently. Pad
                % the final accepted block, then discard the dummy values.
                moved=X(accept,:); count=size(moved,1);
                moved=[moved;repmat(moved(end,:),mod(-count,2048),1)];
                updated=bayesian_functions(moved,model,testFunction);
                values(accept,:)=updated(1:count,:);
            end
            stepMean=mean(values,1);
            if t>T0,accumulated=accumulated+stepMean;end
            if recordHistory,history(t,:)=stepMean;end
        end
    end

    %% Return estimates and diagnostics
    out=struct('average',accumulated/(T-T0),...
        'acceptance',accepted/(N*T),'maxCollisions',collisionMax,...
        'precisionChangedFraction',precisionChanges/max(1,precisionComparisons),...
        'precisionComparisons',precisionComparisons,'precisionReferenceBits',[20 28], ...
        'hilbertBits',bits,'testNames',{testNames},'stratificationPassed',true);
    if precisionComparisons==0,out.precisionChangedFraction=[NaN NaN];end
    if recordHistory,out.particleMeans=history;end
end
