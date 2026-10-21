function out=gaussian_mwg_replication(N,T0,T,p,stepScale,method,seed,substreamStart,options)
%GAUSSIAN_MWG_REPLICATION Systematic-scan MWG for N_p(0,I_p).
% Array method sorts COMPLETE particles by the active coordinate, then
% matches a fresh scrambled 3D net by its first-coordinate stratum.
% Coordinate r owns substream substreamStart+r-1 throughout the replication.
%
% For the array method, whole particle rows are sorted before each scalar update.
% Projection is formed from coordinate means after complete sweeps.
% ReturnState optionally exposes the final particles and marginal averages.

    %% Declare particle, sweep and randomization inputs
    arguments
        N (1,1) double {mustBeInteger,mustBePositive}
        T0 (1,1) double {mustBeInteger,mustBeNonnegative}
        T (1,1) double {mustBeInteger,mustBePositive}
        p (1,1) double {mustBeInteger,mustBePositive}
        stepScale (1,1) double {mustBePositive,mustBeFinite}
        method (1,:) char {mustBeMember(method,{'iid_mwg','array_mwg'})}
        seed (1,1) double {mustBeInteger,mustBeNonnegative}
        substreamStart (1,1) double {mustBeInteger,mustBePositive}
        options.Initialization (1,:) char {mustBeMember(options.Initialization,{'zero','fixed','stationary'})} = 'zero'
        options.InitialState double = []
        options.ValidateStratification (1,1) logical = false
        options.ReturnState (1,1) logical = false
    end

    %% Validate the grid and allocate coordinate streams
    assert(T>T0 && log2(N)==fix(log2(N)),'Require T>BurnIn and N=2^m');
    assert(seed<2^32 && substreamStart+p-1<=2^51,'Random stream index out of range');
    isArray=strcmp(method,'array_mwg');
    previous=RandStream.getGlobalStream;
    cleanup=onCleanup(@()RandStream.setGlobalStream(previous));
    streams=cell(1,p);
    for r=1:p
        streams{r}=RandStream('mrg32k3a','Seed',seed);
        streams{r}.Substream=substreamStart+r-1;
    end

    %% Initialize the particle array
    X=zeros(N,p); delta=2^-53;
    if strcmp(options.Initialization,'stationary')
        assert(isempty(options.InitialState),'InitialState must be empty for stationary initialization');
        for r=1:p
            RandStream.setGlobalStream(streams{r});
            if isArray
                u=net(scramble(sobolset(1),'MatousekAffineOwen'),N);
                X(:,r)=-sqrt(2)*erfcinv(2*min(max(u,delta),1-delta));
            else
                X(:,r)=randn(streams{r},N,1);
            end
        end
    elseif strcmp(options.Initialization,'fixed')
        x0=options.InitialState;
        if isempty(x0),x0=zeros(1,p);end
        if isscalar(x0),x0=repmat(x0,1,p);end
        assert(isequal(size(x0),[1 p]) && all(isfinite(x0)) && isreal(x0), ...
            'InitialState must be a finite scalar or 1-by-p row');
        X=repmat(x0,N,1);
    else
        assert(isempty(options.InitialState),'InitialState must be empty for zero initialization');
    end
    if isArray,base=sobolset(3,'Skip',0,'Leap',0,'PointOrder','graycode');end
    accumulated=zeros(p,1); accepted=zeros(1,p);

    %% Perform systematic coordinate scans
    for sweep=1:T
        for r=1:p
            stream=streams{r};
            if isArray
                [~,order]=sort(X(:,r),'ascend');
                X=X(order,:); % Preserve pairings: move all coordinates together.
                RandStream.setGlobalStream(stream);
                W=net(scramble(base,'MatousekAffineOwen'),N);
                rank=floor(N*W(:,1))+1;
                if options.ValidateStratification
                    assert(all(rank>=1 & rank<=N) && all(accumarray(rank,1,[N 1])==1), ...
                        'First Sobol coordinate is not exactly N-stratified');
                end
                U=zeros(N,2); U(rank,:)=W(:,2:3);
                innovation=-sqrt(2)*erfcinv(2*min(max(U(:,1),delta),1-delta));
                uaccept=U(:,2);
            else
                innovation=randn(stream,N,1); uaccept=rand(stream,N,1);
            end

            %% Accept or reject the active scalar proposal
            old=X(:,r); proposal=old+stepScale*innovation;
            logalpha=min(-0.5*(proposal.^2-old.^2),0);
            accept=isfinite(logalpha) & log(uaccept)<=logalpha;
            old(accept)=proposal(accept); X(:,r)=old;
            accepted(r)=accepted(r)+nnz(accept);
        end
        if sweep>T0
            % Evaluate after the COMPLETE sweep; averaging is linear.
            accumulated=accumulated+mean(X,1)';
        end
    end

    %% Normalize retained sweep averages and assemble outputs
    marginal=accumulated/(T-T0);
    out.average=sum(marginal)/sqrt(p);
    out.acceptanceByCoordinate=accepted/(N*T);
    out.acceptance=mean(out.acceptanceByCoordinate);
    if options.ReturnState
        out.finalState=X;
        out.marginalAverage=marginal;
    end
end
