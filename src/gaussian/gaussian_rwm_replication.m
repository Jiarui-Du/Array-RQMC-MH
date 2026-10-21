function out = gaussian_rwm_replication(N,T0,T,mu,Sigma,stepScale,method,seed,options)
%GAUSSIAN_RWM_REPLICATION Multivariate extension of the 1d Gaussian RWM.
% mu is a row vector; Sigma is SPD. Proposal covariance is
% stepScale^2*Sigma. The sole test is the standardized linear projection.
% method: 'hilbert' or 'iid'. Independent scramble every step.
%
% Particle states are N-by-d. Output tests are retained time averages.
% Acceptance is counted over all N*T proposals. Cached Hilbert
% keys are refreshed only for particles accepted at the preceding transition.

    %% Declare particle, target and initialization inputs
    arguments
        N (1,1) double {mustBeInteger,mustBePositive}
        T0 (1,1) double {mustBeInteger,mustBeNonnegative}
        T (1,1) double {mustBeInteger,mustBePositive}
        mu (1,:) double {mustBeFinite,mustBeReal}
        Sigma (:,:) double {mustBeFinite,mustBeReal}
        stepScale (1,1) double {mustBePositive,mustBeFinite}
        method (1,:) char
        seed (1,1) double {mustBeInteger,mustBeNonnegative}
        options.Initialization (1,:) char = 'stationary'
        options.InitialState double = []
        options.HilbertBits double = []
        options.SortCoordinates (1,:) char = 'whitened'
    end

    %% Validate target geometry and construct the random stream
    d = numel(mu);
    bits=resolve_hilbert_bits(d,options.HilbertBits);
    assert(T>T0 && log2(N)==fix(log2(N)), 'Require T>T0 and N=2^m.');
    assert(isequal(size(Sigma),[d,d]) && norm(Sigma-Sigma','fro')<=1e-12*max(1,norm(Sigma,'fro')), 'Sigma must be symmetric d-by-d.');
    L = chol(Sigma,'lower');
    assert(ismember(method,{'hilbert','iid'}),'Unknown method.');
    assert(ismember(options.SortCoordinates,{'whitened','marginal'}),'Unknown sort coordinates.');
    assert(ismember(options.Initialization,{'stationary','zero','fixed'}),'Unknown initialization.');
    stream = RandStream('mt19937ar','Seed',seed);
    previous = RandStream.getGlobalStream;
    cleanup = onCleanup(@() RandStream.setGlobalStream(previous));
    RandStream.setGlobalStream(stream);
    delta = 2^-53;

    %% Initialize particles and log densities
    if strcmp(options.Initialization,'stationary')
        if strcmp(method,'iid')
            Z = randn(stream,N,d);
        else
            U = net(scramble(sobolset(d,'Skip',0,'Leap',0),'MatousekAffineOwen'),N);
            Z = -sqrt(2)*erfcinv(2*min(max(U,delta),1-delta));
        end
        X = mu+Z*L';
    else
        x0 = options.InitialState;
        if strcmp(options.Initialization,'zero')
            assert(isempty(x0),'InitialState must be empty when Initialization is zero');
            x0=zeros(1,d);
        elseif isempty(x0)
            x0=mu;
        end
        assert(isequal(size(x0),[1,d]) && all(isfinite(x0)), 'InitialState must be a finite row.');
        X = repmat(x0,N,1);
        Z = (X-mu)/L';
    end
    logp = -0.5*sum(Z.^2,2);
    base = sobolset(d+2,'Skip',0,'Leap',0,'PointOrder','graycode');
    accumulated = 0; acceptedTotal = 0;
    collisionMax = 0; collisionTotal = 0;
    sortKeys = []; changed = true(N,1);

    %% Sort particles and generate transition inputs
    for s = 1:T
        if strcmp(method,'hilbert')
            if strcmp(options.SortCoordinates,'whitened')
                coordinates = (X-mu)/L';
            else
                coordinates = (X-mu)./sqrt(diag(Sigma))';
            end
            Ustate = 0.5*erfc(-coordinates/sqrt(2));
            if isempty(sortKeys)
                sortKeys = hilbert_keys(Ustate,bits);
            elseif any(changed)
                sortKeys(changed,:) = hilbert_keys(Ustate(changed,:),bits);
            end
            [order,c] = hilbert_order_from_keys(sortKeys,bits);
            sortKeys = sortKeys(order,:);
            X = X(order,:); logp = logp(order);
            collisionMax = max(collisionMax,c); collisionTotal = collisionTotal+c;
        end
        if strcmp(method,'iid')
            innovations = randn(stream,N,d); uaccept = rand(stream,N,1);
        else
            W = net(scramble(base,'MatousekAffineOwen'),N);
            rank = floor(N*W(:,1))+1;
            assert(all(rank>=1 & rank<=N) && all(accumarray(rank,1,[N,1])==1),...
                'Sobol coordinate 1 must be N-stratified.');
            U = zeros(N,d); uaccept = zeros(N,1);
            U(rank,:) = W(:,2:d+1); uaccept(rank) = W(:,d+2);
            innovations = -sqrt(2)*erfcinv(2*min(max(U,delta),1-delta));
        end

        %% Update all coordinates with a joint Metropolis proposal
        Y = X+stepScale*innovations*L';
        ZY = (Y-mu)/L'; logq = -0.5*sum(ZY.^2,2);
        logalpha = min(logq-logp,0);
        accept = isfinite(logalpha) & log(uaccept)<=logalpha;
        X(accept,:) = Y(accept,:); logp(accept) = logq(accept);
        changed = accept;
        acceptedTotal = acceptedTotal+nnz(accept);

        %% Accumulate only the projection on the retained window
        if s>T0
            Z = (X-mu)/L';
            projection = sum(Z,2)/sqrt(d);
            stepMean = mean(projection,1);
            accumulated = accumulated+stepMean;
        end
    end

    %% Return estimates, acceptance and collision counts
    out.average = accumulated/(T-T0);
    out.acceptance = acceptedTotal/(N*T);
    out.collisionMax = collisionMax; out.collisionMean = collisionTotal/T;
    out.testNames = {'projection'};
    out.reference = 0;
    out.hilbertBits = bits;
end
