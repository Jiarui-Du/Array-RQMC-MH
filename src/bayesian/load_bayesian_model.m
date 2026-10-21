function model = load_bayesian_model(dataset)
% Use exact author text matrices, preserving published file rounding.
%
% Returns the design, response, MLE center and inverse observed-information
% covariance. L is its lower Cholesky factor; the proposal uses scale*Z*L'.

    %% Load author data and place the intercept first
    root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
    folder=fullfile(root,'data',dataset,'author');
    raw=readmatrix(fullfile(folder,'myX.txt'),'FileType','text');
    y=readmatrix(fullfile(folder,'Y.txt'),'FileType','text'); y=y(:);
    p=size(raw,2); X=raw(:,[p,1:p-1]); % author's initial_inds
    assert(all(raw(:,end)==1) && all(ismember(y,[0 1])) && all(isfinite(X),'all'));

    %% Find the mode by Newton steps with backtracking
    beta=zeros(1,p);
    for it=1:100
        [lp,g,H]=bayesian_logposterior(beta,X,y);
        if norm(g,Inf)<1e-9, break; end
        direction=(-H)\g; step=1;
        while bayesian_logposterior(beta+step*direction',X,y)<lp+1e-4*step*(g'*direction)
            step=step/2; assert(step>2^-40,'Newton line search failed');
        end
        beta=beta+step*direction';
    end
    [lp,g,H]=bayesian_logposterior(beta,X,y);
    assert(norm(g,Inf)<1e-7,'MLE did not converge');

    %% Build proposal geometry and numerical diagnostics
    Sigma=(-H)\eye(p); Sigma=(Sigma+Sigma')/2; L=chol(Sigma,'lower');
    model=struct('dataset',dataset,'X',X,'y',y,'p',p,'center',beta,'Sigma',Sigma,'L',L,...
        'scale',2.38/sqrt(p),'prior','improper_flat','logMode',lp,'gradientNorm',norm(g),...
        'eigenvalues',eig(Sigma),'condition',cond(Sigma),...
        'cholError',norm(L*L'-Sigma,'fro')/norm(Sigma,'fro'),...
        'sourceDirectory',folder,'root',root);
end
