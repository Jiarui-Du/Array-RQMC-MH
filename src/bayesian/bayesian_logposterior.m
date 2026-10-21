function [lp,g,H] = bayesian_logposterior(beta,X,y)
% Flat-prior logistic posterior; rows of beta are parameter vectors.
%
% Returns one log posterior per parameter row, up to an additive constant.
% Optional gradient g and Hessian H are defined only for a single row of beta.

    %% Evaluate large particle arrays in bounded blocks
    if size(beta,1)>2048
        assert(nargout==1,'Derivatives are defined for a single parameter vector');
        lp=zeros(size(beta,1),1);
        for first=1:2048:size(beta,1)
            block=first:min(first+2047,size(beta,1));
            lp(block)=bayesian_logposterior(beta(block,:),X,y);
        end
        return
    end

    %% Compute the log likelihood without exponential overflow
    eta = X*beta';
    signed = (2*y-1).*eta;
    lp = -sum(max(-signed,0)+log1p(exp(-abs(signed))),1)';

    %% Compute derivatives for mode finding
    if nargout>1
        assert(size(beta,1)==1);
        prob = zeros(size(eta)); mask=eta>=0;
        prob(mask)=1./(1+exp(-eta(mask)));
        z=exp(eta(~mask)); prob(~mask)=z./(1+z);
        g=X'*(y-prob);
        H=-(X'*(X.*(prob.*(1-prob))));
    end
end
