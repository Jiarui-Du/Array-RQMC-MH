function [values,names] = bayesian_functions(beta,model,selection)
%BAYESIAN_FUNCTIONS Evaluate selected tests; unselected tests are not computed.
% Projection uses the MLE and diagonal of inverse observed information.
% Paper drivers select projection; other branches support historical results and tests.
% In those legacy branches, beta_1 is the intercept.
    %% Resolve the requested output columns
    if nargin<3,selection='all';end
    available=[{'predictive','projection'},arrayfun(@(j)sprintf('beta_%d',j),1:model.p,'UniformOutput',false)];
    if strcmp(selection,'all'),names=available;else
        assert(ismember(selection,available),'Unknown Bayesian test function');names={char(selection)};
    end
    %% Compute only the selected tests
    values=zeros(size(beta,1),numel(names));
    for j=1:numel(names)
        name=names{j};
        if strcmp(name,'projection')
            values(:,j)=sum((beta-model.center)./sqrt(diag(model.Sigma))',2)/sqrt(model.p);
        elseif startsWith(name,'beta_')
            index=sscanf(name,'beta_%d');values(:,j)=beta(:,index);
        else
            % Keep the established BLAS block width for predictive regression.
            for first=1:2048:size(beta,1)
                block=first:min(first+2047,size(beta,1));
                eta=model.X*beta(block,:)';prob=zeros(size(eta));mask=eta>=0;
                prob(mask)=1./(1+exp(-eta(mask)));e=exp(eta(~mask));prob(~mask)=e./(1+e);
                values(block,j)=mean(prob,1)';
            end
        end
    end
end
