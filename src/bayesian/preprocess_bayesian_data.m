function [X,y,details] = preprocess_bayesian_data(dataset)
% Explicit reconstruction from unstandardized MASS CSV source files.
%
% Returns the reconstructed author design X, binary y and preprocessing metadata.
% This is a consistency check against the author matrices, not a data-file rewrite.

    %% Load the dataset-specific source columns
    root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
    folder=fullfile(root,'data',dataset,'source');
    if strcmp(dataset,'ripley')
        tab=readtable(fullfile(folder,'synth.tr.csv'),'TextType','string');
        covariates={'xs','ys'}; raw=tab{:,covariates}; y=double(tab.yc);
    else
        a=readtable(fullfile(folder,'Pima.tr.csv'),'TextType','string');
        b=readtable(fullfile(folder,'Pima.te.csv'),'TextType','string'); tab=[a;b];
        covariates={'npreg','glu','bp','skin','bmi','ped','age'};
        raw=tab{:,covariates}; y=double(tab.type=="Yes");
    end

    %% Standardize, reorder and reproduce published rounding
    assert(all(isfinite(raw),'all'),'Selected MASS data must already be complete');
    center=mean(raw,1); scale=std(raw,0,1); Z=(raw-center)./scale;
    match=jsondecode(fileread(fullfile(folder,'source_match.json')));
    order=match.source_row_for_author_row; Z=Z(order,:); y=y(order);
    if strcmp(dataset,'pima'), Z=round(Z,5,'significant'); end
    X=[ones(size(Z,1),1),Z];

    %% Verify reconstruction against the author matrices
    model=load_bayesian_model(dataset);
    assert(max(abs(X-model.X),[],'all')<1e-12 && isequal(y,model.y),...
        'Reconstruction differs from author matrix; halt experiments');
    details=struct('covariates',{covariates},'sourceMean',center,'sourceSD',scale,...
        'rowOrder',order,'roundSignificantDigits',5*strcmp(dataset,'pima'),...
        'interceptPosition',1,'deletedRowsAtThisStage',[]);
end
