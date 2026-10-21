function path=select_result_file(dataset,dimensions,exponents,explicitPath,runRoot)
% Select newest readable complete result satisfying the requested grid.
% Explicit paths are validated without fallback; this function never simulates.
    if nargin<5,root=fileparts(fileparts(fileparts(mfilename('fullpath'))));runRoot=fullfile(root,'results');end
    if strcmp(dataset,'mwg'),dataset='gaussian_mwg';end
    folders=struct('gaussian','gaussian','gaussian_mwg','gaussian_mwg','ripley','bayesian','pima','bayesian');
    patterns=struct('gaussian','gaussian_hilbert_*.mat','gaussian_mwg','gaussian_mwg_v2_*.mat', ...
        'ripley','ripley_formal_v3_*.mat','pima','pima_formal_v3_*.mat');
    assert(isfield(folders,dataset),'Unknown dataset');
    if ~isempty(explicitPath)
        validate(explicitPath,dataset,dimensions,exponents);path=explicitPath;
    else
        entries=dir(fullfile(runRoot,folders.(dataset),patterns.(dataset)));
        entries=entries(~contains(lower(string({entries.name})),'checkpoint'));
        [~,order]=sort([entries.datenum],'descend');path='';
        for i=order
            candidate=fullfile(entries(i).folder,entries(i).name);
            try
                validate(candidate,dataset,dimensions,exponents);path=candidate;break;
            catch err
                fprintf('Skipping %s: %s\n',entries(i).name,err.message);
            end
        end
        assert(~isempty(path),'No complete result satisfies the selected dimensions and N under %s',runRoot);
    end
    [ok,info]=fileattrib(path);assert(ok);path=info.Name;
    fprintf('Reading result: %s\n',path);
end

function validate(path,dataset,dimensions,exponents)
    saved=load(path);
    if ismember(dataset,{'ripley','pima'})
        r=saved.result;assert(strcmp(r.model.dataset,dataset),'Dataset mismatch');
        cases=struct('dimension',r.model.p,'N',r.N,'estimates',r.estimates,'completed',r.completed);
    else
        r=saved.results;cases=r.cases;
        if strcmp(dataset,'gaussian_mwg'),assert(isequal(r.methods,{'iid_mwg','array_mwg'}));
        else,assert(isequal(r.methods,{'iid','hilbert'}));end
    end
    if isempty(dimensions),dimensions=[cases.dimension];end
    assert(~isempty(dimensions) && all(ismember(dimensions,[cases.dimension])),'Requested dimensions missing');
    for d=dimensions
        c=cases([cases.dimension]==d);assert(isscalar(c),'Duplicate dimension');
        assert(all(ismember(2.^exponents,c.N)),'Requested N missing');
        % Slopes use the source grid, so all its time-average samples must be complete.
        A=c.estimates(:,:,:,:,1);
        assert(size(A,1)>=2 && size(A,2)==numel(c.N) && size(A,4)==2 && all(isfinite(A),'all'),'Incomplete estimates');
        if isfield(c,'completed'),assert(all(c.completed,'all'),'Incomplete replications');end
    end
    if strcmp(dataset,'gaussian_mwg')
        summarize_gaussian_mwg(r,'Dimensions',dimensions,'NExponents',exponents,'Bootstrap',0);
    else
        read_convergence_statistics(path,'Dimensions',dimensions,'NExponents',exponents,'Bootstrap',0);
    end
end
