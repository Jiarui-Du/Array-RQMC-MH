function [results,matFile]=load_array_rqmc_horizon_result(matFile)
% Select only complete aggregate results, never individual files/checkpoints.
    root=fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
    explicit=~isempty(matFile);
    if explicit
        files={matFile};
    else
        listing=dir(fullfile(root,'results','one_dim_horizon','**','*.mat'));
        [~,order]=sort([listing.datenum],'descend');listing=listing(order);
        files=arrayfun(@(x)fullfile(x.folder,x.name),listing,'UniformOutput',false);
    end
    for k=1:numel(files)
        if contains(files{k},'_checkpoint'),continue;end
        try
            data=load(files{k},'results');r=data.results;
            assert(all(isfield(r,{'variance_array','scaled_variance_array','sigma','target_id','test_id'})));
            assert(all(isfinite(r.variance_array),'all'));
            if ~isfield(r,'individual_results')
                assert(isfield(r,'reproducibility') && isfield(r,'T_vec') && ...
                    isequal(size(r.variance_array),[numel(r.T_vec),numel(r.N_vec)]) && ...
                    all(isfinite(r.scaled_variance_array),'all'), 'Incomplete aggregate result.');
            else
              for j=1:numel(r.individual_results)
                s=r.individual_results{j};
                assert(s.completed_idx==numel(r.T_vec) && all(isfinite(s.est_sobol_all),'all') && all(isfinite(s.est_mcmc_all),'all'));
              end
            end
            if ~explicit
                assert(strcmp(r.target_id,'gaussian_mixture3') && strcmp(r.test_id,'cdf'));
            end
            results=r;matFile=files{k};fprintf('Reading result: %s\n',matFile);return
        catch exception
            if explicit,rethrow(exception);end
        end
    end
    error('array_rqmc:NoHorizonResult','No complete matching horizon result found.');
end
