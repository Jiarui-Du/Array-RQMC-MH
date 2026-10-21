function V=view_variance_reduction(dataset,options)
% Display VRF and independent serial median elapsed seconds for matching kernels.
% Saved tic/toc timings are elapsed time, not process CPU time.
% Legacy coefficient_total uses the covariance trace; paper entries select projection.
%
% Returns one row per selected case/test/N, joining each array estimate with
% its matching IID estimate. CsvFile optionally writes this table; Display only
% controls console output. Missing timing propagates to TAE as NaN.

    %% Declare dataset, test and export options
    arguments
        dataset (1,:) char = 'ripley'
        options.TailCount (1,1) double {mustBeInteger,mustBePositive} = 5
        options.Bootstrap (1,1) double {mustBeInteger,mustBeNonnegative} = 1000
        options.Dimension (1,:) double = 8
        options.NExponents (1,:) double = 10:18
        options.Tests = 'all'
        options.MatFile (1,:) char = ''
        options.Display (1,1) logical = true
        options.CsvFile (1,:) char = ''
    end

    %% Load the selected schema and separate matching methods
    dataset=lower(dataset);
    if strcmp(dataset,'mwg'),dataset='gaussian_mwg';end
    dims=options.Dimension;if ismember(dataset,{'ripley','pima'}),dims=[];end
    source=select_result_file(dataset,dims,options.NExponents,options.MatFile);
    if strcmp(dataset,'gaussian_mwg')
        saved=load(source,'results');
        S=summarize_gaussian_mwg(saved.results,'Dimensions',options.Dimension, ...
            'NExponents',options.NExponents,'Bootstrap',options.Bootstrap);
        S.Properties.VariableNames{strcmp(S.Properties.VariableNames,'p')}='d';
        S.Properties.VariableNames{strcmp(S.Properties.VariableNames,'VRF')}='vrf';
        S.dataset=repmat({dataset},height(S),1);
        S.R=repmat(saved.results.settings.Replications,height(S),1);
        S.estimator=repmat({'time_average'},height(S),1);
        H=S(strcmp(S.method,'array_mwg'),:);I=S(strcmp(S.method,'iid_mwg'),:);
        arrayLabel='Array-MWG';
    else
        dims=[];if strcmp(dataset,'gaussian'),dims=options.Dimension;end
        out=read_convergence_statistics(source,'NExponents',options.NExponents,'Dimensions',dims,'Bootstrap',options.Bootstrap);
        S=out.summary;assert(all(strcmp(S.dataset,dataset)),'MAT dataset does not match selection');
        H=S(strcmp(S.method,'hilbert'),:);I=S(strcmp(S.method,'iid'),:);
        arrayLabel='Hilbert';
    end

    %% Join methods by case identity and compute time-adjusted efficiency
    keys={'dataset','d','N','R','estimator','test'};
    H=H(:,[keys,{'variance','vrf','seconds','vrf_ci_low','vrf_ci_high','vrf_bootstrap_valid'}]);
    H.Properties.VariableNames(numel(keys)+(1:3))={'variance_hilbert','VRF','seconds_hilbert'};
    I=I(:,[keys,{'variance','seconds'}]);I.Properties.VariableNames(end-1:end)={'variance_iid','seconds_iid'};
    V=innerjoin(I,H,'Keys',keys);
    if strcmp(dataset,'gaussian_mwg')
        V.Properties.VariableNames{strcmp(V.Properties.VariableNames,'variance_hilbert')}='variance_mwg';
        V.Properties.VariableNames{strcmp(V.Properties.VariableNames,'seconds_hilbert')}='seconds_mwg';
        timeColumn='seconds_mwg';
    else
        timeColumn='seconds_hilbert';
    end
    V.TAE=V.VRF.*V.seconds_iid./V.(timeColumn);

    %% Filter tests and emit the requested outputs
    wanted=string(options.Tests);
    if ~isequal(wanted,"all")
        assert(all(ismember(wanted,string(V.test))),'Unknown test function');
        V=V(ismember(string(V.test),wanted),:);
    end
    V=attach_slope_intervals(V,source,options.TailCount,options.Bootstrap);
    V(:,contains(lower(string(V.Properties.VariableNames)),'status'))=[];
    V=sortrows(V,{'d','test','N'});
    if options.Display
        fprintf('%s: VRF = Var(matching IID)/Var(%s).\n',dataset,arrayLabel);
        fprintf('seconds: independent serial single-test median elapsed time; NaN means unavailable.\n');
        fprintf('TAE = VRF * seconds_iid / %s; VRF or TAE > 1 favors %s.\n',timeColumn,arrayLabel);
        if any(isnan(V.TAE)),fprintf('Legacy concurrent timings and unmeasured N are not used for TAE.\n');end
        disp(V(:,{'dataset','d','N','test','VRF','vrf_ci_low','vrf_ci_high','seconds_iid',timeColumn,'TAE'}));
        fprintf('Slopes use the largest TailCount N values in the source file, independently of displayed N.\n');
        disp(unique(V(:,{'dataset','d','test','slope_N_min','slope_N_max', ...
            'slope_iid','slope_ci_low_iid','slope_ci_high_iid',...
            'slope_array','slope_ci_low_array','slope_ci_high_array'}),'rows','stable'));
    end
    if ~isempty(options.CsvFile)
        folder=fileparts(options.CsvFile);
        if ~isempty(folder) && ~isfolder(folder),mkdir(folder);end
        writetable(V,options.CsvFile);
    end
end
