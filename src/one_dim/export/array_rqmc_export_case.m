function [row, problem] = array_rqmc_export_case( ...
        file_index, target_id, test_id, opt)
%ARRAY_RQMC_EXPORT_CASE Summarize one target/test result file.
%
% Returns a fixed-schema table row plus a problem message. Missing or invalid
% results are represented explicitly. TAE prefers recomputation from median
% elapsed times and falls back to a positive stored value when needed.


    %% Initialize the row and select a saved candidate
    row = template_row();
    row.TargetID = target_id;
    row.Target = array_rqmc_export_label('target_plain', target_id);
    row.TestFunctionID = test_id;
    row.TestFunction = array_rqmc_export_label('test_plain', test_id);
    row.TargetN = double(opt.TargetN);
    row.TaESource = 'unavailable';
    problem = '';

    [selected, selection_status] = array_rqmc_select_result_file( ...
        file_index, target_id, test_id, opt.RequiredR, opt.TargetN);
    if isempty(selected)
        row.Status = selection_status;
        problem = selection_status;
        return;
    end
    row.ResultFile = selected.Path{1};
    row.Replications = selected.Replications(1);
    try
        loaded = load(row.ResultFile, 'results');
        if ~isfield(loaded, 'results')
            error('Missing results struct.');
        end
        results = loaded.results;
        [N, variance_iid, variance_array] = ...
            array_rqmc_extract_variance_data(results, false);
    catch ME
        row.Status = ['invalid result schema: ', ME.message];
        problem = row.Status;
        return;
    end


    %% Fit the variance slope and retain CDF metadata
    valid_N = N(isfinite(N) & N > 0);
    if ~isempty(valid_N)
        row.NMaxCompleted = max(valid_N);
    end
    fit = array_rqmc_fit_loglog_slope(N, variance_array, opt.TailCount);
    row.ArraySlope = fit.slope;
    row.SlopeSE = fit.standard_error;
    row.TailPointCount = fit.n_points;
    if ~isempty(fit.x_used)
        row.NMinSlope = fit.x_used(1);
        row.NMaxSlope = fit.x_used(end);
    end
    % Recompute from raw macro estimates, including legacy saved results.
    if isfield(results,'est_mcmc') && isfield(results,'est_sobol')
        bs=bootstrap_macro_statistics(cat(4,results.est_mcmc,results.est_sobol),N, ...
            'Bootstrap',opt.Bootstrap,'Windows',opt.TailCount);
        row.ArraySlope=bs.slope(1,2);row.SlopeCILow=bs.slopeLow(1,2);row.SlopeCIHigh=bs.slopeHigh(1,2);
        row.SlopeBootstrapValid=bs.slopeValid(1,2);row.SlopeBootstrapStatus=char(bs.slopeStatus(1,2));
        row.BootstrapCount=bs.count;row.BootstrapSeed=bs.seed;
        if opt.TargetN==max(N)
            row.VRFCILow=bs.vrfLow(1);row.VRFCIHigh=bs.vrfHigh(1);
            row.VRFBootstrapValid=bs.vrfValid(1);row.VRFBootstrapStatus=char(bs.vrfStatus(1));
        else,row.VRFBootstrapStatus='not_max_N';end
    end
    row = read_cdf_metadata(row, results);


    %% Compute VRF and timing efficiency at the requested N
    position = find(abs(N - double(opt.TargetN)) <= ...
        8 * eps(double(opt.TargetN)), 1, 'last');
    issues = cell(0, 1);
    if ~isempty(position) && positive(variance_iid(position)) && ...
            positive(variance_array(position))
        row.VRF = variance_iid(position) / variance_array(position);
        timing = array_rqmc_extract_timing_at_N( ...
            results, double(opt.TargetN), position);
        row.StoredTaE = positive_or_nan(timing.stored_tae);
        row.MedianTimeIIDSeconds = positive_or_nan(timing.median_iid);
        row.MedianTimeArraySeconds = positive_or_nan(timing.median_array);
        row.TimingReplications = timing.R;
        row.TimingScope = timing.scope;
        if positive(row.MedianTimeIIDSeconds) && ...
                positive(row.MedianTimeArraySeconds)
            row.TimeRatioIIDToArray = row.MedianTimeIIDSeconds / ...
                row.MedianTimeArraySeconds;
            row.RecomputedTaE = row.VRF * row.TimeRatioIIDToArray;
            row.TaE = row.RecomputedTaE;
            row.TaESource = timing.source;
        elseif positive(row.StoredTaE)
            row.TaE = row.StoredTaE;
            row.TaESource = [timing.source, ' (stored TaE)'];
        else
            issues{end+1} = 'standalone timing unavailable'; %#ok<AGROW>
        end
        if positive(row.StoredTaE) && positive(row.RecomputedTaE)
            row.TaERelativeDifference = abs( ...
                row.StoredTaE - row.RecomputedTaE) / ...
                max(abs(row.RecomputedTaE), realmin);
        end
    else
        issues{end+1} = sprintf('N=%g unavailable', opt.TargetN); %#ok<AGROW>
    end

    %% Report unavailable quantities without inventing values
    if ~isfinite(row.ArraySlope)
        issues{end+1} = 'slope unavailable'; %#ok<AGROW>
    end
    if isempty(issues)
        row.Status = 'complete';
    else
        row.Status = strjoin(issues, '; ');
        problem = row.Status;
    end
end

%% Local helper: template_row
% Define all export columns with explicit missing-value defaults.
function row = template_row()
    row = struct('TargetID','','Target','','TestFunctionID','', ...
        'TestFunction','','TargetN',NaN,'VRF',NaN,'TaE',NaN, ...
        'BootstrapCount',NaN,'BootstrapSeed',NaN, ...
        'SlopeCILow',NaN,'SlopeCIHigh',NaN,'SlopeBootstrapValid',0,'SlopeBootstrapStatus','raw_unavailable', ...
        'VRFCILow',NaN,'VRFCIHigh',NaN,'VRFBootstrapValid',0,'VRFBootstrapStatus','raw_unavailable', ...
        'ArraySlope',NaN,'SlopeSE',NaN,'TailPointCount',0, ...
        'NMinSlope',NaN,'NMaxSlope',NaN,'NMaxCompleted',NaN, ...
        'Replications',NaN,'TimingReplications',NaN,'StoredTaE',NaN, ...
        'RecomputedTaE',NaN,'MedianTimeIIDSeconds',NaN, ...
        'MedianTimeArraySeconds',NaN,'TimeRatioIIDToArray',NaN, ...
        'TaERelativeDifference',NaN,'TaESource','', ...
        'TimingScope','','CDFProbability',NaN,'CDFThreshold',NaN, ...
        'Status','','ResultFile','');
end

%% Local helper: read_cdf_metadata
% Copy saved CDF probability and threshold metadata when available.
function row = read_cdf_metadata(row, results)
    if isfield(results, 'settings') && isstruct(results.settings)
        settings = results.settings;
        if isfield(settings, 'cdf_probability') && ...
                isnumeric(settings.cdf_probability) && ...
                isscalar(settings.cdf_probability)
            row.CDFProbability = double(settings.cdf_probability);
        end
        if isfield(settings, 'cdf_threshold') && ...
                isnumeric(settings.cdf_threshold) && ...
                isscalar(settings.cdf_threshold)
            row.CDFThreshold = double(settings.cdf_threshold);
        end
    end
end

%% Local export-value checks
% Accept only positive finite values for efficiency calculations and export.
function tf = positive(value), tf = isfinite(value) && value > 0; end
%% Local helper: positive_or_nan
% Convert nonpositive or nonfinite export values to missing values.
function value = positive_or_nan(value)
    value = double(value);
    if ~positive(value), value = NaN; end
end
