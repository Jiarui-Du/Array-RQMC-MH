function results = summarize_array_rqmc_statistics( ...
        N_vec, m_vec, data, settings, reference_values, slope_tail_count)
%SUMMARIZE_ARRAY_RQMC_STATISTICS Assemble raw and summary statistics.
%
% data contains R-by-N-by-test estimates and R-by-N acceptance arrays.
% Variances use sample normalization across independent macro replications.
% Bias and MSE are populated only when an explicit finite reference is supplied.


    %% Normalize grid and reference metadata
    N_vec = double(N_vec(:));
    m_vec = double(m_vec(:));
    reference_values = reshape(double(reference_values), 1, []);
    N_count = numel(N_vec);
    test_count = size(data.est_mcmc, 3);
    if numel(reference_values) ~= test_count
        error('array_rqmc:ReferenceValues', ...
            'ReferenceValues must have one entry per test function.');
    end


    %% Allocate raw and derived result fields
    results = struct;
    results.N_vec = N_vec;
    results.m_vec = m_vec;
    results.completed_idx = N_count;
    results.test_ids = settings.test_ids;
    results.test_labels = settings.test_labels;
    results.reference_values = reference_values;
    results.est_mcmc = data.est_mcmc;
    results.est_sobol = data.est_sobol;
    results.mean_mcmc = nan(N_count, test_count);
    results.mean_sobol = nan(N_count, test_count);
    results.var_mcmc = nan(N_count, test_count);
    results.var_sobol = nan(N_count, test_count);
    results.vrf_sobol = nan(N_count, test_count);
    results.acceptance_mcmc = data.accept_mcmc;
    results.acceptance_sobol = data.accept_sobol;
    results.mean_acceptance_mcmc = mean(data.accept_mcmc, 1).';
    results.mean_acceptance_sobol = mean(data.accept_sobol, 1).';
    results.bias_mcmc = nan(N_count, test_count);
    results.bias_sobol = nan(N_count, test_count);
    results.empirical_mse_mcmc = nan(N_count, test_count);
    results.empirical_mse_sobol = nan(N_count, test_count);

    empty_fit = regression_slope([], [], slope_tail_count);
    results.slopes = struct( ...
        'tail_count_requested', slope_tail_count, ...
        'tail_count_used', min(slope_tail_count, N_count), ...
        'mcmc', repmat(empty_fit, 1, test_count), ...
        'sobol', repmat(empty_fit, 1, test_count));


    %% Summarize each test across replications and fit variance slopes
    for test_index = 1:test_count
        iid = data.est_mcmc(:, :, test_index);
        array = data.est_sobol(:, :, test_index);
        results.mean_mcmc(:, test_index) = mean(iid, 1).';
        results.mean_sobol(:, test_index) = mean(array, 1).';
        results.var_mcmc(:, test_index) = var(iid, 0, 1).';
        results.var_sobol(:, test_index) = var(array, 0, 1).';
        results.vrf_sobol(:, test_index) = safe_ratio( ...
            results.var_mcmc(:, test_index), ...
            results.var_sobol(:, test_index));
        results.slopes.mcmc(test_index) = regression_slope( ...
            N_vec, results.var_mcmc(:, test_index), slope_tail_count);
        results.slopes.sobol(test_index) = regression_slope( ...
            N_vec, results.var_sobol(:, test_index), slope_tail_count);

        reference = reference_values(test_index);
        if isfinite(reference)
            results.bias_mcmc(:, test_index) = ...
                results.mean_mcmc(:, test_index) - reference;
            results.bias_sobol(:, test_index) = ...
                results.mean_sobol(:, test_index) - reference;
            results.empirical_mse_mcmc(:, test_index) = ...
                mean((iid - reference).^2, 1).';
            results.empirical_mse_sobol(:, test_index) = ...
                mean((array - reference).^2, 1).';
        end
    end


    %% Bootstrap paired tests, independent methods and independent N blocks
    B=1000;if isfield(settings,'bootstrap_count'),B=settings.bootstrap_count;end
    bs=bootstrap_macro_statistics(cat(4,data.est_mcmc,data.est_sobol),N_vec, ...
        'Bootstrap',B,'Windows',slope_tail_count);
    results.bootstrap=bs;
    results.vrf_ci_low=nan(N_count,test_count);results.vrf_ci_high=nan(N_count,test_count);
    results.vrf_ci_low(end,:)=bs.vrfLow;results.vrf_ci_high(end,:)=bs.vrfHigh;
    for j=1:test_count
        for m=1:2
            names={'mcmc','sobol'};name=names{m};
            results.slopes.(name)(j).slope=bs.slope(j,m);
            results.slopes.(name)(j).bootstrap_low=bs.slopeLow(j,m);
            results.slopes.(name)(j).bootstrap_high=bs.slopeHigh(j,m);
            results.slopes.(name)(j).bootstrap_valid=bs.slopeValid(j,m);
            results.slopes.(name)(j).bootstrap_status=bs.slopeStatus(j,m);
        end
    end

    %% Attach settings and the established random-stream provenance
    results.settings = settings;
    results.reproducibility = struct( ...
        'schema_version', 3, 'generator', mfilename, ...
        'matlab_version', version, 'computer', computer, ...
        'stream_type', 'mrg32k3a', ...
        'statistical_method_offsets', [11, 37], ...
        'substream_formula', ...
            'm * 1e9 + m * 1e6 + replication');
end

%% Local helper: safe_ratio
% Divide valid variance pairs; retain Inf for positive numerator over zero denominator.
function values = safe_ratio(numerator, denominator)
    values = nan(size(numerator));
    regular = isfinite(numerator) & numerator >= 0 & ...
        isfinite(denominator) & denominator > 0;
    values(regular) = numerator(regular) ./ denominator(regular);
    values(denominator == 0 & numerator > 0) = Inf;
end

%% Local helper: regression_slope
% Adapt the common log-log fit fields to the historical result schema.
function fit = regression_slope(N_vec, variances, tail_count)
    generic = array_rqmc_fit_loglog_slope(N_vec, variances, tail_count);
    fit = struct('slope', generic.slope, ...
        'intercept', generic.intercept, ...
        'standard_error', generic.standard_error, ...
        'r_squared', generic.r_squared, ...
        'n_points', generic.n_points, 'N_used', generic.x_used);
end
