function results = summarize_array_rqmc_chain_length(raw, context)
%SUMMARIZE_ARRAY_RQMC_CHAIN_LENGTH Compute reusable chain-length summaries.

    results = context;
    results.effective_T = context.T_vec - context.T0;
    results.est_mcmc_all = raw.est_mcmc_all;
    results.est_sobol_all = raw.est_sobol_all;
    results.acc_mcmc_all = raw.acc_mcmc_all;
    results.acc_sobol_all = raw.acc_sobol_all;
    results.elapsed_sec = raw.elapsed_sec;
    results.completed_idx = raw.completed_idx;

    results.mean_mcmc = mean(raw.est_mcmc_all, 1).';
    results.mean_sobol = mean(raw.est_sobol_all, 1).';
    results.var_mcmc = var(raw.est_mcmc_all, 0, 1).';
    results.var_sobol = var(raw.est_sobol_all, 0, 1).';
    results.vrf_vec = results.var_mcmc ./ results.var_sobol;
    results.mean_acc_mcmc = mean(raw.acc_mcmc_all, 1).';
    results.mean_acc_sobol = mean(raw.acc_sobol_all, 1).';

    fit_iid = array_rqmc_fit_loglog_slope( ...
        results.effective_T, results.var_mcmc, numel(context.T_vec));
    fit_array = array_rqmc_fit_loglog_slope( ...
        results.effective_T, results.var_sobol, numel(context.T_vec));
    results.fit_iid = [fit_iid.intercept; fit_iid.slope];
    results.fit_arqmc = [fit_array.intercept; fit_array.slope];
    results.slope_iid = fit_iid.slope;
    results.slope_arqmc = fit_array.slope;
    results.slope_details_iid = fit_iid;
    results.slope_details_arqmc = fit_array;

    reference = context.stationary_reference;
    results.bias_mcmc = nan(size(results.mean_mcmc));
    results.bias_sobol = nan(size(results.mean_sobol));
    results.empirical_mse_mcmc = nan(size(results.mean_mcmc));
    results.empirical_mse_sobol = nan(size(results.mean_sobol));
    if isfinite(reference)
        results.bias_mcmc = results.mean_mcmc - reference;
        results.bias_sobol = results.mean_sobol - reference;
        results.empirical_mse_mcmc = ...
            mean((raw.est_mcmc_all - reference).^2, 1).';
        results.empirical_mse_sobol = ...
            mean((raw.est_sobol_all - reference).^2, 1).';
    end

    alpha = 1 - context.confidence_level;
    critical = tinv(1 - alpha / 2, context.R - 1);
    iid_half_width = critical .* sqrt(results.var_mcmc ./ context.R);
    array_half_width = critical .* sqrt(results.var_sobol ./ context.R);
    results.ci_mcmc = [results.mean_mcmc - iid_half_width, ...
        results.mean_mcmc + iid_half_width];
    results.ci_sobol = [results.mean_sobol - array_half_width, ...
        results.mean_sobol + array_half_width];

    results.summary_table = table(context.T_vec(:), ...
        results.effective_T(:), results.mean_mcmc, results.var_mcmc, ...
        results.mean_sobol, results.var_sobol, results.vrf_vec, ...
        results.bias_mcmc, results.bias_sobol, ...
        results.empirical_mse_mcmc, results.empirical_mse_sobol, ...
        results.ci_mcmc(:,1), results.ci_mcmc(:,2), ...
        results.ci_sobol(:,1), results.ci_sobol(:,2), ...
        results.mean_acc_mcmc, results.mean_acc_sobol, raw.elapsed_sec(:), ...
        'VariableNames', {'T','EffectiveLength','Mean_IID','Var_IID', ...
        'Mean_Array','Var_Array','VRF','Bias_IID','Bias_Array', ...
        'MSE_IID','MSE_Array','CI_Lower_IID','CI_Upper_IID', ...
        'CI_Lower_Array','CI_Upper_Array','MeanAcceptance_IID', ...
        'MeanAcceptance_Array','ElapsedSeconds'});
    results.saved_at = datetime('now');
    results.matlab_version = version;
    results.computer = computer;
    results.reproducibility = struct('schema_version', 3, ...
        'stream_type', 'mrg32k3a', 'method_offsets', [11, 37], ...
        'substream_formula', ...
            'm * 1e9 + T index * 1e6 + replication');
end
