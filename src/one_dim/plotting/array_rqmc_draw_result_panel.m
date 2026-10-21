function [rows, problems] = array_rqmc_draw_result_panel( ...
        ax, file_index, target_id, test_id, opt, rows, problems, ...
        show_y_label, panel_mode)
%ARRAY_RQMC_DRAW_RESULT_PANEL Draw one target/test variance panel.
%
% Selects one saved case and draws positive paired variances plus reference lines.
% Appends two fitted-slope rows on success; missing data produce an annotated
% panel and a diagnostic row rather than a replacement simulation.


    %% Select and validate saved variance data
    [selected, selection_status] = array_rqmc_select_result_file( ...
        file_index, target_id, test_id, opt.RequiredR, []);
    if isempty(selected)
        [rows, problems] = missing_result(ax, rows, problems, ...
            target_id, test_id, selection_status, panel_mode);
        return;
    end

    file_path = selected.Path{1};
    try
        loaded = load(file_path, 'results');
        if ~isfield(loaded, 'results')
            error('Missing results struct.');
        end
        [N, variance_iid, variance_array, array_name] = ...
            array_rqmc_extract_variance_data(loaded.results, true);
    catch ME
        [rows, problems] = missing_result(ax, rows, problems, ...
            target_id, test_id, ME.message, panel_mode);
        return;
    end
    if isempty(N)
        [rows, problems] = missing_result(ax, rows, problems, ...
            target_id, test_id, 'No positive paired variances', panel_mode);
        return;
    end


    %% Draw empirical curves, tail fits and anchored reference guides
    fit_iid = array_rqmc_fit_loglog_slope(N, variance_iid, opt.TailCount);
    fit_array = array_rqmc_fit_loglog_slope(N, variance_array, opt.TailCount);
    loglog(ax, N, variance_iid, '-o', 'LineWidth', 1.8, ...
        'MarkerSize', 6, ...
        'DisplayName', sprintf('IID, slope %.3f', fit_iid.slope));
    hold(ax, 'on');
    loglog(ax, N, variance_array, '-s', 'LineWidth', 1.8, ...
        'MarkerSize', 6, ...
        'DisplayName', sprintf('%s, slope %.3f', ...
        array_name, fit_array.slope));
    loglog(ax, N, variance_iid(1) .* (N ./ N(1)).^(-1), '--', ...
        'LineWidth', 1.1, 'DisplayName', '$N^{-1}$');
    anchor = max(1, ceil(numel(N) / 2));
    loglog(ax, N, variance_array(anchor) .* ...
        (N ./ N(anchor)).^(-4/3), ':', 'LineWidth', 1.3, ...
        'DisplayName', '$N^{-4/3}$');
    grid(ax, 'on'); box(ax, 'on');
    xlabel(ax, '$N$', 'Interpreter', 'latex', 'FontSize', 16);
    if show_y_label
        ylabel(ax, 'Sample variance', 'FontSize', 16);
    end
    if strcmpi(panel_mode, 'combined')
        title(ax, sprintf('%s; %s', ...
            array_rqmc_plot_label('target', target_id), ...
            array_rqmc_plot_label('test', test_id)), ...
            'Interpreter', 'latex', 'FontSize', 11);
        font_size = 8;
    else
        title(ax, array_rqmc_plot_label('test', test_id), ...
            'Interpreter', 'latex', 'FontSize', 16);
        font_size = 10;
    end
    legend(ax, 'Location', 'southwest', ...
        'Interpreter', 'latex', 'FontSize', font_size);


    %% Append the plotted data provenance and fitted slopes
    R = selected.Replications(1);
    rows = append_fit_row(rows, target_id, test_id, 'IID-MCMC', ...
        fit_iid, file_path, selection_status, R);
    rows = append_fit_row(rows, target_id, test_id, array_name, ...
        fit_array, file_path, selection_status, R);
end

%% Local helper: append_fit_row
% Record fitted slopes together with their N range and selected source file.
function rows = append_fit_row(rows, target_id, test_id, method, ...
        fit, file_path, selection_status, R)
    if isempty(fit.x_used)
        N_min = NaN; N_max = NaN;
    else
        N_min = fit.x_used(1); N_max = fit.x_used(end);
    end
    rows(end+1, :) = {target_id, test_id, method, fit.slope, ...
        fit.standard_error, fit.r_squared, fit.n_points, N_min, N_max, ...
        file_path, selection_status, R, 'complete'};
end

%% Local helper: missing_result
% Render a missing-data panel and append explicit diagnostic metadata.
function [rows, problems] = missing_result( ...
        ax, rows, problems, target_id, test_id, message, panel_mode)
    axis(ax, 'off');
    if strcmpi(panel_mode, 'combined')
        title(ax, sprintf('%s; %s', ...
            array_rqmc_plot_label('target', target_id), ...
            array_rqmc_plot_label('test', test_id)), 'Interpreter', 'latex');
    else
        title(ax, array_rqmc_plot_label('test', test_id), ...
            'Interpreter', 'latex');
    end
    text(ax, 0.5, 0.5, message, 'HorizontalAlignment', 'center', ...
        'Interpreter', 'none');
    problems{end+1, 1} = sprintf('%s/%s: %s', ...
        target_id, test_id, message);
    rows(end+1, :) = {target_id, test_id, 'Unavailable', ...
        NaN, NaN, NaN, 0, NaN, NaN, '', message, NaN, 'missing'};
end
