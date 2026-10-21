function [figure_handle, output_stem] = ...
        plot_array_rqmc_chain_length_scaling_results( ...
        results, figure_stem, varargin)
%PLOT_ARRAY_RQMC_CHAIN_LENGTH_SCALING_RESULTS Plot empirical multi-N scaling.
%
% Panel 1 shows Var(mu_hat_{N,T}) versus L=T-T0 for each N.
% Panel 2 shows N^(4/3)*Var(mu_hat_{N,T}) versus L. Similar scaled curves
% and slopes near -1 are empirical support only, not a theorem.

    if nargin < 2 || isempty(figure_stem)
        root=fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
        figure_stem=fullfile(root,'figures',['variance_' results.target_id '_T']);
    end
    figure_stem = char(figure_stem);

    ip = inputParser;
    addParameter(ip, 'Visible', 'on', ...
        @(x) any(strcmpi(char(x), {'on','off'})));
    addParameter(ip, 'Save', true, @(x) islogical(x) && isscalar(x));
    parse(ip, varargin{:});
    opt = ip.Results;

    validate_results(results);
    L = double(results.effective_T(:));
    variance_array = double(results.variance_array);
    scaled_variance = double(results.scaled_variance_array);
    m_vec = double(results.m_vec(:));
    colors = lines(numel(m_vec));

    figure_handle = figure('Color', 'w', ...
        'Position', [70, 80, 1180, 500], ...
        'Visible', char(opt.Visible));
    layout = tiledlayout(figure_handle, 1, 2, ...
        'TileSpacing', 'compact', 'Padding', 'compact');

    ax_variance = nexttile(layout);
    set(ax_variance, 'XScale', 'log', 'YScale', 'log');
    draw_curves(ax_variance, L, variance_array, m_vec, ...
        results.slope_array, colors);
    draw_inverse_length_reference(ax_variance, L, variance_array);
    grid(ax_variance, 'on'); box(ax_variance, 'on');
    ax_variance.XMinorGrid = 'on';
    ax_variance.YMinorGrid = 'on';
    xlabel(ax_variance, '$L=T-T_0$', 'Interpreter', 'latex','FontSize',12);
    ylabel(ax_variance, ...
        '$\mathrm{Var}(\widehat{\mu}_{N,T})$', 'Interpreter', 'latex','FontSize',12);
    title(ax_variance, 'Array-RQMC variance', 'Interpreter', 'latex','FontSize',12);
    legend(ax_variance, 'Location', 'best', 'Interpreter', 'latex','FontSize',12);

    ax_scaled = nexttile(layout);
    set(ax_scaled, 'XScale', 'log', 'YScale', 'log');
    draw_curves(ax_scaled, L, scaled_variance, m_vec, ...
        results.slope_array, colors);
    draw_inverse_length_reference(ax_scaled, L, scaled_variance);
    grid(ax_scaled, 'on'); box(ax_scaled, 'on');
    ax_scaled.XMinorGrid = 'on';
    ax_scaled.YMinorGrid = 'on';
    xlabel(ax_scaled, '$L=T-T_0$', 'Interpreter', 'latex','FontSize',12);
    ylabel(ax_scaled, ...
        '$N^{4/3}\,\mathrm{Var}(\widehat{\mu}_{N,T})$', ...
        'Interpreter', 'latex','FontSize',12);
    title(ax_scaled, 'Scaled variance', 'Interpreter', 'latex','FontSize',12);
    legend(ax_scaled, 'Location', 'best', 'Interpreter', 'latex','FontSize',12);
    output_stem = figure_stem;
    if opt.Save
        ensure_parent_directory(output_stem);
        savefig(figure_handle, [output_stem, '.fig']);
        exportgraphics(figure_handle, [output_stem, '.pdf'], ...
            'ContentType', 'vector');
    end
end

function ensure_parent_directory(figure_stem)
    output_dir = fileparts(figure_stem);
    if ~isempty(output_dir) && ~isfolder(output_dir)
        [ok, message] = mkdir(output_dir);
        if ~ok
            error('array_rqmc:CreateFigureDirectory', ...
                'Could not create %s: %s', output_dir, message);
        end
    end
end

function draw_curves(ax, L, values, m_vec, slopes, colors)
    hold(ax, 'on');
    for N_index = 1:numel(m_vec)
        loglog(ax, L, values(:, N_index), '-o', ...
            'Color', colors(N_index, :), ...
            'LineWidth', 1.8, 'MarkerSize', 5.5, ...
            'DisplayName', sprintf( ...
                '$N=2^{%d}$ (slope %.2f)', ...
                m_vec(N_index), slopes(N_index)));
    end
end

function draw_inverse_length_reference(ax, L, values)
    valid_row = find(any(isfinite(values) & values > 0, 2), 1, 'first');
    if isempty(valid_row)
        return;
    end
    anchors = values(valid_row, :);
    anchors = anchors(isfinite(anchors) & anchors > 0);
    if isempty(anchors)
        return;
    end
    reference = median(anchors) .* (L ./ L(valid_row)).^(-1);
    loglog(ax, L, reference, '--k', 'LineWidth', 1.2, ...
        'DisplayName', '$L^{-1}$ reference');
end

function validate_results(results)
    required = {'effective_T','N_vec','m_vec','variance_array', ...
        'scaled_variance_array','slope_array','target_label','test_label'};
    if ~isstruct(results) || ~all(isfield(results, required))
        error('array_rqmc:InvalidScalingResults', ...
            'Scaling results are missing required fields.');
    end
    expected_size = [numel(results.effective_T), numel(results.N_vec)];
    if ~isequal(size(results.variance_array), expected_size) || ...
            ~isequal(size(results.scaled_variance_array), expected_size)
        error('array_rqmc:ScalingResultSize', ...
            'Variance matrices must have size numel(T) by numel(N).');
    end
end
