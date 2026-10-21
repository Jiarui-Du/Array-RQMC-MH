function write_array_rqmc_latex_table( ...
        file_path, summary, target_ids, test_ids, opt)
%WRITE_ARRAY_RQMC_LATEX_TABLE Write the paper VRF/TaE/slope table.
%
% Writes one row per target with VRF, TAE and slope columns for each test.
% Nonfinite or missing values are rendered as dashes; this function formats
% provided summaries and does not recompute experimental statistics.


    %% Open the TeX destination and write the caption and column headings
    fid = fopen(file_path, 'w');
    if fid < 0
        error('array_rqmc:OpenTeX', 'Could not open %s.', file_path);
    end
    cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid, '\\begin{table}[t]\n\\centering\n');
    fprintf(fid, ['\\caption{Variance-reduction factors and median-time-', ...
        'adjusted efficiency at $N=2^{%d}$, with Array-RQMC log--log ', ...
        'variance slopes fitted over %s. Timing uses separate serial ', ...
        'single-function runs.}\n'], round(log2(double(opt.TargetN))), ...
        fit_range_text(summary, opt.TailCount));
    fprintf(fid, '\\label{tab:all-vrf-tae-slopes}\n');
    fprintf(fid, '\\resizebox{\\textwidth}{!}{%%\n');
    fprintf(fid, '\\begin{tabular}{%s}\n\\toprule\n', ...
        ['l', repmat('rrr', 1, numel(test_ids))]);
    fprintf(fid, 'Target');
    for test_index = 1:numel(test_ids)
        fprintf(fid, ' & \\multicolumn{3}{c}{%s}', ...
            array_rqmc_export_label('test_tex', test_ids{test_index}));
    end
    fprintf(fid, ' \\\\\n');
    for test_index = 1:numel(test_ids)
        first = 2 + 3 * (test_index - 1);
        fprintf(fid, '\\cmidrule(lr){%d-%d}', first, first + 2);
    end
    fprintf(fid, '\n');
    for test_index = 1:numel(test_ids)
        fprintf(fid, ' & VRF & TaE & Slope');
    end
    fprintf(fid, ' \\\\\n\\midrule\n');


    %% Render target/test rows with missing-value placeholders
    for target_index = 1:numel(target_ids)
        target_id = target_ids{target_index};
        fprintf(fid, '%s', array_rqmc_export_label('target_tex', target_id));
        for test_index = 1:numel(test_ids)
            match = find(strcmp(summary.TargetID, target_id) & ...
                strcmp(summary.TestFunctionID, test_ids{test_index}), 1);
            if isempty(match)
                fprintf(fid, ' & -- & -- & --');
            else
                fprintf(fid, ' & %s & %s & %s', ...
                    integer_or_dash(summary.VRF(match)), ...
                    integer_or_dash(summary.TaE(match)), ...
                    slope_or_dash(summary.ArraySlope(match)));
            end
        end
        fprintf(fid, ' \\\\\n');
    end
    fprintf(fid, '\\bottomrule\n\\end{tabular}%%\n}\n\\end{table}\n');
end

%% Local helper: fit_range_text
% Describe a shared power-of-two fit interval, or use the requested tail-count label.
function text = fit_range_text(summary, tail_count)
    valid = isfinite(summary.NMinSlope) & summary.NMinSlope > 0 & ...
        isfinite(summary.NMaxSlope) & summary.NMaxSlope > 0;
    if any(valid)
        N_min = summary.NMinSlope(valid);
        N_max = summary.NMaxSlope(valid);
        if all(N_min == N_min(1)) && all(N_max == N_max(1))
            m_min = log2(N_min(1)); m_max = log2(N_max(1));
            if abs(m_min-round(m_min)) < 1e-12 && ...
                    abs(m_max-round(m_max)) < 1e-12
                text = sprintf('$N=2^{%d},\\ldots,2^{%d}$', ...
                    round(m_min), round(m_max));
                return;
            end
        end
    end
    text = sprintf('the %d largest particle counts', round(tail_count));
end

%% Local helper: integer_or_dash
% Format finite efficiency values as integers; mark unavailable values with dashes.
function text = integer_or_dash(value)
    if isfinite(value), text = sprintf('%.0f', value); else, text = '--'; end
end
%% Local helper: slope_or_dash
% Format finite slopes to two decimals; mark unavailable values with dashes.
function text = slope_or_dash(value)
    if isfinite(value), text = sprintf('%.2f', value); else, text = '--'; end
end
