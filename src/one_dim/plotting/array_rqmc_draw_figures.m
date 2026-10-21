function [figures, rows, problems] = array_rqmc_draw_figures( ...
        file_index, target_ids, test_ids, opt, output_dir)
%ARRAY_RQMC_DRAW_FIGURES Dispatch supported suite layouts.
%
% Creates per-target, combined or paired-target layouts from a saved-file index.
% Returns figure handles, slope rows and missing-data messages collected by panels.
% SaveFigures writes FIG/PDF; CloseAfterSave controls figure cleanup independently.


    %% Allocate panel metadata and dispatch the selected layout
    rows = cell(0, 13);
    problems = cell(0, 1);
    switch lower(char(opt.FigureMode))
        case 'per-target'
            figures = gobjects(numel(target_ids), 1);
            for target_index = 1:numel(target_ids)
                target_id = target_ids{target_index};
                fig = figure('Color', 'w', 'Position', [60,100, ...
                    max(620, round(1520*numel(test_ids)/3)),480], ...
                    'Visible', char(opt.FigureVisible), ...
                    'Name', ['Variance: ', target_id]);
                figures(target_index) = fig;
                layout = tiledlayout(fig, 1, numel(test_ids), ...
                    'TileSpacing', 'compact', 'Padding', 'compact');
                for test_index = 1:numel(test_ids)
                    [rows, problems] = array_rqmc_draw_result_panel( ...
                        nexttile(layout), file_index, target_id, ...
                        test_ids{test_index}, opt, rows, problems, ...
                        test_index == 1, 'per-target');
                end
                finish_figure(fig, fullfile(output_dir, ...
                    ['variance_', target_id]), opt);
            end
        case 'combined'
            fig = figure('Color', 'w', 'Position', [40,40, ...
                max(1100,500*numel(test_ids)), ...
                max(700,280*numel(target_ids))], ...
                'Visible', char(opt.FigureVisible), ...
                'Name', 'Array-RQMC variance suite');
            figures = fig;
            layout = tiledlayout(fig, numel(target_ids), numel(test_ids), ...
                'TileSpacing', 'compact', 'Padding', 'compact');
            for target_index = 1:numel(target_ids)
                for test_index = 1:numel(test_ids)
                    [rows, problems] = array_rqmc_draw_result_panel( ...
                        nexttile(layout), file_index, ...
                        target_ids{target_index}, test_ids{test_index}, ...
                        opt, rows, problems, test_index == 1, 'combined');
                end
            end
            finish_figure(fig, fullfile(output_dir, 'variance_suite'), opt);
        case 'paper-pairs'
            [figures, rows, problems] = draw_paper_pairs( ...
                file_index, target_ids, test_ids, opt, output_dir, ...
                rows, problems);
        otherwise
            error('array_rqmc:FigureMode', ...
                'Unsupported FigureMode: %s', char(opt.FigureMode));
    end
end

%% Local helper: draw_paper_pairs
% Lay out at most two targets per figure and collect their panel diagnostics.
function [figures, rows, problems] = draw_paper_pairs( ...
        file_index, target_ids, test_ids, opt, output_dir, rows, problems)
    figure_count = ceil(numel(target_ids) / 2);
    figures = gobjects(figure_count, 1);
    for figure_index = 1:figure_count
        first = (figure_index - 1) * 2 + 1;
        local_targets = target_ids(first:min(first + 1, numel(target_ids)));
        fig = figure('Color', 'w', 'Position', [40,40, ...
            max(1100,500*numel(test_ids)),420*numel(local_targets)], ...
            'Visible', char(opt.FigureVisible));
        figures(figure_index) = fig;
        layout = tiledlayout(fig, numel(local_targets), numel(test_ids), ...
            'TileSpacing', 'compact', 'Padding', 'loose');
        for row_index = 1:numel(local_targets)
            for test_index = 1:numel(test_ids)
                [rows, problems] = array_rqmc_draw_result_panel( ...
                    nexttile(layout), file_index, local_targets{row_index}, ...
                    test_ids{test_index}, opt, rows, problems, ...
                    test_index == 1, 'paper-pair');
            end
        end
        finish_figure(fig, fullfile(output_dir, sprintf( ...
            'variance_suite_part%d', figure_index)), opt);
    end
end

%% Local helper: finish_figure
% Export requested FIG/PDF files and close the figure only when requested.
function finish_figure(fig, stem, opt)
    if opt.SaveFigures
        output_dir = fileparts(stem);
        if ~isfolder(output_dir)
            [ok, message] = mkdir(output_dir);
            if ~ok
                error('array_rqmc:CreateFigureDirectory', ...
                    'Could not create %s: %s', output_dir, message);
            end
        end
        savefig(fig, [stem, '.fig']);
        exportgraphics(fig, [stem, '.pdf'], 'ContentType', 'vector');
    end
    if opt.CloseAfterSave
        close(fig);
    end
end
