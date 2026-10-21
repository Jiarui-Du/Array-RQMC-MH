function artifacts = save_array_rqmc_chain_length_scaling_results(results)
%SAVE_ARRAY_RQMC_CHAIN_LENGTH_SCALING_RESULTS Save aggregate MAT and CSV.

    requested_mat = results.result_file;
    mat_file = array_rqmc_unique_output_path(requested_mat);
    [folder, stem] = fileparts(mat_file);
    if ~isfolder(folder)
        [ok, message] = mkdir(folder);
        if ~ok
            error('array_rqmc:CreateDirectory', ...
                'Could not create %s: %s', folder, message);
        end
    end

    results.result_file = mat_file;
    save(mat_file, 'results', '-v7.3');
    artifacts = struct( ...
        'mat', mat_file, ...
        'summary_csv', fullfile(folder, [stem, '_summary.csv']), ...
        'slopes_csv', fullfile(folder, [stem, '_slopes.csv']), ...
        'figure_stem', fullfile(folder, stem));
    writetable(results.summary_table, artifacts.summary_csv);
    writetable(results.slope_table, artifacts.slopes_csv);
end
