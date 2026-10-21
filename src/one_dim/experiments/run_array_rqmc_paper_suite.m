function manifest = run_array_rqmc_paper_suite(varargin)
%RUN_ARRAY_RQMC_PAPER_SUITE Run joint statistics and separate timing.
%
% Returns one manifest row per target/test pair and writes individual MAT files.
% Tests share statistical trajectories within a target; timing benchmarks each
% test separately. Completed target checkpoints are removed after result export.


    %% Resolve catalog selections, output paths and timing grid
    opt = array_rqmc_parse_suite_options(varargin{:});
    catalog = array_rqmc_case_catalog();
    targets = array_rqmc_select_catalog_entries( ...
        catalog.targets, opt.TargetIDs, 'target');
    tests = array_rqmc_select_catalog_entries( ...
        catalog.tests, opt.TestFunctionIDs, 'test function');
    output_dir = char(opt.OutputDir);
    if ~isfolder(output_dir)
        [ok, message] = mkdir(output_dir);
        if ~ok
            error('array_rqmc:CreateDirectory', ...
                'Could not create %s: %s', output_dir, message);
        end
    end

    m_vec = sort(unique(reshape(double(opt.mVec), 1, [])));
    if isempty(opt.TimingM)
        timing_m = max(m_vec);
    else
        timing_m = double(opt.TimingM);
    end
    if opt.RunTiming && ~ismember(timing_m, m_vec)
        error('array_rqmc:TimingMNotInGrid', ...
            'TimingM must belong to mVec.');
    end


    %% Allocate the manifest and record run settings
    manifest = repmat(manifest_template(), numel(targets) * numel(tests), 1);
    test_ids = {tests.id};
    test_labels = {tests.label};
    run_settings = opt;
    run_settings.resolved_timing_m = timing_m;
    requested_manifest = fullfile(output_dir, 'array_rqmc_suite_manifest.mat');
    run_settings.manifest_file = ...
        array_rqmc_unique_output_path(requested_manifest);

    fprintf('============================================================\n');
    fprintf('Array-RQMC suite: %d targets, %d tests, output %s\n', ...
        numel(targets), numel(tests), output_dir);
    fprintf('Joint statistics: yes; separate timing: %d\n', opt.RunTiming);
    fprintf('============================================================\n');


    suiteTimer=tic;
    %% Run joint statistics for each selected target
    for target_index = 1:numel(targets)
        targetTimer=tic;
        target = targets(target_index);
        [test_functions, ~, ~, reference_values] = ...
            array_rqmc_resolve_tests(tests, target, double(opt.InitialState));
        target_status = 'completed';
        target_message = 'Completed successfully.';
        timing_status = 'not requested';
        combined_results = [];
        target_timing = [];
        checkpoint_file = fullfile(output_dir, sprintf( ...
            'array_rqmc_checkpoint__target-%s.mat', target.id));

        fprintf('\n[%d/%d] target=%s; tests=%s\n', ...
            target_index, numel(targets), target.id, strjoin(test_ids, ', '));
        try
            combined_results = compare_array_rqmc_gaussian_rwm( ...
                'mVec', m_vec, 'T0', opt.T0, 'T', opt.T, 'R', opt.R, ...
                'sigma', opt.sigma, 'LogPi', target.log_pi, ...
                'TestFunctions', test_functions, ...
                'TestFunctionIDs', test_ids, ...
                'TestFunctionLabels', test_labels, ...
                'InitialState', opt.InitialState, ...
                'ReferenceValues', reference_values, ...
                'MasterSeed', opt.MasterSeed, ...
                'UseParallel', opt.UseParallel, ...
                'AutoStartPool', opt.AutoStartPool, ...
                'SlopeTailCount', opt.SlopeTailCount, ...
                'Bootstrap', opt.Bootstrap, ...
                'CheckpointFile', checkpoint_file, ...
                'ResumeFromCheckpoint', true);
            combined_results.settings = add_target_metadata( ...
                combined_results.settings, target, output_dir);
        catch ME
            target_status = 'failed';
            target_message = getReport(ME, 'extended', 'hyperlinks', 'off');
            warning('array_rqmc:TargetFailed', ...
                'Target %s failed: %s', target.id, ME.message);
            if ~opt.ContinueOnError
                rethrow(ME);
            end
        end


        %% Benchmark complete single-test runs separately
        if strcmp(target_status, 'completed') && opt.RunTiming
            try
                target_timing = benchmark_array_rqmc_single_test_times( ...
                    'm', timing_m, 'T0', opt.T0, 'T', opt.T, ...
                    'R', opt.TimingR, 'sigma', opt.sigma, ...
                    'LogPi', target.log_pi, ...
                    'TestFunctions', test_functions, ...
                    'TestFunctionIDs', test_ids, ...
                    'TestFunctionLabels', test_labels, ...
                    'InitialState', opt.InitialState, ...
                    'MasterSeed', opt.MasterSeed, ...
                    'WarmupM', opt.TimingWarmupM);
                timing_status = 'completed';
            catch ME
                timing_status = ['failed: ', ME.message];
                warning('array_rqmc:TimingFailed', ...
                    'Timing for %s failed: %s', target.id, ME.message);
                if ~opt.ContinueOnError
                    rethrow(ME);
                end
            end
        end


        %% Export individual test results and update the manifest
        for test_index = 1:numel(tests)
            case_index = (target_index - 1) * numel(tests) + test_index;
            test = tests(test_index);
            requested_file = array_rqmc_auto_filename( ...
                output_dir, target.id, test.id, m_vec, opt.T0, opt.T, ...
                opt.R, opt.sigma, opt.InitialState, opt.MasterSeed, test_ids);
            result_file = array_rqmc_unique_output_path(requested_file);
            manifest(case_index) = build_manifest_row( ...
                case_index, target, test, result_file, ...
                reference_values(test_index), target_status, ...
                timing_status, target_message, opt, timing_m);

            if strcmp(target_status, 'completed')
                results = select_array_rqmc_test_result( ...
                    combined_results, test_index);
                results.settings = add_test_metadata( ...
                    results.settings, target, test, output_dir);
                if strcmp(timing_status, 'completed')
                    results = array_rqmc_attach_single_test_timing( ...
                        results, target_timing, test_index);
                else
                    results.settings.timing_status = timing_status;
                end
                save(result_file, 'results', '-v7');
                fprintf('Saved %s/%s: %s\n', ...
                    target.id, test.id, result_file);
            end
        end
        if strcmp(target_status, 'completed')
            array_rqmc_checkpoint('delete', checkpoint_file, []);
        end
        array_rqmc_save_manifest(output_dir, manifest, run_settings);
        fprintf('Target total elapsed time: %.2f s\n',toc(targetTimer));
    end
    fprintf('Total elapsed time (suite): %.2f s\n',toc(suiteTimer));
end

%% Local helper: manifest_template
% Define the fixed manifest schema, including timing and failure status.
function row = manifest_template()
    row = struct('case_index', 0, 'target_id', '', 'test_id', '', ...
        'target_label', '', 'test_label', '', 'result_file', '', ...
        'reference_value', NaN, 'cdf_probability', NaN, ...
        'cdf_threshold', NaN, 'timing_N', NaN, 'timing_R', NaN, ...
        'status', '', 'timing_status', '', 'message', '');
end

%% Local helper: build_manifest_row
% Record a target/test artifact path, reference metadata and execution status.
function row = build_manifest_row(index, target, test, result_file, ...
        reference, status, timing_status, message, opt, timing_m)
    row = manifest_template();
    row.case_index = index;
    row.target_id = target.id;
    row.test_id = test.id;
    row.target_label = target.label;
    row.test_label = test.label;
    row.result_file = result_file;
    row.reference_value = reference;
    row.status = status;
    row.timing_status = timing_status;
    row.message = message;
    if opt.RunTiming
        row.timing_N = 2^timing_m;
        row.timing_R = double(opt.TimingR);
    end
    if strcmp(test.id, 'cdf')
        row.cdf_probability = target.cdf_probability;
        row.cdf_threshold = target.cdf_threshold;
    end
end

%% Local helper: add_target_metadata
% Attach target identity and CDF metadata to the saved experiment settings.
function settings = add_target_metadata(settings, target, output_dir)
    settings.target_id = target.id;
    settings.target_label = target.label;
    settings.target_latex_label = target.latex_label;
    settings.cdf_probability = target.cdf_probability;
    settings.cdf_threshold = target.cdf_threshold;
    settings.suite_output_dir = output_dir;
end

%% Local helper: add_test_metadata
% Extend target settings with the selected single-test identity.
function settings = add_test_metadata(settings, target, test, output_dir)
    settings = add_target_metadata(settings, target, output_dir);
    settings.test_function_id = test.id;
    settings.test_function_label = test.label;
end
