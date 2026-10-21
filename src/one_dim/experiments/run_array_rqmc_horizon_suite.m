function [results,artifacts]=run_array_rqmc_horizon_suite(options)
% Fixed-N horizon diagnostic using the same kernels as the one-dimensional suite.
% No independent runtime benchmark: elapsed times describe statistical work only.
    arguments
        options.TargetID (1,:) char = 'gaussian_mixture3'
        options.TestFunctionID (1,:) char = 'cdf'
        options.mVec double = [12 14 16]
        options.T0 (1,1) double {mustBeInteger,mustBeNonnegative} = 128
        options.TVec double = 2.^(8:16)
        options.R (1,1) double {mustBeInteger,mustBePositive} = 50
        options.sigma (1,1) double {mustBePositive} = 2.38
        options.InitialState (1,1) double = 0
        options.MasterSeed (1,1) double {mustBeInteger,mustBeNonnegative} = 2026
        options.UseParallel (1,1) logical = true
        options.AutoStartPool (1,1) logical = true
        options.SlopeTailCount (1,1) double {mustBeInteger,mustBePositive} = 5
        options.Bootstrap (1,1) double {mustBeInteger,mustBeNonnegative} = 1000
        options.OutputDir (1,:) char = ''
        options.ResumeDir (1,:) char = ''
    end
    %% Validate ordered grids and reserve a separate run directory
    assert(isvector(options.mVec) && ~isempty(options.mVec) && all(isfinite(options.mVec)) ...
        && all(options.mVec>=1 & options.mVec==fix(options.mVec)) && all(diff(options.mVec)>0));
    assert(isvector(options.TVec) && numel(options.TVec)>=2 && all(isfinite(options.TVec)) ...
        && all(options.TVec==fix(options.TVec)) && all(diff(options.TVec)>0) && all(options.TVec>options.T0));
    assert(options.R>=2 && options.SlopeTailCount>=2 && isfinite(options.sigma) && isfinite(options.InitialState));
    root=fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
    folder=options.OutputDir;if isempty(folder),folder=fullfile(root,'runs','one_dim_horizon');end
    if isempty(options.ResumeDir)
        folder=fullfile(folder,['run_' char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'))]);
        assert(~isfolder(folder),'Run directory already exists; retry.');mkdir(folder);
    else
        folder=options.ResumeDir;assert(isfolder(folder),'Resume directory does not exist.');
    end
    fprintf('Run directory: %s\n',folder);
    %% Retain the original N/T/replication substream mapping
    results=run_array_rqmc_chain_length_scaling_experiment( ...
        'TargetID',options.TargetID,'TestFunctionID',options.TestFunctionID, ...
        'NVec',2.^options.mVec,'T0',options.T0,'TVec',options.TVec,'R',options.R, ...
        'sigma',options.sigma,'InitialState',options.InitialState,'MasterSeed',options.MasterSeed, ...
        'UseParallel',options.UseParallel,'AutoStartPool',options.AutoStartPool, ...
        'SingleNOutputDir',fullfile(folder,'individual'),'OutputDir',folder, ...
        'SlopeTailCount',options.SlopeTailCount,'ResumeFromCheckpoint',true,'ReuseCompletedResults',true);
    %% Bootstrap independent horizons and methods, using the shared degeneracy rule
    results.bootstrap=cell(numel(results.N_vec),1);
    for k=1:numel(results.N_vec)
        single=results.individual_results{k};
        A=cat(4,single.est_mcmc_all,single.est_sobol_all);
        bs=bootstrap_macro_statistics(A,results.effective_T,'Bootstrap',options.Bootstrap, ...
            'Seed',options.MasterSeed+round(log2(results.N_vec(k)))*10000, ...
            'Windows',options.SlopeTailCount,'PairedN',false,'PairedMethods',false);
        bs.gridMeaning='effective horizon L=T-T0; N axis in shared helper denotes L';
        results.bootstrap{k}=bs;
        results.slope_table.ArraySlopeLow(k)=bs.slopeLow(1,2);
        results.slope_table.ArraySlopeHigh(k)=bs.slopeHigh(1,2);
        results.slope_table.VRFAtMaxL(k)=bs.vrf;
        results.slope_table.VRFLow(k)=bs.vrfLow;
        results.slope_table.VRFHigh(k)=bs.vrfHigh;
    end
    results.settings=options;
    artifacts=save_array_rqmc_chain_length_scaling_results(results);
    results.result_file=artifacts.mat;
    fprintf('Saved result: %s\n',artifacts.mat);
end
