% Fixed-N, increasing-horizon mixture experiment. Running starts simulation.
% T includes burn-in; each N/T pair uses independent macro replications.

%% Resolve repository paths
root=fileparts(mfilename('fullpath'));addpath(genpath(fullfile(root,'src')));

%% Select target, test and simulation settings
targetID='gaussian_mixture3';testID='cdf';
mVec=[12 14 16];T0=128;TVec=2.^(8:16);
R=50;sigma=2.38;x0=0;seed=2026;slopeTailCount=5;

%% Configure workers and output
useParallel=true;numWorkers=4;
outputDir=fullfile(root,'runs','one_dim_horizon');
% Empty starts a fresh run. For resumption, specify its printed run directory.
resumeDir='';
if useParallel
    pool=gcp('nocreate');
    if isempty(pool),parpool('Processes',numWorkers);end
end

%% Run statistics and save MAT/CSV outputs
[results,artifacts]=run_array_rqmc_horizon_suite( ...
    'TargetID',targetID,'TestFunctionID',testID,'mVec',mVec, ...
    'T0',T0,'TVec',TVec,'R',R,'sigma',sigma,'InitialState',x0, ...
    'MasterSeed',seed,'UseParallel',useParallel,'AutoStartPool',false, ...
    'SlopeTailCount',slopeTailCount,'OutputDir',outputDir,'ResumeDir',resumeDir);
