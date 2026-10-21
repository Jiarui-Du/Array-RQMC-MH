% Four one-dimensional paper targets. Edit parameters here; running starts simulation.
%
% Entry script: runs all selected targets and writes a manifest plus one MAT file
% per target/test pair. T includes T0; retained observations are T0+1:T.
% The settings below are full experiment defaults; reduce them for smoke tests.

%% Resolve repository paths
root=fileparts(mfilename('fullpath')); addpath(genpath(fullfile(root,'src')));

%% Select targets, tests and simulation settings
targetIDs={'normal','double_well_e1','gaussian_mixture3','student_t3'};
testIDs={'cdf','x','x2'};
mVec=10:18; T0=8; T=128; R=200; sigma=2.38; x0=0; seed=2026;

%% Configure workers and separate timing
useParallel=true; numWorkers=4;
runTiming=true; timingR=20;
outputDir=fullfile(root,'results','one_dim');

%% Start the requested pool
if useParallel
    pool=gcp('nocreate');
    if isempty(pool),parpool('Processes',numWorkers);end
end

%% Run joint statistics and export individual test results
manifest=run_array_rqmc_suite('TargetIDs',targetIDs,'TestFunctionIDs',testIDs, ...
    'mVec',mVec,'T0',T0,'T',T,'R',R,'sigma',sigma,'InitialState',x0, ...
    'MasterSeed',seed,'UseParallel',useParallel,'AutoStartPool',false, ...
    'SlopeTailCount',5,'RunTiming',runTiming,'TimingM',max(mVec), ...
    'TimingR',timingR,'TimingWarmupM',12,'OutputDir',outputDir,'ContinueOnError',false);
