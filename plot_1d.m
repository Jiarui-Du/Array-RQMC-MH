% Read saved one-dimensional results only; never rerun experiments.
%
% Entry script: selects saved MAT files and produces one figure per target.
% RequiredR filters the saved replication count; this script starts no simulation.

%% Resolve source and result folders
root=fileparts(mfilename('fullpath')); addpath(genpath(fullfile(root,'src')));
resultsDir=fullfile(root,'results','one_dim'); % Fresh results from run_1d.

%% Select panels and figure visibility
targetIDs={'normal','double_well_e1','gaussian_mixture3','student_t3'};
testIDs={'cdf','x','x2'};
saveFigures=true; visible='on';

%% Read results, draw panels and optionally save figures
[figures,slopes,index]=plot_array_rqmc_suite_results(resultsDir, ...
    'TargetIDs',targetIDs,'TestFunctionIDs',testIDs,'RequiredR',200,'TailCount',5, ...
    'FigureMode','per-target','SaveFigures',saveFigures,'FigureVisible',visible, ...
    'OutputDir',fullfile(root,'figures'), ...
    'MetadataOutputDir','');
