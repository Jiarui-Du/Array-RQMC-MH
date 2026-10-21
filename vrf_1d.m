% Display one-dimensional VRF, standalone timing efficiency and slopes.
%
% Entry script: displays saved variance ratios, separate-timing efficiency and
% empirical slopes. RequiredR and TargetN select the intended archived run.

%% Resolve saved result paths
root=fileparts(mfilename('fullpath')); addpath(genpath(fullfile(root,'src')));
resultsDir=fullfile(root,'results','one_dim');

%% Display the summary; set nonempty OutputCSV/OutputTeX paths to export
summary=export_array_rqmc_vrf_tae_slopes_v2(resultsDir, ...
    'TargetN',2^18,'TailCount',5,'RequiredR',200,'Strict',true, ...
    'OutputCSV','', ...
    'OutputTeX','');
