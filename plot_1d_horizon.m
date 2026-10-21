% Read a completed horizon result and export the supplementary figure.
%% Resolve paths and choose a result
root=fileparts(mfilename('fullpath'));addpath(genpath(fullfile(root,'src')));
matFile=''; % Empty selects the newest complete mixture/CDF result.
saveFigures=true;showFigures=true;
%% Load and export; no experiments are run
[results,matFile]=load_array_rqmc_horizon_result(matFile);
visible='off';if showFigures,visible='on';end
[figureHandle,figureStem]=plot_array_rqmc_chain_length_scaling_results(results, ...
    fullfile(root,'figures',['variance_' results.target_id '_T']), ...
    'Visible',visible,'Save',saveFigures);
