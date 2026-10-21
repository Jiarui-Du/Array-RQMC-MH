% Read and display variance reduction; no simulation or plotting.
%
% Entry script: reads saved replications for the selected dataset and tests.
% VRF uses the matching IID kernel; TAE additionally accounts for elapsed time.

%% Resolve repository paths
driver=mfilename('fullpath');root=fileparts(driver);
addpath(genpath(fullfile(root,'src')));

%% Select dataset, dimensions and particle counts
dataset='ripley'; % 'ripley', 'pima', 'gaussian', or 'gaussian_mwg' ('mwg')
dimension=[2 4 8];      % For gaussian / gaussian_mwg; select available dimensions
exponents=18;      % Use 18 to show only N=2^18
tailCount=5; bootstrap=1000; % Slope fit uses the source grid, even if exponents=18
tests='projection';  % Paper results: projection for Gaussian, MWG and Bayesian experiments
% Display VRF and slope intervals, serial median seconds, and TAE (NaN if not timed).
matFile='';          % Empty checks newest to oldest under results/; reads the first complete matching dimensions/N.
csvFile='';          % Empty disables CSV export; otherwise save the displayed statistics to this path.

%% Read and display the requested variance-reduction table
vrf=view_variance_reduction(dataset, ...
    'Dimension',dimension, ...
    'NExponents',exponents, ...
    'TailCount',tailCount,'Bootstrap',bootstrap, ...
    'Tests',tests, ...
    'MatFile',matFile, ...
    'CsvFile',csvFile);
