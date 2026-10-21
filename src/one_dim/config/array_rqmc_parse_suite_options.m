function opt = array_rqmc_parse_suite_options(varargin)
%ARRAY_RQMC_PARSE_SUITE_OPTIONS Parse paper-suite configuration.
%
% Returns validated name/value settings used by the paper-suite runner.
% mVec contains base-2 particle exponents; TimingM selects a separate timing N.


    %% Declare target, simulation, timing and output options
    ip = inputParser;
    ip.FunctionName = 'run_array_rqmc_suite';
    addParameter(ip, 'TargetIDs', 'all', @id_selection);
    addParameter(ip, 'TestFunctionIDs', 'all', @id_selection);
    addParameter(ip, 'mVec', 10:18, @m_vector);
    addParameter(ip, 'T0', 8, @nonnegative_integer);
    addParameter(ip, 'T', 128, @positive_integer);
    addParameter(ip, 'R', 200, @replication_count);
    addParameter(ip, 'sigma', 2.38, @positive_scalar);
    addParameter(ip, 'InitialState', 0, @real_scalar);
    addParameter(ip, 'MasterSeed', 2026, @nonnegative_integer);
    addParameter(ip, 'UseParallel', true, @logical_scalar);
    addParameter(ip, 'AutoStartPool', true, @logical_scalar);
    addParameter(ip, 'Bootstrap', 1000, @nonnegative_integer);
    addParameter(ip, 'SlopeTailCount', 5, @(x) positive_integer(x) && x >= 2);
    addParameter(ip, 'RunTiming', true, @logical_scalar);
    addParameter(ip, 'TimingM', [], @(x) isempty(x) || m_scalar(x));
    addParameter(ip, 'TimingR', 20, @replication_count);
    addParameter(ip, 'TimingWarmupM', 12, @m_scalar);
    addParameter(ip, 'OutputDir', 'array_rqmc_results', @text_scalar);
    addParameter(ip, 'ContinueOnError', true, @logical_scalar);
    parse(ip, varargin{:});
    opt = ip.Results;


    %% Require a nonempty retained observation window
    if opt.T0 >= opt.T
        error('array_rqmc:TimeIndices', 'Require 0 <= T0 < T.');
    end
end

%% Local input validators
% Predicates used by the option parser; they do not alter supplied values.
function tf = id_selection(x), tf = ischar(x) || isstring(x) || iscell(x); end
function tf = m_vector(x), tf = isnumeric(x) && isvector(x) && ...
    ~isempty(x) && all(isfinite(x)) && all(x == floor(x)) && ...
    all(x >= 1) && all(x <= 30); end
function tf = nonnegative_integer(x), tf = isnumeric(x) && isscalar(x) && ...
    isfinite(x) && x == floor(x) && x >= 0; end
function tf = positive_integer(x), tf = nonnegative_integer(x) && x >= 1; end
function tf = m_scalar(x), tf = positive_integer(x) && x <= 30; end
function tf = replication_count(x), tf = positive_integer(x) && x >= 2; end
function tf = positive_scalar(x), tf = isnumeric(x) && isreal(x) && ...
    isscalar(x) && isfinite(x) && x > 0; end
function tf = real_scalar(x), tf = isnumeric(x) && isreal(x) && ...
    isscalar(x) && isfinite(x); end
function tf = logical_scalar(x), tf = islogical(x) && isscalar(x); end
function tf = text_scalar(x), tf = ischar(x) || ...
    (isstring(x) && isscalar(x)); end
