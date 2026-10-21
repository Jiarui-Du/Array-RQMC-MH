function opt = array_rqmc_parse_timing_options(varargin)
%ARRAY_RQMC_PARSE_TIMING_OPTIONS Parse standalone timing configuration.
%
% Returns settings for complete serial runs of each individual test.
% WarmupM is a particle exponent; warmup runs are excluded from recorded timings.


    %% Declare benchmark inputs and metadata
    ip = inputParser;
    ip.FunctionName = 'benchmark_array_rqmc_single_test_times';
    addParameter(ip, 'm', 20, @m_scalar);
    addParameter(ip, 'T0', 8, @nonnegative_integer);
    addParameter(ip, 'T', 128, @positive_integer);
    addParameter(ip, 'R', 20, @replication_count);
    addParameter(ip, 'sigma', 2.38, @positive_scalar);
    addParameter(ip, 'LogPi', @(x) -0.5 .* x.^2, ...
        @(x) isa(x, 'function_handle'));
    addParameter(ip, 'TestFunctions', {@(x) x}, @test_functions);
    addParameter(ip, 'TestFunctionIDs', {}, @text_list);
    addParameter(ip, 'TestFunctionLabels', {}, @text_list);
    addParameter(ip, 'InitialState', 0, @real_scalar);
    addParameter(ip, 'MasterSeed', 2026, @nonnegative_integer);
    addParameter(ip, 'WarmupM', [], @(x) isempty(x) || m_scalar(x));
    parse(ip, varargin{:});
    opt = ip.Results;

    %% Require a nonempty retained observation window
    if opt.T0 >= opt.T
        error('array_rqmc:TimeIndices', 'Require 0 <= T0 < T.');
    end
end

%% Local input validators
% Predicates used by the option parser; they do not alter supplied values.
function tf = nonnegative_integer(x), tf = isnumeric(x) && isscalar(x) && ...
    isfinite(x) && x == floor(x) && x >= 0; end
function tf = positive_integer(x), tf = nonnegative_integer(x) && x >= 1; end
function tf = m_scalar(x), tf = positive_integer(x) && x <= 30; end
function tf = replication_count(x), tf = positive_integer(x) && x >= 2; end
function tf = positive_scalar(x), tf = isnumeric(x) && isreal(x) && ...
    isscalar(x) && isfinite(x) && x > 0; end
function tf = real_scalar(x), tf = isnumeric(x) && isreal(x) && ...
    isscalar(x) && isfinite(x); end
function tf = test_functions(x), tf = isa(x, 'function_handle') || ...
    (iscell(x) && ~isempty(x) && ...
    all(cellfun(@(h) isa(h, 'function_handle'), x(:)))); end
function tf = text_list(x), tf = isempty(x) || ischar(x) || ...
    isstring(x) || (iscell(x) && all(cellfun(@(v) ischar(v) || ...
    (isstring(v) && isscalar(v)), x(:)))); end
