function opt = array_rqmc_parse_chain_options(varargin)
%ARRAY_RQMC_PARSE_CHAIN_OPTIONS Parse the chain-length experiment settings.

    ip = inputParser;
    ip.FunctionName = 'run_array_rqmc_chain_length_experiment';
    addParameter(ip, 'TargetID', 'gaussian_mixture3', @text_scalar);
    addParameter(ip, 'TestFunctionID', 'cdf', @text_scalar);
    addParameter(ip, 'N', 2^14, @positive_integer);
    addParameter(ip, 'T0', 2^7, @nonnegative_integer);
    addParameter(ip, 'TVec', 2.^(8:16), @integer_vector);
    addParameter(ip, 'R', 50, @replication_count);
    addParameter(ip, 'sigma', 2.38, @positive_scalar);
    addParameter(ip, 'InitialState', 0, @real_scalar);
    addParameter(ip, 'MasterSeed', 2026, @nonnegative_integer);
    addParameter(ip, 'UseParallel', true, @logical_scalar);
    addParameter(ip, 'AutoStartPool', true, @logical_scalar);
    addParameter(ip, 'OutputDir', fullfile(pwd, 'results', 'one_dim_horizon', 'individual'), @text_scalar);
    addParameter(ip, 'ResumeFromCheckpoint', true, @logical_scalar);
    addParameter(ip, 'ConfidenceLevel', 0.95, ...
        @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0 && x < 1);
    parse(ip, varargin{:});
    opt = ip.Results;
    opt.TargetID = char(opt.TargetID);
    opt.TestFunctionID = char(opt.TestFunctionID);
    opt.OutputDir = char(opt.OutputDir);
    opt.TVec = reshape(double(opt.TVec), 1, []);
    if any(opt.TVec <= opt.T0)
        error('array_rqmc:TimeIndices', ...
            'Every chain length must satisfy T > T0.');
    end
end

function tf = text_scalar(x), tf = ischar(x) || ...
    (isstring(x) && isscalar(x)); end
function tf = nonnegative_integer(x), tf = isnumeric(x) && isscalar(x) && ...
    isfinite(x) && x == floor(x) && x >= 0; end
function tf = positive_integer(x), tf = nonnegative_integer(x) && x >= 1; end
function tf = integer_vector(x), tf = isnumeric(x) && isvector(x) && ...
    ~isempty(x) && all(isfinite(x)) && all(x == floor(x)) && all(x >= 1); end
function tf = replication_count(x), tf = positive_integer(x) && x >= 2; end
function tf = positive_scalar(x), tf = isnumeric(x) && isreal(x) && ...
    isscalar(x) && isfinite(x) && x > 0; end
function tf = real_scalar(x), tf = isnumeric(x) && isreal(x) && ...
    isscalar(x) && isfinite(x); end
function tf = logical_scalar(x), tf = islogical(x) && isscalar(x); end
