function stream = array_rqmc_make_method_stream( ...
        master_seed, method_offset, substream)
%ARRAY_RQMC_MAKE_METHOD_STREAM Construct the established method stream.
%
% The seed formula and mrg32k3a substream assignment are part of the
% experiment's reproducibility contract and must not be changed.
%
% Returns a local mrg32k3a stream for one method and substream index.
% The caller passes this stream explicitly or installs it temporarily for scrambling.


    %% Apply the fixed seed offset and select the requested substream
    modulus = 2^32 - 1;
    seed = mod( ...
        double(master_seed) + 104729 * double(method_offset), modulus);
    stream = RandStream('mrg32k3a', 'Seed', seed);
    stream.Substream = double(substream);
end
