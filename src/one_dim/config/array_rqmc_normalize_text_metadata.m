function metadata = array_rqmc_normalize_text_metadata(value, count, prefix)
%ARRAY_RQMC_NORMALIZE_TEXT_METADATA Normalize IDs or labels.
%
% Returns one text item per test, generating prefix-plus-index names only
% when metadata are omitted. A length mismatch is rejected.


    %% Normalize metadata and enforce the test count
    if isempty(value)
        metadata = arrayfun(@(k) sprintf('%s%d', prefix, k), 1:count, ...
            'UniformOutput', false);
    elseif ischar(value)
        metadata = {value};
    elseif isstring(value)
        metadata = cellstr(value(:)).';
    else
        metadata = value(:).';
        metadata = cellfun(@char, metadata, 'UniformOutput', false);
    end
    if numel(metadata) ~= count
        error('array_rqmc:TestMetadata', ...
            'Test metadata must match the test-function count.');
    end
end
