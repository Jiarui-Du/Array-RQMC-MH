function ids = array_rqmc_normalize_id_selection(value)
%ARRAY_RQMC_NORMALIZE_ID_SELECTION Normalize text IDs to a row cell array.
%
% Accepts char, string or cell text and preserves selection order. The result
% is a row cell array of char values; this helper does not validate catalog IDs.


    %% Normalize supported text containers
    if ischar(value)
        ids = {value};
    elseif isstring(value)
        ids = cellstr(value(:)).';
    elseif iscell(value)
        ids = value(:).';
        valid = cellfun(@(x) ischar(x) || ...
            (isstring(x) && isscalar(x)), ids);
        if ~all(valid)
            error('array_rqmc:IDSelection', ...
                'Every selected ID must be text.');
        end
        ids = cellfun(@char, ids, 'UniformOutput', false);
    else
        error('array_rqmc:IDSelection', ...
            'ID selections must be text or a cell array of text values.');
    end
end
