function entries = array_rqmc_select_catalog_entries( ...
        all_entries, requested, kind)
%ARRAY_RQMC_SELECT_CATALOG_ENTRIES Select catalog entries by stable ID.
%
% Matches IDs case-insensitively and preserves the requested ordering.
% A sole all selection returns the catalog unchanged; unknown IDs raise an error.


    %% Resolve the all selector or explicit catalog IDs
    ids = array_rqmc_normalize_id_selection(requested);
    if numel(ids) == 1 && strcmpi(ids{1}, 'all')
        entries = all_entries;
        return;
    end

    all_ids = {all_entries.id};
    indices = zeros(1, numel(ids));
    for k = 1:numel(ids)
        match = find(strcmpi(ids{k}, all_ids), 1);
        if isempty(match)
            error('array_rqmc:UnknownID', ...
                'Unknown %s ID ''%s''. Available IDs: %s.', ...
                kind, ids{k}, strjoin(all_ids, ', '));
        end
        indices(k) = match;
    end
    entries = all_entries(indices);
end
