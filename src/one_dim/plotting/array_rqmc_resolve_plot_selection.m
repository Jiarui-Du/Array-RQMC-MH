function ids = array_rqmc_resolve_plot_selection( ...
        requested, discovered, canonical, kind)
%ARRAY_RQMC_RESOLVE_PLOT_SELECTION Resolve all/explicit plot selections.
%
% For all, keeps discovered canonical IDs first, then appends other discovered
% IDs. Explicit missing requests are retained with a warning so the panel layer
% can display their missing-data diagnostics.


    %% Resolve requested IDs against discovered result metadata
    ids = array_rqmc_normalize_id_selection(requested);
    if numel(ids) == 1 && strcmpi(ids{1}, 'all')
        ids = canonical(ismember(canonical, discovered));
        ids = [ids, discovered(~ismember(discovered, ids)).'];
        return;
    end
    missing = ids(~ismember(ids, discovered));
    if ~isempty(missing)
        warning('array_rqmc:MissingPlotSelection', ...
            'Requested %s IDs not found: %s.', kind, strjoin(missing, ', '));
    end
end
