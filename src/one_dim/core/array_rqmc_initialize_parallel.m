function use_parallel = array_rqmc_initialize_parallel(requested, auto_start)
%ARRAY_RQMC_INITIALIZE_PARALLEL Resolve the requested execution mode.
%
% Returns the actual execution mode after checking the request and pool.
% Unavailable tooling or a failed pool start falls back to serial with a warning.


    %% Check whether parallel execution is requested and available
    use_parallel = false;
    if ~requested
        return;
    end
    if isempty(ver('parallel'))
        warning('array_rqmc:NoParallelToolbox', ...
            'Parallel Computing Toolbox is unavailable; running serially.');
        return;
    end


    %% Reuse a pool or start one only when allowed
    pool = gcp('nocreate');
    if isempty(pool)
        if ~auto_start
            warning('array_rqmc:NoParallelPool', ...
                'No parallel pool exists; running serially.');
            return;
        end
        try
            parpool;
        catch ME
            warning('array_rqmc:PoolStartFailed', ...
                'Could not start a parallel pool: %s. Running serially.', ...
                ME.message);
            return;
        end
    end
    use_parallel = true;
end
