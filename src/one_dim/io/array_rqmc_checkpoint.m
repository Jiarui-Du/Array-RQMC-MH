function output = array_rqmc_checkpoint(action, checkpoint_file, payload)
%ARRAY_RQMC_CHECKPOINT Save, load, or delete one experiment checkpoint.
%
% save writes a temporary MAT beside the destination before replacing it.
% load returns empty for absent/unusable checkpoints; delete removes only the
% explicitly supplied checkpoint path. Callers decide compatibility.


    %% Dispatch the requested checkpoint operation
    checkpoint_file = char(checkpoint_file);
    output = [];
    switch lower(char(action))
        case 'save'
            parent = fileparts(checkpoint_file);
            if ~isempty(parent) && ~isfolder(parent)
                [ok, message] = mkdir(parent);
                if ~ok
                    error('array_rqmc:CreateCheckpointDirectory', ...
                        'Could not create %s: %s', parent, message);
                end
            end
            checkpoint = payload;
            temporary_parent = parent;
            if isempty(temporary_parent)
                temporary_parent = pwd;
            end
            temporary_file = [tempname(temporary_parent), '.mat'];
            temporary_cleanup = onCleanup( ...
                @() delete_if_present(temporary_file));
            save(temporary_file, 'checkpoint', '-v7.3');
            [ok, message] = movefile( ...
                temporary_file, checkpoint_file, 'f');
            if ~ok
                error('array_rqmc:ReplaceCheckpoint', ...
                    'Could not replace checkpoint %s: %s', ...
                    checkpoint_file, message);
            end
            clear temporary_cleanup;
            output = checkpoint_file;
        case 'load'
            if isfile(checkpoint_file)
                try
                    loaded = load(checkpoint_file, 'checkpoint');
                catch ME
                    warning('array_rqmc:CheckpointLoadFailed', ...
                        'Ignoring unreadable checkpoint %s: %s', ...
                        checkpoint_file, ME.message);
                    loaded = struct;
                end
                if isfield(loaded, 'checkpoint') && ...
                        isstruct(loaded.checkpoint)
                    output = loaded.checkpoint;
                end
            end
        case 'delete'
            if isfile(checkpoint_file)
                delete(checkpoint_file);
            end
        otherwise
            error('array_rqmc:CheckpointAction', ...
                'Unknown checkpoint action: %s', action);
    end
end

%% Local helper: delete_if_present
% Remove only the explicitly named temporary file when cleanup runs.
function delete_if_present(file)
    if isfile(file)
        delete(file);
    end
end
