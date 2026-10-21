function file_path = array_rqmc_unique_output_path(requested_path)
%ARRAY_RQMC_UNIQUE_OUTPUT_PATH Avoid overwriting an existing artifact.
%
% Returns the requested path if unused; otherwise appends a millisecond tag.
% This helper creates no file and does not reserve the returned path.


    %% Preserve an unused path or add a run timestamp
    file_path = char(requested_path);
    if ~isfile(file_path)
        return;
    end
    [folder, stem, extension] = fileparts(file_path);
    timestamp = char(datetime('now', 'Format', 'yyyyMMdd''T''HHmmssSSS'));
    file_path = fullfile(folder, sprintf('%s__run-%s%s', ...
        stem, timestamp, extension));
end
