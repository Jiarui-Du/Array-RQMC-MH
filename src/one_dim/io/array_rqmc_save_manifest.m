function array_rqmc_save_manifest(output_dir, manifest, settings)
%ARRAY_RQMC_SAVE_MANIFEST Incrementally save the suite manifest.
%
% Writes manifest and manifest_settings to the resolved run manifest file.
% The suite calls this after each target so partial progress remains inspectable.


    %% Resolve the destination and save current suite progress
    manifest_settings = settings; %#ok<NASGU>
    if isfield(settings, 'manifest_file')
        manifest_file = settings.manifest_file;
    else
        manifest_file = fullfile(output_dir, 'array_rqmc_suite_manifest.mat');
    end
    save(manifest_file, 'manifest', 'manifest_settings', '-v7');
end
