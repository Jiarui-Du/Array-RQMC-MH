function project_root = array_rqmc_setup()
%ARRAY_RQMC_SETUP Add the project's source modules to the MATLAB path.
%
% Call this once from a top-level script. Compatibility entry points also
% call it so existing workflows that start in the project root keep working.
%
% Returns the repository root and adds the one-dimensional modules and shared helpers.
% Paths are resolved from this file so compatibility wrappers can be called
% without relying on the current working directory.


    %% Resolve the package root and register source modules
    function_root = fileparts(mfilename('fullpath'));
    project_root = fileparts(fileparts(function_root));
    addpath(function_root);
    addpath(fullfile(project_root,'src','shared'));
    module_names = { ...
        'core', 'config', 'experiments', 'analysis', ...
        'io', 'plotting', 'export'};

    for module_index = 1:numel(module_names)
        module_dir = fullfile(function_root, module_names{module_index});
        if isfolder(module_dir)
            addpath(module_dir);
        end
    end
end
