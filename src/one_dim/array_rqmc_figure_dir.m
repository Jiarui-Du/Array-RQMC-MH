function figure_dir = array_rqmc_figure_dir(category)
%ARRAY_RQMC_FIGURE_DIR Return and create the central figure directory.
%
%   ARRAY_RQMC_FIGURE_DIR() returns PROJECT_ROOT/figures.
%   ARRAY_RQMC_FIGURE_DIR(CATEGORY) returns a category subdirectory, such
%   as PROJECT_ROOT/figures/paper. Category names are sanitized so plot
%   artifacts cannot escape the central figure directory.
%
% Returns an absolute folder and creates it when absent.


    %% Resolve the central figure folder
    function_root = fileparts(mfilename('fullpath'));
    project_root = fileparts(fileparts(function_root));
    figure_dir = fullfile(project_root, 'figures');


    %% Sanitize an optional category
    if nargin >= 1 && ~isempty(category)
        if ~(ischar(category) || (isstring(category) && isscalar(category)))
            error('array_rqmc:FigureCategory', ...
                'Figure category must be a text scalar.');
        end
        category = regexprep(strtrim(char(category)), ...
            '[^A-Za-z0-9_-]+', '_');
        category = regexprep(category, '_+', '_');
        category = regexprep(category, '^_+|_+$', '');
        if isempty(category)
            error('array_rqmc:FigureCategory', ...
                'Figure category must contain a letter or number.');
        end
        figure_dir = fullfile(figure_dir, category);
    end


    %% Create the selected destination
    if ~isfolder(figure_dir)
        [ok, message] = mkdir(figure_dir);
        if ~ok
            error('array_rqmc:CreateFigureDirectory', ...
                'Could not create %s: %s', figure_dir, message);
        end
    end
end
