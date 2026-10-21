function file_path = array_rqmc_auto_filename( ...
        output_dir, target_id, test_id, m_vec, T0, T, R, sigma, ...
        initial_state, seed, joint_test_ids)
%ARRAY_RQMC_AUTO_FILENAME Construct a readable result filename.
%
% JOINT_TEST_IDS is optional. When several test functions are accumulated
% in the same target run, their IDs are included in the filename so that a
% later run with a different test collection does not silently overwrite it.
%
% Example:
% array_rqmc__target-normal__f-cdf__joint-cdf-x-x2...
%     __Nexp-10to18__T0-8__T-128__R-200...
%
% Returns a full MAT path and creates the output directory when needed.
% Name encoding describes settings; collision avoidance is handled separately
% by array_rqmc_unique_output_path.


    %% Resolve optional joint-test IDs and output directory
    if nargin < 11
        joint_test_ids = {};
    end

    output_dir = char(output_dir);
    if isempty(output_dir)
        output_dir = '.';
    end
    if ~exist(output_dir, 'dir')
        [ok, message] = mkdir(output_dir);
        if ~ok
            error('array_rqmc_auto_filename:CreateDirectory', ...
                'Could not create output directory ''%s'': %s', ...
                output_dir, message);
        end
    end


    %% Encode settings into filesystem-safe filename components
    target_tag = sanitize_component(target_id);
    test_tag = sanitize_component(test_id);
    m_tag = encode_m_vector(m_vec);
    sigma_tag = encode_real_number(sigma);
    x0_tag = encode_real_number(initial_state);
    joint_tag = encode_id_list(joint_test_ids);

    if isempty(joint_tag)
        joint_component = '';
    else
        joint_component = ['__joint-', joint_tag];
    end


    %% Assemble the descriptive result path
    filename = sprintf([ ...
        'array_rqmc__target-%s__f-%s%s__Nexp-%s', ...
        '__T0-%d__T-%d__R-%d__sigma-%s__x0-%s', ...
        '__seed-%d.mat'], ...
        target_tag, test_tag, joint_component, m_tag, ...
        round(T0), round(T), round(R), sigma_tag, x0_tag, ...
        round(seed));

    file_path = fullfile(output_dir, filename);
end


%% Local helper: encode_id_list
% Join sanitized test IDs while preserving their requested order.
function tag = encode_id_list(values)
    if isempty(values)
        tag = '';
        return;
    end

    if ischar(values)
        values = {values};
    elseif isstring(values)
        values = cellstr(values(:)).';
    else
        values = values(:).';
    end

    pieces = cellfun(@sanitize_component, values, 'UniformOutput', false);
    tag = strjoin(pieces, '-');
end


%% Local helper: encode_m_vector
% Encode a contiguous exponent range compactly, or list a noncontiguous grid.
function tag = encode_m_vector(m_vec)
    m = sort(unique(reshape(double(m_vec), 1, [])));
    if isempty(m)
        error('array_rqmc_auto_filename:EmptyMVec', ...
            'm_vec must be nonempty.');
    end

    if numel(m) > 1 && all(diff(m) == 1)
        tag = sprintf('%dto%d', m(1), m(end));
    elseif numel(m) == 1
        tag = sprintf('%d', m);
    else
        pieces = arrayfun(@(z) sprintf('%d', z), m, ...
            'UniformOutput', false);
        tag = strjoin(pieces, '_');
    end
end


%% Local helper: encode_real_number
% Replace signs and decimal points with filename-safe number tokens.
function tag = encode_real_number(x)
    tag = sprintf('%.12g', double(x));
    tag = strrep(tag, '-', 'm');
    tag = strrep(tag, '+', '');
    tag = strrep(tag, '.', 'p');
end


%% Local helper: sanitize_component
% Normalize an ID to a nonempty filename-safe token.
function tag = sanitize_component(value)
    tag = lower(strtrim(char(value)));
    tag = regexprep(tag, '[^a-z0-9_-]+', '-');
    tag = regexprep(tag, '-+', '-');
    tag = regexprep(tag, '(^-|-$)', '');

    if isempty(tag)
        error('array_rqmc_auto_filename:EmptyID', ...
            'IDs must contain letters or digits.');
    end
end
