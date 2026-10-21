function [file_path, result_tag] = array_rqmc_chain_scaling_filename( ...
        output_dir, target_id, test_id, N_vec, T0, T_vec, R, sigma, ...
        initial_state, seed)
%ARRAY_RQMC_CHAIN_SCALING_FILENAME Name one multi-N chain-scaling result.

    m_vec = round(log2(double(N_vec(:).')));
    if any(double(N_vec(:).') ~= 2.^m_vec)
        error('array_rqmc:ParticlePowersOfTwo', ...
            'Every particle count must be an exact power of two.');
    end

    target_tag = sanitize_component(target_id);
    test_tag = sanitize_component(test_id);
    N_tag = encode_integer_vector(m_vec);
    T_tag = encode_time_vector(double(T_vec(:).'));
    result_tag = sprintf([ ...
        'array_rqmc_chain_scaling__target-%s__f-%s', ...
        '__Nexp-%s__T0-%d__%s__R-%d__sigma-%s__x0-%s__seed-%d'], ...
        target_tag, test_tag, N_tag, round(T0), T_tag, round(R), ...
        encode_real_number(sigma), encode_real_number(initial_state), ...
        round(seed));
    file_path = fullfile(char(output_dir), [result_tag, '.mat']);
end

function tag = encode_time_vector(values)
    exponents = log2(values);
    if all(values > 0) && all(exponents == round(exponents))
        tag = ['Texp-', encode_integer_vector(round(exponents))];
    else
        tag = ['T-', encode_integer_vector(values)];
    end
end

function tag = encode_integer_vector(values)
    values = double(values(:).');
    if numel(values) > 1 && all(diff(values) == 1)
        tag = sprintf('%dto%d', values(1), values(end));
    elseif isscalar(values)
        tag = sprintf('%d', values);
    else
        pieces = arrayfun(@(x) sprintf('%d', x), values, ...
            'UniformOutput', false);
        tag = strjoin(pieces, '_');
    end
end

function tag = encode_real_number(value)
    tag = sprintf('%.12g', double(value));
    tag = strrep(tag, '-', 'm');
    tag = strrep(tag, '+', '');
    tag = strrep(tag, '.', 'p');
end

function tag = sanitize_component(value)
    tag = lower(strtrim(char(value)));
    tag = regexprep(tag, '[^a-z0-9_-]+', '-');
    tag = regexprep(tag, '-+', '-');
    tag = regexprep(tag, '(^-|-$)', '');
    if isempty(tag)
        error('array_rqmc:EmptyID', ...
            'Target and test IDs must contain letters or digits.');
    end
end
