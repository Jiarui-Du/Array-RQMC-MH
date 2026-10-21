function timing = array_rqmc_extract_timing_at_N( ...
        results, target_N, target_position)
%ARRAY_RQMC_EXTRACT_TIMING_AT_N Read current or legacy timing fields.
%
% Prefers separate benchmark timings when their N matches the requested value.
% Otherwise reads legacy in-simulation fields by target_position. Missing
% values remain NaN so exports can report unavailable efficiency.


    %% Read a matching standalone benchmark when present
    timing = struct('stored_tae', NaN, 'median_iid', NaN, ...
        'median_array', NaN, 'R', NaN, 'scope', '', ...
        'source', 'timing unavailable');
    if isfield(results, 'timing') && isstruct(results.timing)
        value = results.timing;
        if scalar_number(value, 'N') && ...
                abs(double(value.N) - target_N) <= 8 * eps(target_N)
            timing.median_iid = number_or_nan(value, 'median_time_mcmc');
            timing.median_array = number_or_nan(value, 'median_time_sobol');
            timing.stored_tae = number_or_nan( ...
                value, 'time_adjusted_efficiency_sobol');
            timing.R = number_or_nan(value, 'R');
            timing.scope = text_or_empty(value, 'scope');
            timing.source = 'separate standalone benchmark';
            if ~isfinite(timing.median_iid) && isfield(value, 'time_mcmc')
                timing.median_iid = median(double(value.time_mcmc(:)));
            end
            if ~isfinite(timing.median_array) && isfield(value, 'time_sobol')
                timing.median_array = median(double(value.time_sobol(:)));
            end
            return;
        end
    end


    %% Fall back to historical timing fields
    timing.source = 'legacy in-simulation timing';
    if isfield(results, 'settings') && isstruct(results.settings)
        timing.scope = text_or_empty(results.settings, 'timing_scope');
    end
    timing.stored_tae = vector_value( ...
        results, 'time_adjusted_efficiency_sobol', target_position);
    timing.median_iid = legacy_median( ...
        results, 'median_time_mcmc', 'time_mcmc', target_position);
    timing.median_array = legacy_median( ...
        results, 'median_time_sobol', 'time_sobol', target_position);
end

%% Local helper: legacy_median
% Prefer a stored median, otherwise calculate it from the selected raw timing column.
function value = legacy_median(results, median_name, raw_name, position)
    value = vector_value(results, median_name, position);
    if isfinite(value) || ~isfield(results, raw_name)
        return;
    end
    raw = double(results.(raw_name));
    if ismatrix(raw) && size(raw, 2) >= position
        value = median(raw(:, position));
    end
end

%% Local helper: vector_value
% Read the requested legacy vector element, leaving absent values as NaN.
function value = vector_value(results, name, position)
    value = NaN;
    if isfield(results, name)
        values = double(results.(name)(:));
        if numel(values) >= position
            value = values(position);
        end
    end
end

%% Local helper: scalar_number
% Recognize an existing finite numeric scalar field.
function tf = scalar_number(s, name)
    tf = isfield(s, name) && isnumeric(s.(name)) && ...
        isscalar(s.(name)) && isfinite(s.(name));
end

%% Local helper: number_or_nan
% Extract a finite scalar field or preserve the missing-value marker.
function value = number_or_nan(s, name)
    value = NaN;
    if scalar_number(s, name)
        value = double(s.(name));
    end
end

%% Local helper: text_or_empty
% Extract scalar text metadata or return an empty label.
function value = text_or_empty(s, name)
    value = '';
    if isfield(s, name) && (ischar(s.(name)) || ...
            (isstring(s.(name)) && isscalar(s.(name))))
        value = char(s.(name));
    end
end
