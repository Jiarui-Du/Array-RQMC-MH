function fit = array_rqmc_fit_loglog_slope(x_values, y_values, tail_count)
%ARRAY_RQMC_FIT_LOGLOG_SLOPE Fit log10(y) on log10(x).
%
% Fixes the last tail_count pairs before checking positive finite variances.
% The caller controls ordering (normally increasing N). Fewer than two points
% leave the fit undefined; standard error additionally requires more than two.


    %% Initialize missing-value defaults and filter usable pairs
    fit = struct( ...
        'slope', NaN, 'intercept', NaN, 'standard_error', NaN, ...
        'r_squared', NaN, 'n_points', 0, 'x_used', []);
    x_values = double(x_values(:));
    y_values = double(y_values(:));
    if numel(x_values) < 2
        return;
    end


    %% Fit the selected tail by ordinary least squares
    count = min(numel(x_values), max(2, round(tail_count)));
    x_values = x_values(end-count+1:end);
    y_values = y_values(end-count+1:end);
    fit.n_points=count;fit.x_used=x_values;
    if any(~isfinite(x_values) | x_values<=0 | ~isfinite(y_values) | y_values<=0)
        return;
    end
    x = log10(x_values);
    y = log10(y_values);
    design = [ones(count, 1), x];
    coefficient = design \ y;
    residual = y - design * coefficient;


    %% Record the fit and residual diagnostics
    fit.intercept = coefficient(1);
    fit.slope = coefficient(2);
    fit.n_points = count;
    fit.x_used = x_values;
    total_ss = sum((y - mean(y)).^2);
    residual_ss = sum(residual.^2);
    if total_ss > 0
        fit.r_squared = 1 - residual_ss / total_ss;
    end
    if count > 2
        centered_x_ss = sum((x - mean(x)).^2);
        if centered_x_ss > 0
            fit.standard_error = sqrt( ...
                (residual_ss / (count - 2)) / centered_x_ss);
        end
    end
end
