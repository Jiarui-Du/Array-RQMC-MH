function reference = array_rqmc_stationary_reference(target, test_id)
%ARRAY_RQMC_STATIONARY_REFERENCE Return pi(f) when known reliably.

    switch char(test_id)
        case 'cdf'
            reference = double(target.cdf_probability);
        case 'x'
            if target.is_symmetric
                reference = 0;
            elseif strcmp(target.id, 'gaussian_mixture3')
                [reference, ~] = gaussian_mixture_moments();
            else
                reference = NaN;
            end
        case 'x2'
            switch target.id
                case 'normal'
                    reference = 1;
                case 'student_t3'
                    reference = 3;
                case 'gaussian_mixture3'
                    [~, reference] = gaussian_mixture_moments();
                case 'double_well_e1'
                    density = @(x) exp(-(x.^2 - 1).^2);
                    normalizer = integral(density, -Inf, Inf);
                    reference = integral(@(x) x.^2 .* density(x), ...
                        -Inf, Inf) / normalizer;
                otherwise
                    reference = NaN;
            end
        otherwise
            reference = NaN;
    end
end

function [first_moment, second_moment] = gaussian_mixture_moments()
    centers = [-1.5, 0, 1.5];
    precisions = [2, 3, 4];
    coefficients = [1, 0.8, 0.6];
    weights = coefficients ./ sqrt(precisions);
    weights = weights ./ sum(weights);
    first_moment = sum(weights .* centers);
    variances = 1 ./ (2 .* precisions);
    second_moment = sum(weights .* (centers.^2 + variances));
end
