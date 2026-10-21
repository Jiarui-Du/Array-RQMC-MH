function label = array_rqmc_plot_label(kind, id)
%ARRAY_RQMC_PLOT_LABEL Return stable labels for figures and exports.
%
% Maps stable IDs to labels accepted by the LaTeX figure interpreter.
% Unknown target/test IDs have underscores escaped for readable fallback labels.


    %% Resolve target or test display text
    kind = lower(char(kind));
    id = char(id);
    switch kind
        case 'target'
            switch id
                case 'normal', label = 'Standard normal';
                case 'double_well_e1', label = 'Double well';
                case 'gaussian_mixture3', label = 'Gaussian mixture';
                case 'student_t3', label = 'Student-$t_3$';
                otherwise, label = strrep(id, '_', '\_');
            end
        case 'test'
            switch id
                case 'cdf', label = '$f_{\mathrm{CDF}}(x)$';
                case 'x', label = '$f(x)=x$';
                case 'x2', label = '$f(x)=x^2$';
                otherwise, label = strrep(id, '_', '\_');
            end
        otherwise
            label = id;
    end
end
