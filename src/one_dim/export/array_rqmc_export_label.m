function label = array_rqmc_export_label(kind, id)
%ARRAY_RQMC_EXPORT_LABEL Return plain-text or LaTeX case labels.
%
% Maps stable target/test IDs to display labels for CSV or TeX.
% Unknown IDs are passed through; TeX labels reuse the plotting label helper.


    %% Resolve the requested label family
    kind = lower(char(kind));
    id = char(id);
    switch kind
        case 'target_plain'
            switch id
                case 'normal', label = 'Standard normal';
                case 'double_well_e1', label = 'Double well';
                case 'gaussian_mixture3', label = 'Gaussian mixture';
                case 'student_t3', label = 'Student-t3';
                otherwise, label = id;
            end
        case 'target_tex'
            label = array_rqmc_plot_label('target', id);
        case 'test_plain'
            switch id
                case 'cdf', label = 'f_CDF(x)=1{x<=a_pi}';
                case 'x', label = 'f(x)=x';
                case 'x2', label = 'f(x)=x^2';
                otherwise, label = id;
            end
        case 'test_tex'
            label = array_rqmc_plot_label('test', id);
        otherwise
            label = id;
    end
end
