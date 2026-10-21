function print_experiment_header(title,tests,exponents,T0,T,R,parallel)
% Consistent statistical-stage banner; no random or simulation state changes.
    fprintf('============================================================\n');
    fprintf('%s\n',title);
    fprintf('Tests: %s\n',strjoin(tests,', '));
    fprintf('m values: %s; T0=%d; T=%d; R=%d\n',mat2str(exponents),T0,T,R);
    fprintf('Parallel execution: %d\n',parallel);
    fprintf('============================================================\n');
end
