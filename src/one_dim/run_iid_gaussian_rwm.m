function [estimate, acceptance_rate] = ...
    run_iid_gaussian_rwm( ...
        N, T0, T, sigma, x0, logp0, log_pi, test_functions, stream)
%RUN_IID_GAUSSIAN_RWM Run N independent Gaussian RWM chains.
%
% All supplied test functions are accumulated along the same particle
% trajectories.  This function performs no timing.  Standalone timing is
% handled externally by BENCHMARK_ARRAY_RQMC_SINGLE_TEST_TIMES, which calls
% this function once per test function.
%
% Inputs log_pi and test handles operate on an N-vector. Returns averages
% over particles and transitions T0+1:T; acceptance includes all N*T proposals.
% All random draws use the supplied stream rather than the global stream.


    %% Initialize particles, log densities and test accumulators
    test_functions = reshape(test_functions, 1, []);
    number_of_tests = numel(test_functions);

    X = repmat(x0, N, 1);
    logp = repmat(logp0, N, 1);
    accumulated = zeros(1, number_of_tests);
    accepted_total = 0;


    %% Apply independent Gaussian proposals and MH decisions
    for step = 1:T
        proposal = X + sigma .* randn(stream, N, 1);
        logp_proposal = reshape(log_pi(proposal), [], 1);

        u_accept = rand(stream, N, 1);
        log_acceptance = min(logp_proposal - logp, 0);

        % The MH acceptance event is U <= alpha. Positive equality is
        % accepted. The alpha=0 boundary is rejected explicitly.
        accepted = isfinite(log_acceptance) & ...
            log(u_accept) <= log_acceptance;

        X(accepted) = proposal(accepted);
        logp(accepted) = logp_proposal(accepted);
        accepted_total = accepted_total + nnz(accepted);


        %% Accumulate the retained observation window
        if step > T0
            for test_index = 1:number_of_tests
                values = reshape(test_functions{test_index}(X), [], 1);
                accumulated(test_index) = accumulated(test_index) + ...
                    mean(values);
            end
        end
    end


    %% Normalize estimates and acceptance counts
    estimate = accumulated ./ (T - T0);
    acceptance_rate = accepted_total / (N * T);
end
