function [estimate, acceptance_rate] = ...
    run_array_rqmc_gaussian_rwm( ...
        N, T0, T, sigma, x0, logp0, log_pi, test_functions, ...
        sobol_base, stream)
%RUN_ARRAY_RQMC_GAUSSIAN_RWM Run one Array-RQMC Gaussian RWM replication.
%
% All supplied test functions are accumulated along the same Array-RQMC
% particle evolution.  This function performs no timing.  Standalone timing
% is handled externally by BENCHMARK_ARRAY_RQMC_SINGLE_TEST_TIMES.
%
% At every transition, a fresh MatousekAffineOwen randomization of a
% three-dimensional Sobol' point set is generated. Exact first-coordinate
% stratification assigns point ranks in O(N) operations via
%
%   rank(k) = floor(N * V_k) + 1.
%
% Inputs log_pi and each test handle must accept an N-vector. Returns one
% time-and-particle average per test over T0+1:T plus the acceptance fraction
% over all N*T proposals. Particle ordering is restored before the next step.


    %% Initialize test accumulators and constant particles
    test_functions = reshape(test_functions, 1, []);
    number_of_tests = numel(test_functions);

    X = repmat(x0, N, 1);  % The constant initial array is sorted.
    logp = repmat(logp0, N, 1);
    accumulated = zeros(1, number_of_tests);
    accepted_total = 0;


    %% Install the replication stream with automatic restoration
    % SCRAMBLE uses MATLAB's current global random stream.
    old_stream = RandStream.getGlobalStream;
    cleanup = onCleanup(@() RandStream.setGlobalStream(old_stream)); %#ok<NASGU>
    RandStream.setGlobalStream(stream);

    delta = 2^-53;


    %% Generate a fresh net and perform the Metropolis transition
    for step = 1:T
        randomized_set = scramble(sobol_base, 'MatousekAffineOwen');
        W = net(randomized_set, N);

        % Exact first-coordinate N-stratification implies that RANK is a
        % permutation of 1:N and is equivalent to sorting by W(:,1).
        rank = floor(N .* W(:, 1)) + 1;

        u_proposal = zeros(N, 1);
        u_accept = zeros(N, 1);
        u_proposal(rank) = W(:, 2);
        u_accept(rank) = W(:, 3);

        % Only the proposal coordinate is moved into the open unit interval,
        % because Phi^{-1}(0) and Phi^{-1}(1) are infinite.
        u_proposal = min(max(u_proposal, delta), 1 - delta);

        proposal = X + sigma .* standard_normal_icdf(u_proposal);
        logp_proposal = reshape(log_pi(proposal), [], 1);

        log_acceptance = min(logp_proposal - logp, 0);

        % The MH acceptance event is U <= alpha. Positive equality is
        % accepted. The alpha=0 boundary is rejected explicitly.
        accepted = isfinite(log_acceptance) & ...
            log(u_accept) <= log_acceptance;

        X(accepted) = proposal(accepted);
        logp(accepted) = logp_proposal(accepted);
        accepted_total = accepted_total + nnz(accepted);


        %% Accumulate retained observations and prepare the next ordering
        if step > T0
            for test_index = 1:number_of_tests
                values = reshape(test_functions{test_index}(X), [], 1);
                accumulated(test_index) = accumulated(test_index) + ...
                    mean(values);
            end
        end

        if step < T
            [X, order] = sort(X, 'ascend');
            logp = logp(order);
        end
    end


    %% Normalize estimates and acceptance counts
    estimate = accumulated ./ (T - T0);
    acceptance_rate = accepted_total / (N * T);
end


%% Local helper: standard_normal_icdf
% Map open-unit-interval uniforms to standard-normal innovations.
function z = standard_normal_icdf(u)
%STANDARD_NORMAL_ICDF Evaluate Phi^{-1}(u) without NORMINV.
    z = -sqrt(2) .* erfcinv(2 .* u);
end
