function [selected, status] = array_rqmc_select_result_file( ...
        file_index, target_id, test_id, required_R, target_N)
%ARRAY_RQMC_SELECT_RESULT_FILE Select the best file for one target/test pair.
%
% Returns one candidate row and a status string, or an empty table on failure.
% Preference order is requested N, largest completed N, valid variance pairs,
% replication count, then modification time; RequiredR is a hard filter.


    %% Resolve optional filters and select matching target/test IDs
    if nargin < 4
        required_R = [];
    end
    if nargin < 5
        target_N = [];
    end

    selected = file_index([],:);
    matches = strcmp(file_index.TargetID, char(target_id)) & ...
        strcmp(file_index.TestFunctionID, char(test_id));
    candidates = file_index(matches, :);

    if isempty(candidates)
        status = 'no matching result file';
        return;
    end


    %% Enforce the requested replication count
    if ~isempty(required_R)
        keep = isfinite(candidates.Replications) & ...
            candidates.Replications == double(required_R);
        candidates = candidates(keep, :);
        if isempty(candidates)
            status = sprintf('no file with R=%d', round(required_R));
            return;
        end
    end


    %% Rank candidates by coverage, completeness and recency
    has_target_N = zeros(height(candidates), 1);
    if ~isempty(target_N)
        for row = 1:height(candidates)
            N = candidates.NValues{row};
            has_target_N(row) = any(abs(N - double(target_N)) <= ...
                8 * eps(double(target_N)));
        end
    end

    score = table( ...
        -has_target_N, ...
        -replace_nan(candidates.NMaxCompleted), ...
        -candidates.ValidVariancePairs, ...
        -replace_nan(candidates.Replications), ...
        -candidates.ModifiedDatenum, ...
        (1:height(candidates)).', ...
        'VariableNames', {'TargetN', 'NMax', 'ValidPairs', ...
                          'R', 'Modified', 'OriginalRow'});
    [~, order] = sortrows(score, ...
        {'TargetN', 'NMax', 'ValidPairs', 'R', 'Modified'});
    selected = candidates(order(1), :);
    status = 'selected';
end


%% Local helper: replace_nan
% Rank unavailable coverage/count values below finite candidate scores.
function values = replace_nan(values)
    values = double(values);
    values(~isfinite(values)) = -Inf;
end
