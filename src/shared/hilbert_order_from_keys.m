function [order,collisions] = hilbert_order_from_keys(keys,bits)
%HILBERT_ORDER_FROM_KEYS Stable order of keys produced by hilbert_keys.
% Keys must be in the current particle order, including after MH updates.
%
% Returns the row permutation; optional collisions counts adjacent equal cells
% after sorting. In one dimension the permutation uses exact coordinates while
% the collision diagnostic still uses the requested finite-bit cells.

    %% Dispatch sorting and collision counts by key representation
    if isa(keys,'double')
        [sorted,order]=sort(keys(:,1));
        if nargout>1
            cells=floor(min(sorted,1-eps)*2^bits);
            collisions=nnz(cells(2:end)==cells(1:end-1));
        end
    elseif isa(keys,'uint64')
        [sorted,order]=sort(keys);
        if nargout>1
            collisions=nnz(sorted(2:end)==sorted(1:end-1));
        end
    else
        [sorted,order]=sortrows(keys);
        if nargout>1
            collisions=nnz(all(sorted(2:end,:)==sorted(1:end-1,:),2));
        end
    end
end
