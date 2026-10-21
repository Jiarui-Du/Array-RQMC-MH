function [order,collisions] = hilbert_sort_order(U,bits)
%HILBERT_SORT_ORDER Skilling Hilbert order; equal cells retain input order.
% Encoding is separate from ordering to allow MH to reuse unchanged keys.
% d*bits<=64 uses uint64 + sort; larger keys use uint32 words + sortrows.
% d=1 keeps the historical exact ascending order and finite-bit collisions.
%
% U contains one point per row; order indexes those rows in Hilbert order.
% Requesting only order avoids computing the optional collision diagnostic.

    %% Encode coordinates and sort the resulting keys
    if nargin<2 || isempty(bits),bits=resolve_hilbert_bits(size(U,2));end
    keys=hilbert_keys(U,bits);
    if nargout>1
        [order,collisions]=hilbert_order_from_keys(keys,bits);
    else
        order=hilbert_order_from_keys(keys,bits);
    end
end
