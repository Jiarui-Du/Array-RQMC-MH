function bits=resolve_hilbert_bits(p,requested)
%RESOLVE_HILBERT_BITS Default to an exact scalar uint64 key in multiple dimensions.
% Empty requested selects min(20,floor(64/p)). Explicit depths are retained
% for reference comparisons; hilbert_keys then supports multiword keys.
    validateattributes(p,{'double'},{'scalar','integer','positive','<=',64});
    if nargin<2 || isempty(requested)
        bits=min(20,floor(64/p));
    else
        validateattributes(requested,{'double'},{'scalar','integer','>=',1,'<=',30});
        bits=requested;
    end
end
