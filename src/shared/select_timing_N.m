function N=select_timing_N(grid,requestedExponents)
%SELECT_TIMING_N Default to the largest statistical N; never extrapolate costs.
    if isempty(requestedExponents),N=max(grid);else,N=2.^requestedExponents;end
    assert(all(ismember(N,grid)) && all(diff(N)>0), ...
        'TimingNExponents must be an increasing subset of the statistical grid');
end
