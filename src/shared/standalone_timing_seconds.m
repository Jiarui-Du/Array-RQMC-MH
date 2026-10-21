function seconds=standalone_timing_seconds(source,N,methods)
%STANDALONE_TIMING_SECONDS Read only independent serial median elapsed times.
% Missing timing, legacy concurrent timings, and unmeasured N remain NaN.
    seconds=nan(numel(N),numel(methods));
    if ~isfield(source,'timing') || isempty(source.timing),return;end
    t=source.timing;
    if ~isfield(t,'protocol') || ~strcmp(t.protocol,'serial-single-test-v1'),return;end
    [present,idx]=ismember(N,t.N);[found,mi]=ismember(methods,t.methods);
    assert(all(found),'Timing methods do not match statistical methods');
    for k=find(present)
        raw=t.seconds(:,idx(k),mi);
        assert(all(isfinite(raw) & raw>0,'all'),'Invalid standalone timing samples');
        seconds(k,:)=reshape(median(raw,1),1,[]);
    end
end
