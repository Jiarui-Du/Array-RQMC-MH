function seed=standalone_timing_seed(master,caseID,N,methodID,r,isWarmup)
%STANDALONE_TIMING_SEED Separate timing namespace, recorded by each runner.
% At most 499 samples keeps warmup/sample blocks disjoint. Runners also check
% the actual timing seed set against their saved statistical seeds.
    assert(r>=1 && r<=499);
    seed=mod(master+2^31+caseID*1000000+log2(N)*10000+methodID*1000+r+500*isWarmup,2^32);
end
