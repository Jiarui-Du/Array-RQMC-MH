function report=validate_1d_horizon
% Small migration check only: N<=32, T<=16, R=4. Never starts a formal run.
    root=fileparts(fileparts(mfilename('fullpath')));addpath(genpath(fullfile(root,'src')));
    folder=tempname;mkdir(folder);cleanup=onCleanup(@()rmdir(folder,'s'));
    args={'mVec',[4 5],'T0',2,'TVec',[4 8 16],'R',4,'sigma',2.38, ...
        'SlopeTailCount',3,'Bootstrap',20,'OutputDir',folder,'AutoStartPool',false};
    [a,art]=run_array_rqmc_horizon_suite(args{:},'UseParallel',false);
    assert(a.sigma==2.38 && all(isfinite(a.variance_array),'all'));
    assert(all(isfile({art.mat,art.summary_csv,art.slopes_csv})));
    [loaded,~]=load_array_rqmc_horizon_result(art.mat);
    assert(isequal(loaded.variance_array,a.variance_array));
    pool=gcp('nocreate');owned=isempty(pool);
    if owned,pool=parpool('Processes',2);end
    poolCleanup=onCleanup(@()close_owned_pool(pool,owned));
    b=run_array_rqmc_horizon_suite(args{:},'UseParallel',true);
    for k=1:2
        assert(isequal(a.individual_results{k}.est_mcmc_all,b.individual_results{k}.est_mcmc_all));
        assert(isequal(a.individual_results{k}.est_sobol_all,b.individual_results{k}.est_sobol_all));
        assert(isequaln(a.bootstrap{k},b.bootstrap{k}));
    end
    c=run_array_rqmc_horizon_suite(args{:},'UseParallel',false,'ResumeDir',fileparts(art.mat));
    assert(isequaln(a.slope_table,c.slope_table));
    rejected=false;
    try
        run_array_rqmc_horizon_suite(args{:},'UseParallel',false,'sigma',2.4,'ResumeDir',fileparts(art.mat));
    catch e
        rejected=strcmp(e.identifier,'array_rqmc:IncompatibleChainResult');
    end
    assert(rejected,'Changed sigma must not reuse old results.');
    [fig,stem]=plot_array_rqmc_chain_length_scaling_results(a,fullfile(folder,'figure'),'Visible','off');
    assert(all(isfile({[stem '.pdf'],[stem '.fig']})));close(fig);
    report=struct('passed',true,'sigma',2.38,'N',[16 32],'T',[4 8 16],'R',4, ...
        'checks',{{'serial/parallel exact estimates','deterministic bootstrap', ...
        'completed-result resume','sigma mismatch rejected','MAT CSV exports','PDF FIG exports'}});
    disp(report);
end

function close_owned_pool(pool,owned)
    if owned && isvalid(pool),delete(pool);end
end
