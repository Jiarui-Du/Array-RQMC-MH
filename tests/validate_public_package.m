function report=validate_public_package
% Bounded release check: saved-result readers plus N<=32, T=8, R=4 runs.
% All newly generated outputs go to a temporary folder, never bundled results.
    root=fileparts(fileparts(mfilename('fullpath')));addpath(genpath(fullfile(root,'src')));
    temp=tempname;mkdir(temp);cleanup=onCleanup(@()rmdir(temp,'s'));
    %% Check every bundled main result using the public readers
    for ds={'gaussian','gaussian_mwg','ripley','pima'}
        v=view_variance_reduction(ds{1},'Dimension',[2 4 8],'NExponents',18, ...
            'Tests','projection','Bootstrap',20,'Display',false);
        assert(all(isfinite(v.VRF)) && all(isfinite(v.TAE)));
    end
    s=export_array_rqmc_vrf_tae_slopes_v2(fullfile(root,'results','one_dim'), ...
        'TargetN',2^18,'RequiredR',200,'TailCount',5,'Strict',true, ...
        'OutputCSV',fullfile(temp,'one_dim.csv'),'OutputTeX',fullfile(temp,'one_dim.tex'));
    assert(~isempty(s));
    assert(isempty(dir(fullfile(root,'results','**','*.csv'))));
    %% Small statistical and independent timing runs
    common={'NExponents',3:5,'Replications',4,'T',8,'BurnIn',2, ...
        'TailCount',3,'Bootstrap',20,'RunTiming',true,'TimingReplications',2, ...
        'TimingWarmupExponent',3,'UseParallel',false};
    run_gaussian_hilbert_experiment(common{:},'Dimensions',[2 4 8],'Target','standard', ...
        'Initialization','zero','OutputDir',fullfile(temp,'joint'));
    run_gaussian_mwg_experiment(common{:},'Dimensions',[2 4 8],'OutputDir',fullfile(temp,'mwg'));
    for ds={'ripley','pima'}
        b=run_bayesian_single_parfor(common{:},'Dataset',ds{1},'OutputDir',fullfile(temp,ds{1}));
        assert(isequal(b.testNames,{'projection'}));
    end
    run_array_rqmc_suite('TargetIDs',{'normal','double_well_e1','gaussian_mixture3','student_t3'},'TestFunctionIDs',{'cdf','x','x2'}, ...
        'mVec',3:5,'T0',2,'T',8,'R',4,'UseParallel',false,'AutoStartPool',false, ...
        'RunTiming',true,'TimingM',5,'TimingR',2,'TimingWarmupM',3,'Bootstrap',20, ...
        'OutputDir',fullfile(temp,'one_dim'),'ContinueOnError',false);
    %% Regenerate every main figure from saved results without changing originals
    p=plot_paper_convergence('all','ShowFigures',false,'OutputDir',temp);
    close(p.gaussian.figure);close(p.bayesian.figure);
    for ds={'ripley','pima'}
        q=plot_paper_convergence(ds{1},'ShowFigures',false,'OutputDir',temp);close(q.(ds{1}).figure);
    end
    plot_array_rqmc_suite_results(fullfile(root,'results','one_dim'), ...
        'RequiredR',200,'FigureVisible','off','SaveFigures',true,'OutputDir',temp, ...
        'MetadataOutputDir',temp);close all;
    r=load_array_rqmc_horizon_result('');
    [fig,~]=plot_array_rqmc_chain_length_scaling_results(r,fullfile(temp,'horizon'),'Visible','off');close(fig);
    report=struct('passed',true,'matlab',version,'N',2.^(3:5),'T',8,'R',4, ...
        'formalRunStarted',false,'checks',{{'bundled result readers','one-dimensional run', ...
        'Gaussian joint/MWG runs','Ripley/Pima projection-only runs','independent timing', ...
        'main and single-dataset figure exports','horizon figure export'}});
    disp(report);
end
