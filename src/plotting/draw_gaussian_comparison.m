function item=draw_gaussian_comparison(mwgPath,jointPath,options)
% Common linear projection: MWG in row 1, joint Hilbert RWM in row 2.
% Keep the two IID baselines and reference exponents specific to each kernel.
%
% Inputs are saved MWG and joint-RWM paths plus the parent plotting options.
% Checks ensure matching observation settings and standard Gaussian cases.
% The returned item retains the source tables used for the six panels.

    %% Load and verify comparable experiments
    saved=load(mwgPath,'results');mwg=saved.results;
    saved=load(jointPath,'results');joint=saved.results;
    for key={'T','BurnIn','Replications','Initialization'}
        assert(isequal(mwg.settings.(key{1}),joint.settings.(key{1})), ...
            'Gaussian comparison setting differs: %s',key{1});
    end
    for source={mwg,joint}
        for p=[2 4 8]
            c=source{1}.cases([source{1}.cases.dimension]==p);
            assert(isscalar(c) && isequal(c.mu,zeros(1,p)) && isequal(c.Sigma,eye(p)), ...
                'Comparison requires standard Gaussian cases p=2,4,8');
        end
    end

    %% Select the common linear projection and normalize table columns
    [M,F]=summarize_gaussian_mwg(mwg,'Dimensions',[2 4 8], ...
        'NExponents',options.NExponents,'TailCount',options.TailCount,'Bootstrap',0);
    H=read_convergence_statistics(jointPath,'Dimensions',[2 4 8], ...
        'NExponents',options.NExponents);
    M=M(strcmp(M.test,'projection'),:);H=H.summary(strcmp(H.summary.test,'projection'),:);
    S=M(:,{'p','N','test','method','variance'});S.Properties.VariableNames{1}='d';
    S.kernel=repmat({'mwg'},height(S),1);
    J=H(:,{'d','N','test','method','variance'});J.kernel=repmat({'joint_hilbert'},height(J),1);

    %% Lay out two kernels across three dimensions
    visible='off';if options.ShowFigures,visible='on';end
    fig=figure('Visible',visible,'Color','w','Name','Gaussian MWG and joint Hilbert', ...
        'Position',[60 60 1520 960]);
    layout=tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
    dims=[2 4 8];
    for row=1:2
        style=options;style.ShowSlope=false;
        if row==1
            source=S;rowLabel='MWG';
            style.IIDMethod='iid_mwg';style.ArrayMethod='array_mwg';
            style.IIDLabel='IID-MWG';style.ArrayLabel='Array-RQMC-MWG';
            style.ReferenceExponent=-4/3;style.ReferenceLabel='$N^{-4/3}$';
        else
            source=J;rowLabel='Joint Hilbert';
            style.IIDLabel='IID-RWM';style.ArrayLabel='Array-RQMC-Hilbert';
        end
        for col=1:3
            ax=nexttile(layout,(row-1)*3+col);
            draw_variance_panel(ax,source(source.d==dims(col),:),'projection',style);
            title(ax,sprintf('$p=%d$',dims(col)),'Interpreter','latex','FontSize',options.LabelFontSize);
            xlabel(ax,'$N$','Interpreter','latex','FontSize',options.LabelFontSize);
            if col==1,ylabel(ax,{rowLabel,'Sample variance'},'FontSize',options.LabelFontSize);end
        end
    end

    %% Return the figure and the data behind it
    item=struct('figure',fig,'summary',[S;J],'mwgSummary',M,'jointSummary',H, ...
        'mwgSlopes',F(strcmp(F.test,'projection'),:), ...
        'sources',{{mwgPath,jointPath}},'files',{{}});
end
