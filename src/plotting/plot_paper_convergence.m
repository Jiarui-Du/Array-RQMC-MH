function output=plot_paper_convergence(selection,options)
% Paper figures: selected Bayesian test; 2-by-3 MWG/joint Gaussian projection.
%
% Returns figure handles, selected summary data and exported filenames.
% All dataset paths may be overridden; empty paths select complete matching results under results/.
% Saving exports PDF/FIG and does not alter the saved simulation results.

    %% Declare figure selection, sources and style
    arguments
        selection (1,:) char = 'bayesian'
        options.NExponents (1,:) double = 10:18
        options.ShowFigures (1,1) logical = true
        options.SaveFigures (1,1) logical = true
        options.OutputDir (1,:) char = ''
        options.MwgMatFile (1,:) char = ''
        options.GaussianMatFile (1,:) char = ''
        options.RipleyMatFile (1,:) char = ''
        options.PimaMatFile (1,:) char = ''
        options.LineWidth (1,1) double = 1.8
        options.MarkerSize (1,1) double = 6
        options.AxisFontSize (1,1) double = 11
        options.LabelFontSize (1,1) double = 16
        options.LegendFontSize (1,1) double = 13
        options.BayesianTest (1,:) char = 'projection'
        options.TailCount (1,1) double = 5
    end

    %% Resolve result files and output folder
    assert(ismember(selection,{'bayesian','ripley','pima','gaussian','all'}),'Choose bayesian, ripley, pima, gaussian, or all');
    files=struct();root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
    if ismember(selection,{'gaussian','all'})
        files.gaussian_mwg=select_result_file('gaussian_mwg',[2 4 8],options.NExponents,options.MwgMatFile);
        files.gaussian=select_result_file('gaussian',[2 4 8],options.NExponents,options.GaussianMatFile);
    end
    if ismember(selection,{'bayesian','ripley','all'})
        files.ripley=select_result_file('ripley',[],options.NExponents,options.RipleyMatFile);
    end
    if ismember(selection,{'bayesian','pima','all'})
        files.pima=select_result_file('pima',[],options.NExponents,options.PimaMatFile);
    end
    folder=options.OutputDir;if isempty(folder),folder=fullfile(root,'figures');end
    groups={selection};if strcmp(selection,'all'),groups={'bayesian','gaussian'};end
    output=struct();options.ShowSlope=false;

    %% Render each requested dataset group
    for k=1:numel(groups)
        group=groups{k};panels={};labels={};
        if strcmp(group,'gaussian')
            item=draw_gaussian_comparison(files.gaussian_mwg,files.gaussian,options);
            if options.SaveFigures
                if ~isfolder(folder),mkdir(folder);end
                target=fullfile(folder,'gaussian_mwg353_joint983_projection_p2_p4_p8');
                savefig(item.figure,[target '.fig']);
                exportgraphics(item.figure,[target '.pdf'],'ContentType','vector');
                item.files={[target '.fig'],[target '.pdf']};
                fprintf('Saved Gaussian projection figure: %s.pdf\n',target);
            end
            output.gaussian=item;
            continue
        end

        %% Assemble only the requested Bayesian dataset panels
        test=options.BayesianTest;
        if strcmp(group,'bayesian')
            datasets={'ripley','pima'};labels={'Ripley','Pima'};
            % Include the actual test in the combined figure filename.
            stem=['bayesian_' test];
        else
            datasets={group};labels={[upper(group(1)) group(2:end)]};
            stem=['bayesian_' group '_' test];
        end
        for dataset=datasets
            raw=read_convergence_statistics(files.(dataset{1}),'NExponents',options.NExponents);
            panel=raw.summary(strcmp(raw.summary.test,test),:);
            assert(~isempty(panel),'Test %s is unavailable in %s',test,dataset{1});
            panels{end+1}=panel;
        end
        visible='off';if options.ShowFigures,visible='on';end
        fig=figure('Visible',visible,'Color','w','Name',[group ' / ' test],...
            'Position',[60 100 round(1520*numel(panels)/3) 480]);
        layout=tiledlayout(fig,1,numel(panels),'TileSpacing','compact','Padding','compact');
        for j=1:numel(panels)
            ax=nexttile(layout,j);draw_variance_panel(ax,panels{j},test,options);
            xlabel(ax,'$N$','Interpreter','latex','FontSize',options.LabelFontSize);
            if j==1,ylabel(ax,'Sample variance','FontSize',options.LabelFontSize);end
            title(ax,labels{j},'Interpreter','latex','FontSize',options.LabelFontSize);
        end

        %% Record figure data and optionally export artifacts
        item=struct('figure',fig,'summary',vertcat(panels{:}),'files',{{}});
        if options.SaveFigures
            if ~isfolder(folder),mkdir(folder);end
            target=fullfile(folder,stem);
            savefig(fig,[target '.fig']);
            exportgraphics(fig,[target '.pdf'],'ContentType','vector');
            item.files={[target '.fig'],[target '.pdf']};
            fprintf('Saved %s figure (test=%s): %s.pdf\n',group,test,target);
        end
        output.(group)=item;
    end
end
