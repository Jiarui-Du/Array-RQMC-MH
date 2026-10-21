function draw_variance_panel(ax,S,test,options)
% Shared empirical curves and reference lines for exploratory and paper plots.
%
% S contains one dimension and paired method/test/N variance rows. Options
% can override kernel labels and the reference exponent. Reference lines are
% visual guides anchored to observed values, not fitted or verified rate claims.

    %% Select matching variance series and empirical slopes
    d=unique(S.d);assert(isscalar(d));requestedN=unique(S.N)';hold(ax,'on');
                iidMethod='iid';arrayMethod='hilbert';
                if isfield(options,'IIDMethod'),iidMethod=options.IIDMethod;end
                if isfield(options,'ArrayMethod'),arrayMethod=options.ArrayMethod;end
                vi=values(S,test,iidMethod,requestedN);
                vh=values(S,test,arrayMethod,requestedN);
                slopeI=tail_slope(requestedN,vi,options.TailCount);
                slopeH=tail_slope(requestedN,vh,options.TailCount);
                labelI='IID';labelH='Array-RQMC';
                if isfield(options,'IIDLabel'),labelI=options.IIDLabel;end
                if isfield(options,'ArrayLabel'),labelH=options.ArrayLabel;end
                if ~isfield(options,'ShowSlope') || options.ShowSlope
                    labelI=sprintf('IID, slope %.3f',slopeI);
                    labelH=sprintf('Array-RQMC, slope %.3f',slopeH);
                end

                %% Draw data curves
                loglog(ax,requestedN,vi,'-o','Color',[0 .447 .741],...
                    'LineWidth',options.LineWidth,'MarkerSize',options.MarkerSize,...
                    'DisplayName',labelI);
                loglog(ax,requestedN,vh,'-s','Color',[.85 .325 .098],...
                    'LineWidth',options.LineWidth,'MarkerSize',options.MarkerSize,...
                    'DisplayName',labelH);

                %% Add anchored reference lines and format the axes
                % Match the 1d reference anchors: first IID point, middle
                % Hilbert point. Only the dimension-dependent exponent differs.
                loglog(ax,requestedN,vi(1)*(requestedN/requestedN(1)).^-1,'--',...
                    'Color',[.929 .694 .125],'LineWidth',1.1,'DisplayName','$N^{-1}$');
                anchor=max(1,ceil(numel(requestedN)/2));
                exponent=-1-1/(d+2);referenceLabel=sprintf('$N^{-1-1/%d}$',d+2);
                if isfield(options,'ReferenceExponent'),exponent=options.ReferenceExponent;end
                if isfield(options,'ReferenceLabel'),referenceLabel=options.ReferenceLabel;end
                ref=vh(anchor)*(requestedN/requestedN(anchor)).^exponent;
                loglog(ax,requestedN,ref,':','Color',[.494 .184 .556],'LineWidth',1.3,...
                    'DisplayName',referenceLabel);
                set(ax,'XScale','log','YScale','log','FontSize',options.AxisFontSize);
                if numel(requestedN)>1,xlim(ax,[requestedN(1),requestedN(end)]);end
                grid(ax,'on');box(ax,'on');
                legend(ax,'Location','southwest','Interpreter','latex','FontSize',options.LegendFontSize);
end

%% Local helper: values
% Align one test/method series to requested N and reject missing or nonpositive variances.
function v=values(S,test,method,N)
    rows=S(strcmp(S.test,test) & strcmp(S.method,method),:);
    assert(height(rows)==numel(N),'Missing or duplicate variance rows');
    [ok,idx]=ismember(N,rows.N);assert(all(ok),'Missing N in statistics');
    v=rows.variance(idx);
    assert(all(isfinite(v) & v>0),'Log plots require positive finite variances');
end

%% Local helper: tail_slope
% Fit the last requested number of positive finite variance points.
function slope=tail_slope(N,v,count)
    valid=find(isfinite(v) & v>0);
    valid=valid(max(1,numel(valid)-count+1):end);
    slope=NaN;if numel(valid)<2,return;end
    fitted=polyfit(log(N(valid(:))),log(v(valid(:))),1);slope=fitted(1);
end
