% Figure 3: N=80 staggered-grid peaks; original plotting style retained.
clc; close all; clear;
here = fileparts(mfilename('fullpath'));
out_dir = fullfile(here,'output');
if ~exist(out_dir,'dir'), mkdir(out_dir); end
S = load(fullfile(here,'data','Figure3_data.mat'));
Ri = S.Ri(:)'; Ev = S.Ev(:);
mask = logical(S.unstable) & ~logical(S.scan_boundary);
G = S.sigma_star;
Y = S.lambda_over_L ./ Ri.^S.fit_Ri_exponent;
G(~mask) = NaN; Y(~mask) = NaN;
Et = S.Ev_fit(:); Yfit = S.fit_prefactor*Et.^S.fit_Ev_exponent;
ms = S.marker_size; mew = S.marker_edge_width;
fig = figure('Color','w','Units','inches','Position',[1 1 14.8 5.8], ...
    'Renderer','painters','Name','Figure 3: hollow-circle markers','NumberTitle','off');
ax1 = axes('Parent',fig,'Position',[.07 .16 .40 .76]); hold(ax1,'on');
ax2 = axes('Parent',fig,'Position',[.57 .16 .40 .76]); hold(ax2,'on');
set([ax1 ax2],'FontName','Times New Roman','FontSize',14,'Box','on', ...
    'XScale','log','XLim',[min(Ev)*.9,max(Ev)*1.1], ...
    'XGrid','on','YGrid','on','XMinorGrid','on','YMinorGrid','on', ...
    'GridLineStyle',':','MinorGridLineStyle',':','GridAlpha',.20,'MinorGridAlpha',.15, ...
    'Layer','bottom','linewidth',1.2);
set(ax1,'YLim',[0 1.82],'YTick',0:.25:1.75,'tickdir','out');
set(ax2,'YScale','log','tickdir','out');

handles = gobjects(numel(Ri),1);
for j=1:numel(Ri)
    c=S.rgb(j,:);
    plot(ax1,Ev,G(:,j),'-','Color',c,'LineWidth',1.9,'HandleVisibility','off');
    plot(ax1,Ev,G(:,j),'LineStyle','none','Marker','o','MarkerSize',ms, ...
        'MarkerFaceColor','none','MarkerEdgeColor',c,'Color',c, ...
        'LineWidth',mew,'HandleVisibility','off');
    handles(j)=plot(ax1,NaN,NaN,'-o','Color',c,'LineWidth',1.5, ...
        'MarkerSize',ms,'MarkerFaceColor','none','MarkerEdgeColor',c, ...
        'DisplayName',sprintf('Ri=%.2f',Ri(j)));
end
plot(ax2,Et,Yfit,'--','Color',[0 0 0],'LineWidth',2.0,'HandleVisibility','off');
for j=1:numel(Ri)
    c=S.rgb(j,:);
    plot(ax2,Ev,Y(:,j),'LineStyle','none','Marker','o','MarkerSize',ms, ...
        'MarkerFaceColor','none','MarkerEdgeColor',c,'Color',c, ...
        'LineWidth',mew,'HandleVisibility','off');
end
vals=[Y(isfinite(Y));Yfit];
ylim(ax2,[min(vals)*.8,max(vals)*1.35]);
xlabel(ax1,'$E_v$','Interpreter','latex','FontSize',16);
xlabel(ax2,'$E_v$','Interpreter','latex','FontSize',16);
ylabel(ax1,'Growth rate $\mathrm{Re}(\sigma)$','Interpreter','latex','FontSize',16);
ylabel(ax2,'$\lambda_*/(L\,Ri^{0.83})$','Interpreter','latex','FontSize',16);
legend(ax1,handles,'Location','southwest','NumColumns',2,'Box','off', ...
    'FontSize',12,'Interpreter','none');
for j=1:2
    if j==1, ax=ax1; tag='a'; else, ax=ax2; tag='b'; end
    text(ax,-.10,1.01,tag,'Units','normalized','FontSize',22,'FontWeight','bold', ...
        'FontName','Times New Roman','Clipping','off');
end
text(ax2,.03,.96,'Open circles: numerical peaks','Units','normalized', ...
    'VerticalAlignment','top','FontSize',13,'FontName','Times New Roman','Interpreter','none');
text(ax2,.03,.885, ...
    'Dashed line: $\lambda_*/(L\,Ri^{0.83})\approx3.33\,E_v^{1/4}$', ...
    'Units','normalized','VerticalAlignment','top','FontSize',13,'Interpreter','latex');

drawnow;
exportgraphics(fig,fullfile(out_dir,'Figure3_N80.pdf'), ...
    'ContentType','vector','BackgroundColor','white');
exportgraphics(fig,fullfile(out_dir,'Figure3_N80.tiff'), ...
    'Resolution',600,'BackgroundColor','white');
savefig(fig,fullfile(out_dir,'Figure3_N80.fig'));
