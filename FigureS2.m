clc; clear; close all;
% Replot archived numerical data; no eigenvalue solves are performed.
here = fileparts(mfilename('fullpath'));
outputDir = fullfile(here,'output');
if ~exist(outputDir,'dir'), mkdir(outputDir); end
S = load(fullfile(here,'data','FigureS2_data.mat'));
delta=S.delta(:)'; Ev=S.Ev(:); mask=logical(S.unstable);
G=S.sigma_star; W=S.lambda_over_L; G(~mask)=NaN; W(~mask)=NaN;
beta=zeros(numel(delta),1); log10c=beta;
for j=1:numel(delta)
    jj=mask(:,j);
    p=polyfit(log10(Ev(jj)),log10(W(jj,j)),1);
    beta(j)=p(1); log10c(j)=p(2);
end
Et=S.reference_Ev(:);
Wt=S.reference_anchor_lambda_over_L*(Et/S.reference_anchor_Ev).^S.reference_exponent;

fig=figure('Color','w','Units','inches','Position',[1 1 14 6.06944], ...
    'Renderer','painters','Name','Figure S2: delta sensitivity','NumberTitle','off','DefaultTextInterpreter','tex');
ax1=axes('Parent',fig,'Position',[.0675 .1610526 .415 .7010526]); hold(ax1,'on');
ax2=axes('Parent',fig,'Position',[.5675 .1610526 .415 .7010526]); hold(ax2,'on');
set([ax1 ax2],'FontName','Times New Roman','FontSize',14,'Box','on', ...
    'LineWidth',.9,'XScale','log','XLim',[9e-6 1.12e-2], ...
    'XTick',[1e-5 1e-4 1e-3 1e-2], ...
    'XGrid','on','YGrid','on','XMinorGrid','on','YMinorGrid','off', ...
    'GridLineStyle',':','MinorGridLineStyle',':','GridAlpha',.16,'MinorGridAlpha',.12);
set(ax1,'YScale','linear','YLim',[-.025 .71],'YTick',0:.1:.7,'TickDir','out');
set(ax2,'YScale','log','YLim',[.112 .90],'YTick',[.15 .2 .3 .5 .8], ...
    'YTickLabel',{'0.15','0.2','0.3','0.5','0.8'},'TickDir','out');
styles={'-','--','-.'}; h1=gobjects(4,1); h2=gobjects(4,1);
for j=1:numel(delta)
    h1(j)=plot(ax1,Ev,G(:,j),'Color',S.rgb(j,:),'LineStyle',styles{j}, ...
        'LineWidth',1.5,'Marker','o','MarkerSize',6,'MarkerFaceColor','none', ...
        'DisplayName',sprintf('\\delta=%.2f',delta(j)));
end
% A shared zero-cross symbol denotes no detected growth for all delta values.
stableAll=all(~mask,2);
h1(4)=plot(ax1,Ev(stableAll),zeros(nnz(stableAll),1),'x','LineStyle','none', ...
    'Color',[.2 .2 .2],'MarkerSize',5.5,'LineWidth',1.2,'DisplayName','No detected growth');
legend(ax1,h1,'Location','southwest','Box','off','FontSize',13,'Interpreter','tex');
h2(4)=plot(ax2,Et,Wt,'--','Color',[0 0 0],'LineWidth',1.4, ...
    'DisplayName','E_v^{1/4} reference');
for j=1:numel(delta)
    h2(j)=plot(ax2,Ev,W(:,j),'LineStyle','none','Marker','o', ...
        'MarkerSize',6,'MarkerFaceColor','none','MarkerEdgeColor',S.rgb(j,:), ...
        'Color',S.rgb(j,:),'LineWidth',1.15, ...
        'DisplayName',sprintf('\\delta=%.2f, \\beta=%.3f',delta(j),beta(j)));
end
lg=legend(ax2,h2,'Location','northwest','Box','off','FontSize',12.4,'Interpreter','tex');
title(lg,'Log-log slope \beta','Interpreter','tex','FontSize',12.5);
xlabel(ax1,'E_v','Interpreter','tex','FontSize',17);
xlabel(ax2,'E_v','Interpreter','tex','FontSize',17);
ylabel(ax1,'Peak growth rate Re(\sigma_*)','Interpreter','tex','FontSize',17);
ylabel(ax2,'Selected wavelength \lambda_*/L','Interpreter','tex','FontSize',17);
text(ax1,-.105,1.025,'a','Units','normalized','FontName','Times New Roman', ...
    'FontSize',23,'FontWeight','bold','VerticalAlignment','bottom','Clipping','off');
text(ax2,-.105,1.025,'b','Units','normalized','FontName','Times New Roman', ...
    'FontSize',23,'FontWeight','bold','VerticalAlignment','bottom','Clipping','off');
name = 'FigureS2';
set(findall(fig,'-property','FontName'),'FontName','Times New Roman');
drawnow;
exportgraphics(fig,fullfile(outputDir,[name '.pdf']),'ContentType','vector');
exportgraphics(fig,fullfile(outputDir,[name '_600dpi.tif']),'Resolution',600);
savefig(fig,fullfile(outputDir,[name '.fig']));
fprintf('Saved to %s\n',outputDir);
