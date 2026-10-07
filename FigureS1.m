clc; clear; close all;
% Replot archived numerical data; no eigenvalue solves are performed.
here = fileparts(mfilename('fullpath'));
outputDir = fullfile(here,'output');
if ~exist(outputDir,'dir'), mkdir(outputDir); end
S = load(fullfile(here,'data','FigureS1_data.mat'));
ri = S.Ri_values(:); gri = S.grid_Ri(:); N = S.grid_N(:);
Eg = 100*abs(bsxfun(@rdivide,S.grid_sigma,S.grid_reference_sigma(:).')-1);
El = 100*abs(bsxfun(@rdivide,S.grid_lambda,S.grid_reference_lambda(:).')-1);
Pg = 100*abs(S.pair_sigma_fd(:)./S.pair_sigma_cheb(:)-1);
Pl = 100*abs(S.pair_lambda_fd(:)./S.pair_lambda_cheb(:)-1);
fprintf('Paired cases: %d; Ev = %.3g to %.3g\n',numel(Pg),min(S.pair_Ev),max(S.pair_Ev));
fprintf('Maximum staggered N=%d / Chebyshev%d differences: growth %.6f%%; wavelength %.6f%%\n',S.pair_fd_N,S.pair_cheb_N,max(Pg),max(Pl));

fig = figure('Color','w','Units','inches','Position',[1 1 13.2 9.6], ...
    'Renderer','painters','Name','Figure S1: archived numerical sensitivity', ...
    'NumberTitle','off','DefaultTextInterpreter','tex');
positions = [.080 .585 .400 .355; .580 .585 .400 .355; ...
             .080 .085 .400 .355; .580 .085 .400 .355];
ax = gobjects(4,1);
for j=1:4
    ax(j)=axes('Parent',fig,'Position',positions(j,:)); hold(ax(j),'on');
    set(ax(j),'FontName','Times New Roman','FontSize',14,'Box','on', ...
        'LineWidth',.9,'TickDir','out','XGrid','on','YGrid','on', ...
        'GridLineStyle',':','GridAlpha',.20,'Layer','bottom');
end
lw=1.6; ms=6.8;
for j=1:numel(gri)
    [~,ic]=min(abs(ri-gri(j))); color=S.rgb(ic,:);
    plot(ax(1),N,Eg(:,j),'-o','Color',color,'LineWidth',lw,'MarkerSize',ms, ...
        'MarkerFaceColor','none','DisplayName',sprintf('Ri = %.2f',gri(j)));
    plot(ax(2),N,El(:,j),'-o','Color',color,'LineWidth',lw,'MarkerSize',ms, ...
        'MarkerFaceColor','none','HandleVisibility','off');
end
for j=1:2
    set(ax(j),'XScale','log','YScale','linear','XTick',N,'XTickLabel',cellstr(num2str(N)), ...
        'XLim',[min(N)*.88 max(N)*1.12]);
    xlabel(ax(j),'Staggered vertical cells N','FontSize',17);
end
set(ax(1),'YLim',[0 .35],'YTick',0:.05:.35);
set(ax(2),'YLim',[-.2 12],'YTick',0:2:12);
ylabel(ax(1),'Peak-growth difference (%)','FontSize',17);
ylabel(ax(2),'Wavelength difference (%)','FontSize',17);
title(ax(1),'Vertical-resolution sensitivity','FontSize',17,'FontWeight','normal');
title(ax(2),'Vertical-resolution sensitivity','FontSize',17,'FontWeight','normal');
legend(ax(1),'show','Location','northeast','NumColumns',2,'Box','off','FontSize',12);
text(ax(1),.035,.95,sprintf('E_v = 10^{%g}; reference: Chebyshev%d',log10(S.grid_Ev),S.grid_reference_N), ...
    'Units','normalized','FontSize',11,'VerticalAlignment','top');
text(ax(2),.43,.95,sprintf('E_v = 10^{%g}; Chebyshev%d reference',log10(S.grid_Ev),S.grid_reference_N), ...
    'Units','normalized','FontSize',11,'VerticalAlignment','top');

limG=[0 max(2,ceil(max([S.pair_sigma_fd(:);S.pair_sigma_cheb(:)])/.4)*.4)];
limL=[0 max(.50,ceil(max([S.pair_lambda_fd(:);S.pair_lambda_cheb(:)])/.1)*.1)];
plot(ax(3),limG,limG,'k--','LineWidth',1.2,'HandleVisibility','off');
plot(ax(4),limL,limL,'k--','LineWidth',1.2,'HandleVisibility','off');
for j=1:numel(ri)
    use=abs(S.pair_Ri(:)-ri(j))<1e-10;
    plot(ax(3),S.pair_sigma_cheb(use),S.pair_sigma_fd(use),'o', ...
        'LineStyle','none','Color',S.rgb(j,:),'MarkerFaceColor','none', ...
        'MarkerSize',ms,'LineWidth',1.2,'DisplayName',sprintf('Ri = %.2f',ri(j)));
    plot(ax(4),S.pair_lambda_cheb(use),S.pair_lambda_fd(use),'o', ...
        'LineStyle','none','Color',S.rgb(j,:),'MarkerFaceColor','none', ...
        'MarkerSize',ms,'LineWidth',1.2,'HandleVisibility','off');
end
set(ax(3),'XLim',limG,'YLim',limG,'XTick',0:.4:limG(2),'YTick',0:.4:limG(2));
set(ax(4),'XLim',limL,'YLim',limL,'XTick',0:.1:limL(2),'YTick',0:.1:limL(2));
xlabel(ax(3),sprintf('Chebyshev%d: peak growth Re(\\sigma_*)',S.pair_cheb_N),'FontSize',16);
ylabel(ax(3),sprintf('Staggered N=%d: peak growth Re(\\sigma_*)',S.pair_fd_N),'FontSize',16);
xlabel(ax(4),sprintf('Chebyshev%d: \\lambda_*/L',S.pair_cheb_N),'FontSize',17);
ylabel(ax(4),sprintf('Staggered N=%d: \\lambda_*/L',S.pair_fd_N),'FontSize',17);
title(ax(3),'Solver comparison: peak growth','FontSize',17,'FontWeight','normal');
title(ax(4),'Solver comparison: wavelength','FontSize',17,'FontWeight','normal');
legend(ax(3),'show','Location','southeast','NumColumns',2,'Box','off','FontSize',11);
text(ax(3),.035,.95,sprintf('%d paired growing cases\nMaximum difference: %.3f%%',numel(Pg),max(Pg)), ...
    'Units','normalized','FontSize',12,'VerticalAlignment','top');
text(ax(4),.035,.95,sprintf('%d paired growing cases\nMaximum difference: %.2f%%',numel(Pl),max(Pl)), ...
    'Units','normalized','FontSize',12,'VerticalAlignment','top');
for j=1:4
    text(ax(j),-.135,1.035,char('a'+j-1),'Units','normalized', ...
        'FontName','Times New Roman','FontWeight','bold','FontSize',25);
end
name = 'FigureS1';
set(findall(fig,'-property','FontName'),'FontName','Times New Roman');
drawnow;
exportgraphics(fig,fullfile(outputDir,[name '.pdf']),'ContentType','vector');
exportgraphics(fig,fullfile(outputDir,[name '_600dpi.tif']),'Resolution',600);
savefig(fig,fullfile(outputDir,[name '.fig']));
fprintf('Saved to %s\n',outputDir);
