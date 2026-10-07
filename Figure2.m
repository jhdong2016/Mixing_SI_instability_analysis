% Figure 2: growth spectra (a) and fastest-growing SI modes (b-e).
clc; close all; clear;
here = fileparts(mfilename('fullpath'));
out_dir = fullfile(here,'output');
if ~exist(out_dir,'dir'), mkdir(out_dir); end
S = load(fullfile(here,'data','Figure2_data.mat'));
M = load(fullfile(here,'data','Figure2_SI_vertical_modes_data.mat'));
font_name = 'Times New Roman';
Ev_selected = [1e-5 1e-4 1e-3 3e-3];
labels = {'10^{-5}','10^{-4}','10^{-3}','3\times10^{-3}'};
letters = 'bcde';
[lambda,order] = sort(2*pi./S.ell(:),'ascend');
[Y,Z] = meshgrid(M.y,M.z);
B = Z-Y/M.Ri;

fig = figure('Color','w','Units','centimeters','Position',[3 3 25.5 19.5]);
set(fig,'DefaultAxesFontName',font_name,'DefaultTextFontName',font_name, ...
    'DefaultTextInterpreter','tex','DefaultAxesTickLabelInterpreter','tex');

% (a) Original Figure 2, occupying the first row.
ax = axes(fig,'Position',[0.075 0.550 0.81 0.34]);
hold(ax,'on');
plot(ax,lambda,S.inviscid(order),'--','LineWidth',2.0,'DisplayName','Inviscid');
for j = 1:4
    [~,idx] = min(abs(S.Ev(:)-Ev_selected(j)));
    plot(ax,lambda,S.growth(idx,order),'LineWidth',1.8, ...
        'DisplayName',sprintf('E_v = %.0e',Ev_selected(j)));
end
set(ax,'XScale','log','FontSize',8,'LineWidth',0.65,'TickDir','out','Box','on','FontSize',10);
xlim(ax,[min(lambda) max(lambda)]); ylim(ax,[-0.01 0.68]);
xlabel(ax,'Wavelength \lambda/L','FontSize',10);
ylabel(ax,'Growth rate Re(\sigma)','FontSize',10);
% text(ax,0.65,0.96,'Ri = 0.7,  \delta = 0.1,  Pr = 1', ...
%     'Units','normalized','VerticalAlignment','top');
legend(ax,'Location','southeast','Box','off','NumColumns',2);
text(ax,-0.02,1.03,'a','Units','normalized','FontWeight','bold','FontSize',15);

% (b-e) Four modes on the second row, with one shared colorbar.
colormap(fig,M.cmap);
for j = 1:4
    ax = axes(fig,'Position',[0.075+(j-1)*0.212 0.09 0.174 0.35]);
    surf(ax,M.y,M.z,zeros(size(M.W(:,:,j))),M.W(:,:,j)); view(ax,2);
    shading(ax,'interp'); hold(ax,'on');
    contour(ax,M.y,M.z,B,-1.6:0.2:0.6,'LineColor',[0.72 0.72 0.72], ...
        'LineWidth',0.55,'LineStyle','-');
    if j ==1
    set(ax,'FontSize',8,'LineWidth',0.65,'Box','on','Layer','top', ...
        'TickDir','out','CLim',[-1 1],'XGrid','off','YGrid','off','fontsize',10);
    else
        set(ax,'FontSize',8,'LineWidth',0.65,'Box','on','Layer','top', ...
            'TickDir','out','CLim',[-1 1],'XGrid','off','YGrid','off','fontsize',10,'yticklabel','');
    end
    xlim(ax,[-0.5 0.5]); ylim(ax,[-1 0]);
    xticks(ax,[-0.4 0 0.4]); yticks(ax,-1:0.2:0);
    xlabel(ax,'y/L','FontSize',10);
    if j==1, ylabel(ax,'z/H','FontSize',10); end
    title(ax,['E_v = ' labels{j}],'FontSize',10,'FontWeight','normal');
    text(ax,-0.15,1.03,letters(j),'Units','normalized', ...
        'FontWeight','bold','FontSize',15);
end
cb = colorbar(ax,'Location','eastoutside');
set(cb,'Position',[0.912 0.09 0.012 0.35],'Ticks',-1:0.5:1, ...
    'FontName',font_name,'FontSize',10,'LineWidth',0.65);
set(ax,'Position',[0.711 0.09 0.174 0.35]);
cb.Label.String = 'Normalized Re(w)';
cb.Label.Interpreter = 'tex'; cb.Label.FontSize = 9;
set(findall(fig,'-property','FontName'),'FontName',font_name);
drawnow;
exportgraphics(fig,fullfile(out_dir,'Figure2_N80.pdf'), ...
    'Resolution',600);
exportgraphics(fig,fullfile(out_dir,'Figure2_N80.tiff'), ...
    'Resolution',600);
savefig(fig,fullfile(out_dir,'Figure2_N80.fig'));
