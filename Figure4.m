% Figure 4: archived N=80 staggered-grid data; original plotting style.
clc; close all; clear;
here = fileparts(mfilename('fullpath'));
out_dir = fullfile(here,'output');
if ~exist(out_dir,'dir'), mkdir(out_dir); end
S = load(fullfile(here,'data','Figure4_data.mat'));

Ev = S.Ev(:)';
ell = S.ell_star(:)';
mz2 = S.mz2_velocity(:)';
font_name = 'Times New Roman';
ell_symbol = char(8467); % Unicode script ell; MATLAB TeX does not support \ell.

fig = figure('Color','w','Units','centimeters','Position',[3 3 25.5 10.8]);
set(fig,'DefaultAxesFontName',font_name,'DefaultTextFontName',font_name, ...
    'DefaultTextInterpreter','tex','DefaultAxesTickLabelInterpreter','tex');
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% (a) Energy budget
ax1 = nexttile;
hold(ax1,'on');
plot(ax1,Ev,S.shear_rate(:)','o-','LineWidth',1.7,'MarkerSize',4,'DisplayName','Shear production P_S');
plot(ax1,Ev,S.front_rate(:)','s-','LineWidth',1.7,'MarkerSize',4,'DisplayName','Buoyancy production P_F');
plot(ax1,Ev,S.viscous_rate(:)','^-','LineWidth',1.7,'MarkerSize',4,'DisplayName','Viscous loss -\epsilon_v');
plot(ax1,Ev,S.diffusive_rate(:)','v-','LineWidth',1.7,'MarkerSize',4,'DisplayName','Diffusive loss -\epsilon_b');
plot(ax1,Ev,S.sigma_star(:)','D--','LineWidth',1.6,'MarkerSize',4,'DisplayName','Net growth Re(\sigma)');
yline(ax1,0,'-','LineWidth',0.8,'HandleVisibility','off');
set(ax1,'XScale','log','tickdir','out','LineWidth',0.8);
xlabel(ax1,'Vertical Ekman number E_v');
ylabel(ax1,'Energy-budget rate / (2E)');
legend(ax1,'Location','best','Box','off', ...
    'FontName',font_name,'Interpreter','tex');
text(ax1,-0.08,1.03,'a','Units','normalized','FontWeight','bold','FontSize',12);
box(ax1,'on');

% (b) Vertical-gradient mechanism
ax2 = nexttile;
hold(ax2,'on');
plot(ax2,ell,mz2,'o-','LineWidth',1.8,'MarkerSize',5,'DisplayName','Fastest-growing SI modes');

[~,iref] = min(abs(Ev-1e-4));
c2 = mz2(iref)/ell(iref)^2;
xref = logspace(log10(min(ell)),log10(max(ell)),200);
plot(ax2,xref,c2*xref.^2,'--','LineWidth',1.4,'DisplayName',['m_{eff}^2 ' char(8733) ' ' ell_symbol '_*^2']);

targets = [1e-5 1e-4 1e-3 3e-3];
labels = {'10^{-5}','10^{-4}','10^{-3}','3\times10^{-3}'};
for j=1:numel(targets)
    [~,ii] = min(abs(Ev-targets(j)));
    x = ell(ii)*1.03; y = mz2(ii)*1.05; align = 'left';
    if ell(ii)==max(ell) % Keep the upper-right label inside the panel.
        x = ell(ii)/1.03; y = mz2(ii)/1.05; align = 'right';
    end
    text(ax2,x,y,labels{j},'FontSize',8,'FontName',font_name, ...
        'Interpreter','tex','HorizontalAlignment',align);
end

set(ax2,'XScale','log','YScale','log','tickdir','out','LineWidth',0.8);
xlabel(ax2,['Fastest-growing wavenumber ' ell_symbol '_*']);
ylabel(ax2,'Velocity-gradient metric m_{eff}^2');
legend(ax2,'Location','northwest','Box','off', ...
    'FontName',font_name,'Interpreter','tex');
text(ax2,-0.08,1.03,'b','Units','normalized','FontWeight','bold','FontSize',12);
box(ax2,'on');


set(findall(fig,'-property','FontName'),'FontName',font_name);
drawnow;
exportgraphics(fig,fullfile(out_dir,'Figure4_N80.pdf'), ...
    'ContentType','vector','BackgroundColor','white');
exportgraphics(fig,fullfile(out_dir,'Figure4_N80.tiff'), ...
    'Resolution',600,'BackgroundColor','white');
savefig(fig,fullfile(out_dir,'Figure4_N80.fig'));
