% Figure 1: N=80 staggered-grid results. Run this script directly.
clc; close all; clear;
here = fileparts(mfilename('fullpath'));
out_dir = fullfile(here,'output');
if ~exist(out_dir,'dir'), mkdir(out_dir); end
load(fullfile(here,'data','Figure1_data.mat'));
Ri_list = [0.25, 0.70, 0.95];
Ev_list = [0, 1e-3];
climits = [0, 1.7];
pdf_dpi = 600;
tiff_dpi = 300;
pdf_content_type = 'image';       % 'image' (seam-free) or 'vector'
show_si_axis = true;              % phi=0 dashed line in the reference
show_zero_contour = true;         % white zero-growth contour in the reference
font_name = 'Times New Roman';


letters = 'abcdef';
panels = cell(2,3);
[lambda,order] = sort(2*pi./alpha(:));
for row = 1:2
    for col = 1:3
        panels{row,col} = struct('Ri',Ri_list(col),'Ev',Ev_list(row), ...
            'phi',phi(:).','lambda',lambda,'growth',growth(order,:,row,col), ...
            'phi_edges',cell_edges(phi,false), ...
            'lambda_edges',cell_edges(lambda,true));
    end
end
%% Exact canvas proportions used by the reference figure (points)
PW = 2.48*72;
PH = 2.35*72;
top = 22;
row_gap = 24;
foot = 18;
W = 3*PW + 65;
H = top + 2*PH + row_gap + foot;

fig = figure('Color','w', 'Units','inches', ...
    'Position',[0.6,0.6,W/72,H/72], ...
    'Name','Figure 1: vertical mixing and SI', 'NumberTitle','off', ...
    'InvertHardcopy','off');
set(fig, 'PaperUnits','inches', 'PaperSize',[W,H]/72, ...
    'PaperPosition',[0,0,W,H]/72, 'PaperPositionMode','manual');
colormap(fig, cmap);
axes_handles = gobjects(2,3);

%% Six panels
for row = 1:2
    for col = 1:3
        ii = (row-1)*3 + col;
        p = panels{row,col};
        panel_left = (col-1)*PW;
        panel_top = top + (row-1)*(PH + row_gap);
        position = [(panel_left + 0.20*PW)/W, ...
            (H-panel_top-PH+0.18*PH)/H, 0.755*PW/W, 0.72*PH/H];
        ax = axes('Parent',fig, 'Units','normalized', 'Position',position);
        axes_handles(row,col) = ax;
        hold(ax, 'on');
        set(ax, 'Color','w', 'FontName',font_name, 'FontSize',8, ...
            'LineWidth',0.65, 'Box','on', 'Layer','top', ...
            'XColor','k', 'YColor','k', 'TickDir','out', ...
            'YScale','log', 'YDir','normal', 'CLim',climits, ...
            'XGrid','off', 'YGrid','off', 'XMinorTick','off', ...
            'YMinorTick','on');

        draw_nonuniform_cells(ax, p.phi_edges, p.lambda_edges, p.growth, cmap, climits);
        xlim(ax, [-90,90]);
        ylim(ax, [p.lambda(1),p.lambda(end)]);        % small bottom, large top
        set(ax, 'XTick',[-90,-45,0,45,90]);
        yt = [0.1,0.2,0.5,1,2,5,10,30];
        yt = yt(yt >= p.lambda(1) & yt <= p.lambda(end));
        set(ax, 'YTick',yt, ...
            'YTickLabel',arrayfun(@(v) sprintf('%g',v),yt,'UniformOutput',false));

        if show_zero_contour && min(p.growth(:)) < 0 && max(p.growth(:)) > 0
            contour(ax, p.phi, p.lambda, p.growth, [0,0], ...
                'LineColor','w', 'LineWidth',0.42);
        end
        if show_si_axis
            line(ax, [0,0], [p.lambda(1),p.lambda(end)], ...
                'Color','k', 'LineStyle','--', 'LineWidth',0.65);
        end

        title(ax, sprintf('{\\it Ri} = %.2f', p.Ri), ...
            'Interpreter','tex', 'FontName',font_name, ...
            'FontSize',10, 'FontWeight','normal', 'Color','k');
        if row == 2
            xlabel(ax, 'Wave-vector angle \phi (deg)', 'Interpreter','tex', ...
                'FontName',font_name, 'FontSize',8.8, 'Color','k');
        end
        if col == 1
            ylabel(ax, 'Wavelength \lambda/L', 'Interpreter','tex', ...
                'FontName',font_name, 'FontSize',8.8, 'Color','k');
        end

        text_box(fig, [(panel_left+0.13*PW)/W, ...
            (H-panel_top-0.012*PH-15)/H, 18/W, 15/H], ...
            letters(ii), 11, 'bold', 'left', font_name);
    end
end

%% Row headers: centered; no descriptive prefix and no separator bar
text_box(fig, [0.0225,(H-28)/H,3*PW/W,17/H], ...
    'E_{v} = 0', 10.3, 'bold', 'center', font_name);
second_baseline = top + PH + row_gap - 9;
text_box(fig, [0.024,(H-second_baseline-10)/H,3*PW/W,17/H], ...
    'E_{v} = 10^{-3}', 10.3, 'bold', 'center', font_name);

%% Shared colorbar, positioned without moving any panel
cbW = 0.80*72;
cbH = 4.03*72;
cbLeft = 3*PW + 4;
cbTop = top + 34;
cbPosition = [(cbLeft+0.10*cbW)/W, ...
    (H-cbTop-0.95*cbH)/H, 0.22*cbW/W, 0.90*cbH/H];
cbAx = axes('Parent',fig, 'Units','normalized', ...
    'Position',cbPosition, 'Visible','off', 'CLim',climits);

cb = colorbar(cbAx, 'Location','eastoutside');
cb.AxisLocation = 'out';
cb.YAxisLocation = 'right';

set(cb, ...
    'Units','normalized', 'Position',cbPosition, ...
    'Limits',climits, 'Ticks',0:0.2:1.6, 'FontName',font_name, ...
    'FontSize',8, 'Color','k', 'LineWidth',0.65);

cb.Label.String = 'Growth rate Re(\sigma)';
cb.Label.Interpreter = 'tex';
cb.Label.FontName = font_name;
cb.Label.FontSize = 9;
cb.Label.Color = 'k';


drawnow;
exportgraphics(fig,fullfile(out_dir,'Figure1_N80.pdf'), ...
    'ContentType',pdf_content_type,'Resolution',pdf_dpi,'BackgroundColor','white');
exportgraphics(fig,fullfile(out_dir,'Figure1_N80.tiff'), ...
    'Resolution',tiff_dpi,'BackgroundColor','white');
savefig(fig,fullfile(out_dir,'Figure1_N80.fig'));
function edge = cell_edges(x,use_log)
x = x(:);
if use_log, x = log(x); end
edge = [x(1)-(x(2)-x(1))/2; (x(1:end-1)+x(2:end))/2; ...
    x(end)+(x(end)-x(end-1))/2];
if use_log, edge = exp(edge); end
end

function draw_nonuniform_cells(ax,xedge,yedge,g,cmap,climits)
[X,Y] = meshgrid(xedge,yedge);
id = reshape(1:numel(X),size(X));
v1 = id(1:end-1,1:end-1);
v2 = id(1:end-1,2:end);
v3 = id(2:end,2:end);
v4 = id(2:end,1:end-1);
faces = [v1(:),v2(:),v3(:),v4(:)];
rgb = growth_to_rgb(g(:),cmap,climits);
patch('Parent',ax,'Faces',faces,'Vertices',[X(:),Y(:)], ...
    'FaceVertexCData',rgb,'FaceColor','flat','EdgeColor','none', ...
    'FaceLighting','none','Clipping','on');
end

function rgb = growth_to_rgb(values,cmap,climits)
nc = size(cmap,1);
t = (values(:)-climits(1))/(climits(2)-climits(1));
idx = min(nc,max(1,floor(nc*t)+1));
rgb = cmap(idx,:);
rgb(values(:)<0,:) = 1;        % THE ONLY white-cell rule
end

function text_box(fig,position,str,fontsize,weight,align,font_name)
annotation(fig,'textbox',position,'String',str,'Interpreter','tex', ...
    'LineStyle','none','Margin',0,'FitBoxToText','off', ...
    'HorizontalAlignment',align,'VerticalAlignment','top', ...
    'FontName',font_name,'FontSize',fontsize,'FontWeight',weight,'Color','k');
end

