% Figure formatting for Fig. 6.1.
% This file is part of the transducer-calibration project, licensed under
% the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

%#ok<*UNRCH>

clear; clc; close all

saveFolder = "exports";
saveFormat = 'png';
resolution  = '-r600';  % dpi

if not(isfolder(saveFolder))
    mkdir(saveFolder)
end

%--------------------------------------------------------------------------
% Load figures
%--------------------------------------------------------------------------
visibility = "Invisible";
fig1 = openfig("fig/Fig1a.fig",visibility); ax1 = fig1.Children;
fig2 = openfig("fig/Fig1b.fig",visibility); ax2 = fig2.Children;
fig3 = openfig("fig/Fig1c.fig",visibility); ax3 = fig3.Children;
fig4 = openfig("fig/Fig1d.fig",visibility); ax4 = fig4.Children(2);
fig5 = openfig("fig/Fig1e.fig",visibility); ax5 = fig5.Children(2);
fig6 = openfig("fig/Fig1f.fig",visibility); ax6 = fig6.Children(2);

cbar4 = fig4.Children(1);
cbar5 = fig5.Children(1);
cbar6 = fig6.Children(1);

figList = {fig1,fig4,fig2,fig5,fig3,fig6};
axList  = {ax1,ax4,ax2,ax5,ax3,ax6};

saveNames = ["Fig1a","Fig1d","Fig1b","Fig1e","Fig1c","Fig1f"];

dataTypes = ["space","kspace","space","kspace","space","kspace"];

load('CustomColorMap.mat')
colorMap = CustomColorMap;
colorMap(:,1) = colorMap(:,1)/max(colorMap(:,1));
colorMap(:,2) = colorMap(:,2)/max(colorMap(:,2));

XLim = [-30 30];
YLim = [-15 15];

figure_height = 3; % inch
figure_width  = 3; % inch

image_height = 1.3;
image_width  = 1.3;

cbar4.Visible = 'off';
cbar5.Visible = 'off';

if strcmp(saveFormat,'png')
    cbar6.Visible = 'off';
end

for n = 1:length(figList)
    fig = figList{n};
    ax  = axList{n};
    fig.WindowStyle = "normal";
    fig.Position = compute_figure_position(figure_width,figure_height);
    ax.Units = 'pixels';
    ax.Position = compute_axes_position(fig,image_width,image_height);
    
    if strcmp(dataTypes(n),"space")
        ax.XLim = XLim;
        ax.YLim = YLim;
        ax.XLabel.String = 'x (mm)';
        ax.YLabel.String = 'y (mm)';
    else
        ax.XLabel.String = 'k_x (mm^{-1})';
        ax.YLabel.String = 'k_y (mm^{-1})';
    end
    
    fig.Colormap = colorMap;
    ax.FontSize  = 8;
    ax.FontName  = 'arial';
    ax.Title.Visible = 'off';
    
    %======================================================================
    % Resize figure panels
    %======================================================================
    
    if ~(n==1 || n==2)
        ax.YAxis.Label.Visible = 'off';
        ax.YAxis.TickLabels = [];
    end

    margin = 5;
       
    if strcmp(dataTypes(n),"space")
        aspectRatio = (YLim(2)-YLim(1))/(XLim(2)-XLim(1));
    else
        aspectRatio = 1;
    end
    
    ax.Position(4) = aspectRatio*ax.Position(4) + margin*2;
    
    if n==1 || n==2
        ax.Position(1) = 42;
    else
        ax.Position(1) = margin;
    end
    
    if strcmp(dataTypes(n),"space")
        ax.Position(2) = 35 - margin;
    else
        ax.Position(2) = 42 - margin;
    end
    
    fig.Position(3) = ax.Position(1) + ax.Position(3) + margin;
    fig.Position(4) = ax.Position(2) + ax.Position(4);
    
    if n == 5 || n == 6
        fig.Position(3) = fig.Position(3) + 42;
    end
    
    %======================================================================
    % Export images
    %======================================================================

    fprintf("Exporting figure %.f out of %.f ...\n", n, length(figList))
    
    if strcmp(saveFormat,'png')
        % Do not show any labels for png image
        ax.XAxis.Visible = 'off';
        ax.YAxis.Visible = 'off';
        ax.ZAxis.Visible = 'off';
        ax.Title.Visible = 'off';
        print(fig,'-dpng',resolution, saveFolder + filesep + saveNames(n))
    elseif strcmp(saveFormat,'pdf')
        print(fig,'-dpdf',resolution, saveFolder + filesep + saveNames(n))
    elseif strcmp(saveFormat,'svg')
        print(fig,'-dsvg',resolution, saveFolder + filesep + saveNames(n))
    end
end

fprintf("Done.\n")
