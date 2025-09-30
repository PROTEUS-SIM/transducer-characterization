% Figure formatting for Fig. 5.
% This file is part of the transducer-characterization project, licensed
% under the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

%#ok<*UNRCH>

clear; clc; close all

saveFolder = "exports";
plotType   = 'abs';
saveFormat = 'png';

resolution = '-r600';  % dpi

visibility = "Invisible";

if strcmp(plotType,'abs')
    fig = openfig("fig/Fig5a.fig",visibility);
    saveName = "Fig5a";
elseif strcmp(plotType,'angle')
    fig = openfig("fig/Fig5b.fig",visibility);
    saveName = "Fig5b";
end

fig.WindowStyle = "normal";

ax = fig.Children;

ax.DataAspectRatio = [1 10 10];
ax.XLim = [0 4];
ax.YLim = [-20 20];
ax.ZLim = [-15 15];
ax.Title.String = '';

h1 = ax.Children(4);
h2 = ax.Children(3);
h3 = ax.Children(2);
h4 = ax.Children(1);

if not(isfolder(saveFolder))
    mkdir(saveFolder)
end

% Color settings
load('CustomColorMap.mat')
colorMap = CustomColorMap;
colorMap(:,1) = colorMap(:,1)/max(colorMap(:,1));
colorMap(:,2) = colorMap(:,2)/max(colorMap(:,2));

colors % Load custom colors
backgroundcolor = 0.1*blue3/255 + 0.9*[1 1 1];

fig.Colormap = colorMap;
fig.InvertHardcopy = 'off';
fig.Color = 'w';
ax.Color = backgroundcolor;

cbar = colorbar;

h4.Color = red2/255;

% Size settings
figure_height = 2.5; % inch
figure_width  = 2.5; % inch
fig.Position = compute_figure_position(figure_width,figure_height);

h1.SizeData = 1;
h2.SizeData = 1;

ax.FontSize = 8;
ax.FontName = 'Arial';

grid(ax,'off')
box(ax,'on')

cbar.Location = 'southoutside';
cbar.Position = [0.5 0.1 0.4 0.02];

if strcmp(plotType,'abs')
    clim([0 0.014])
    cbar.Ticks = [0 0.007 0.014];
else
    clim([-1 1])
end

fprintf("Exporting figure ...\n")

% Save figure
if strcmp(saveFormat,'png')
    delete(h3)
    ax.LineWidth = 0.01;
    print(fig,'-dpng',resolution, saveFolder + filesep + saveName + ".png")
elseif strcmp(saveFormat,'pdf')
    box(ax,'off')
    delete(h1)
    delete(h2)
    print(fig,'-dpdf',resolution, saveFolder + filesep + saveName + ".pdf")
elseif strcmp(saveFormat,'svg')
    delete(h1)
    delete(h2)
    print(fig,'-dsvg',resolution, saveFolder + filesep + saveName + ".svg")
end

fprintf("Done.\n")
