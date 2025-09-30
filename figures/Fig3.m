% Figure formatting for Fig. 3.
% This file is part of the transducer-characterization project, licensed
% under the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

%#ok<*UNRCH>

clear; clc; close all

saveFolder = "exports";
saveFormat = 'png';
saveName  = "Fig3";

resolution  = '-r600';  % dpi

if not(isfolder(saveFolder))
    mkdir(saveFolder)
end

visibility = "Invisible";
fig = openfig("fig/Fig3.fig",visibility); ax = fig.Children(2);
fig.WindowStyle = "normal";

load('CustomColorMap.mat')
colorMap = CustomColorMap;
colorMap(:,1) = colorMap(:,1)/max(colorMap(:,1));
colorMap(:,2) = colorMap(:,2)/max(colorMap(:,2));

fig.Colormap = colorMap;

figure_height = 3.5; % inch
figure_width  = 3.5; % inch
fig.Position = compute_figure_position(figure_width,figure_height);
ax.FontSize = 8;
ax.FontName = 'Arial';

fprintf("Exporting figure ...\n")

if strcmp(saveFormat,'png')
    % Do not show integration points in png image
    delete(ax.Children(1))

    print(fig,'-dpng',resolution, saveFolder + filesep + saveName)
elseif strcmp(saveFormat,'pdf')
    print(fig,'-dpdf',resolution, saveFolder + filesep + saveName)
elseif strcmp(saveFormat,'svg')
    print(fig,'-dsvg',resolution, saveFolder + filesep + saveName)
end

fprintf("Done.\n")
