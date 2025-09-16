% Figure formatting for Fig. 6.8.
% This file is part of the transducer-characterization project, licensed
% under the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

clear; clc; close all

saveFolder = "exports";
resolution = '-r600';  % dpi

if not(isfolder(saveFolder))
    mkdir(saveFolder)
end

%--------------------------------------------------------------------------
% Load figures
%--------------------------------------------------------------------------
visibility = "Invisible";
fig1 = openfig("fig/Fig8a.fig",visibility); fig1.WindowStyle = "normal";
fig4 = openfig("fig/Fig8d.fig",visibility); fig4.WindowStyle = "normal";

ax1 = fig1.Children(4);
ax2 = fig1.Children(3);
ax3 = fig1.Children(2);
ax4 = fig4.Children(1);
h   = fig1.Children(1);

%--------------------------------------------------------------------------
% Color settings
%--------------------------------------------------------------------------
load('CustomColorMap.mat')
colorMap = CustomColorMap;
colorMap(:,1) = colorMap(:,1)/max(colorMap(:,1));
colorMap(:,2) = colorMap(:,2)/max(colorMap(:,2));

colors % Load custom colors
backgroundcolor = 0.1*blue3/255 + 0.9*[1 1 1];

ax4.Children(1).Color = blue1/255;
ax4.Children(2).Color = red2/255;
ax4.Children(3).Color = blue2/255;

%--------------------------------------------------------------------------
% Labels
%--------------------------------------------------------------------------
legend(ax4,'Old calibration','New calibration','Experiment')

font_size = 8;

fig1.Colormap = colorMap;

ax1.FontName = 'Arial';
ax2.FontName = 'Arial';
ax3.FontName = 'Arial';
ax4.FontName = 'Arial';

ax1.FontSize = font_size;
ax2.FontSize = font_size;
ax3.FontSize = font_size;
ax4.FontSize = font_size;

ax1.YLabel.String = 'y (mm)';
ax2.YLabel.String = 'y (mm)';
ax3.YLabel.String = 'y (mm)';

ax1.XTick = [];
ax2.XTick = [];

ax1.XLabel.String = '';
ax2.XLabel.String = '';
ax3.XLabel.String = 'z (mm)';
ax4.XLabel.String = 'z (mm)';

ax4.XLim = [0 130];

h.Label.FontSize = font_size;

%--------------------------------------------------------------------------
% Size settings
%--------------------------------------------------------------------------
figure_height = 1.7; % inch
figure_width  = 2.9; % inch
fig1.Position = compute_figure_position(figure_width,figure_height);

figure_height = 1.7; % inch
figure_width  = 2.1; % inch
fig4.Position = compute_figure_position(figure_width,figure_height);

%--------------------------------------------------------------------------
% Export figures
%--------------------------------------------------------------------------
fprintf("Exporting figures ...\n")
print(fig1,'-dsvg',resolution, saveFolder + filesep + "Fig8a")
print(fig4,'-dsvg',resolution, saveFolder + filesep + "Fig8d")

h.Visible = 'off';
print(fig1,'-dpng',resolution, saveFolder + filesep + "Fig8a")
fprintf("Done.\n")
