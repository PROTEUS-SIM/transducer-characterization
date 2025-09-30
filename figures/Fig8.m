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

ax4.Children(1).Color = red2/255;
ax4.Children(2).Color = blue2/255;
ax4.Children(3).Color = blue1/255;

%--------------------------------------------------------------------------
% Labels
%--------------------------------------------------------------------------
legend(ax4,'a: PROTEUS I exp.','b: PROTEUS I sim.','c: Current study sim.')

font_size = 8;

fig1.Colormap = colorMap;

ax1.FontName = 'DejaVu Sans';
ax2.FontName = 'DejaVu Sans';
ax3.FontName = 'DejaVu Sans';
ax4.FontName = 'DejaVu Sans';

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

ax4.XLim = ax3.XLim;

h.Label.FontSize = font_size;

%--------------------------------------------------------------------------
% Size settings
%--------------------------------------------------------------------------
figure_height = 2.0; % inch
figure_width  = 3.4; % inch
fig1.Position = compute_figure_position(figure_width,figure_height);

figure_height = 1.7; % inch
figure_width  = 3.4; % inch
fig4.Position = compute_figure_position(figure_width,figure_height);

ax4.Position(1) = ax3.Position(1);
ax4.Position(3) = ax3.Position(3);

%--------------------------------------------------------------------------
% Export figures
%--------------------------------------------------------------------------
fprintf("Exporting figures ...\n")
print(fig1,'-dsvg',resolution, saveFolder + filesep + "Fig8a")
print(fig4,'-dsvg',resolution, saveFolder + filesep + "Fig8d")

h.Visible = 'off';
print(fig1,'-dpng',resolution, saveFolder + filesep + "Fig8a")
fprintf("Done.\n")
