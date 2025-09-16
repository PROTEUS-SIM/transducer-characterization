% Figure formatting for Fig. 6.6.
% This file is part of the transducer-calibration project, licensed under
% the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

clear; clc; close all

saveFolder = "exports";
resolution  = '-r600';  % dpi

if not(isfolder(saveFolder))
    mkdir(saveFolder)
end

%--------------------------------------------------------------------------
% Load figures
%--------------------------------------------------------------------------
visibility = "Invisible";
fig1 = openfig("fig/Fig6a1.fig",visibility); ax1 = fig1.Children;
fig2 = openfig("fig/Fig6a2.fig",visibility); ax2 = fig2.Children;
fig3 = openfig("fig/Fig6a.fig" ,visibility); ax3 = fig3.Children;
fig4 = openfig("fig/Fig6b.fig" ,visibility); ax4 = fig4.Children;
fig5 = openfig("fig/Fig6b1.fig",visibility); ax5 = fig5.Children;
fig6 = openfig("fig/Fig6c.fig" ,visibility); ax6 = fig6.Children;
fig7 = openfig("fig/Fig6d.fig" ,visibility); ax7 = fig7.Children;

fig1.WindowStyle = "normal";
fig2.WindowStyle = "normal";
fig3.WindowStyle = "normal";
fig4.WindowStyle = "normal";
fig5.WindowStyle = "normal";
fig6.WindowStyle = "normal";
fig7.WindowStyle = "normal";

%--------------------------------------------------------------------------
% Color settings
%--------------------------------------------------------------------------
load('CustomColorMap.mat')
colorMap = CustomColorMap;
colorMap(:,1) = colorMap(:,1)/max(colorMap(:,1));
colorMap(:,2) = colorMap(:,2)/max(colorMap(:,2));

colors % Load custom colors
backgroundcolor = 0.1*blue3/255 + 0.9*[1 1 1];

fig3.Colormap = colorMap;
fig4.Colormap = colorMap;

ax3.Children(1).LineWidth = 1;
ax3.Children(1).Color = red2/255;
ax3.Children(1).LineStyle = '-';

ax1.Children(1).Color = red2/255;
ax1.Children(2).Color = blue2/255;

ax2.Children(1).Color = red2/255;
ax2.Children(2).Color = blue2/255;

ax5.Children(1).Color = red2/255;
ax5.Children(1).LineWidth = 1;

ax5.Children(2).MarkerSize = 0.5;
ax5.Children(2).Color = blue2/255;

ax6.Children(1).Color = red2/255;
ax6.Children(2).Color = blue2/255;

ax6.Children(1).LineWidth = 0.75;
ax6.Children(2).LineWidth = 0.75;

ax7.Children(1).Color = blue2/255;
ax7.Children(1).LineWidth = 0.75;

%--------------------------------------------------------------------------
% Labels
%--------------------------------------------------------------------------
ax1.Title.String = '';
ax1.XLabel.String = '';
ax1.YLabel.String = '';

ax2.Title.String = '';
ax2.XLabel.String = '';
ax2.YLabel.String = '';

ax3.Title.String = '';
ax3.XLabel.String = 'x (mm)';
ax3.YLabel.String = 'y (mm)';

ax4.Title.String = '';
ax4.XLabel.String = 'x (mm)';
ax4.YLabel.String = 'y (mm)';

ax6.Title.String = '';
ax6.XLabel.String = 'Time (us)';

ax7.Title.String = '';
ax7.XLabel.String = 'Time (us)';
ax7.YLabel.String = 'kPa/V/us';

ax3.FontSize = 8; ax3.FontName = 'Arial';
ax4.FontSize = 8; ax4.FontName = 'Arial';
ax5.FontSize = 8; ax5.FontName = 'Arial';
ax6.FontSize = 8; ax6.FontName = 'Arial';
ax7.FontSize = 8; ax7.FontName = 'Arial';

ax1.XTick = []; ax1.YTick = [];
ax2.XTick = []; ax2.YTick = [];

ax5.Title.String = '';
ax5.XLabel.String = '';
ax5.YLabel.String = '';
ax5.XTick = [];
ax5.YLim = [-0.2 -0.05];

ax6.XLim = [0 3];
ax7.XLim = [0 3];
ax6.YLim = [-40 40];
ax7.YLim = [-150 150];

%--------------------------------------------------------------------------
% Size settings
%--------------------------------------------------------------------------
figure_height = 2.5; % inch
figure_width  = 2.5; % inch
fig3.Position = compute_figure_position(figure_width,figure_height);
ax3.DataAspectRatio = [1 1 1];

fig1.Position = fig3.Position;
fig2.Position = fig3.Position;

aspectRatio = (ax3.YLim(2)-ax3.YLim(1))/(ax3.XLim(2)-ax3.XLim(1));

ax1.YLim = [0 1.5];
ax2.YLim = [0 1.5];
ax1.Position(3) = ax3.Position(3);
ax2.Position(3) = ax3.Position(3)*aspectRatio;

ax1.Position(4) = 0.25;
ax2.Position(4) = 0.25;

fig4.Position = fig3.Position;
ax4.XLim = ax3.XLim;
ax4.YLim = ax3.YLim;
ax4.DataAspectRatio = ax3.DataAspectRatio;

% Rotate plot
view(ax5,[90 -90])

fig5.Position = fig3.Position;
ax5.XLim = ax3.YLim;
ax5.Position(4) = ax2.Position(3);
ax5.Position(3) = 0.25;

figure_height = 1.29; % inch
figure_width  = 1.90; % inch
fig6.Position = compute_figure_position(figure_width,figure_height);
fig7.Position = fig6.Position;
ax6.Position = ax7.Position;

%--------------------------------------------------------------------------
% Export figures
%--------------------------------------------------------------------------
ax5.Children(1).Visible = 'off';
ax5.Children(2).Visible = 'on';

fprintf("Exporting figures to PNG ...\n")

print(fig1,'-dpng',resolution, saveFolder + filesep + "Fig6a1")
print(fig2,'-dpng',resolution, saveFolder + filesep + "Fig6a2")
print(fig3,'-dpng',resolution, saveFolder + filesep + "Fig6a")
print(fig4,'-dpng',resolution, saveFolder + filesep + "Fig6b")
print(fig5,'-dpng',resolution, saveFolder + filesep + "Fig6b1")

ax5.Children(1).Visible = 'on';
ax5.Children(2).Visible = 'off';

fprintf("Exporting figures to SVG ...\n")

print(fig1,'-dsvg',resolution, saveFolder + filesep + "Fig6a1")
print(fig2,'-dsvg',resolution, saveFolder + filesep + "Fig6a2")
print(fig3,'-dsvg',resolution, saveFolder + filesep + "Fig6a")
print(fig4,'-dsvg',resolution, saveFolder + filesep + "Fig6b")
print(fig5,'-dsvg',resolution, saveFolder + filesep + "Fig6b1")
print(fig6,'-dsvg',resolution, saveFolder + filesep + "Fig6c")
print(fig7,'-dsvg',resolution, saveFolder + filesep + "Fig6d")

fprintf("Done.\n")
