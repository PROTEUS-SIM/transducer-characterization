% Figure formatting for Fig. 6.9.
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
fig1 = openfig("fig/Fig9a.fig",visibility); ax1 = fig1.Children;
fig2 = openfig("fig/Fig9b.fig",visibility); ax2 = fig2.Children;
fig3 = openfig("fig/Fig9c.fig",visibility); ax3 = fig3.Children;
fig4 = openfig("fig/Fig9d.fig",visibility); ax4 = fig4.Children;

fig1.WindowStyle = "normal";
fig2.WindowStyle = "normal";
fig3.WindowStyle = "normal";
fig4.WindowStyle = "normal";

ax3 = flipud(ax3);

%--------------------------------------------------------------------------
% Color settings
%--------------------------------------------------------------------------
load('CustomColorMap.mat')
colorMap = CustomColorMap;
colorMap(:,1) = colorMap(:,1)/max(colorMap(:,1));
colorMap(:,2) = colorMap(:,2)/max(colorMap(:,2));

colors % Load custom colors
backgroundcolor = 0.1*blue3/255 + 0.9*[1 1 1];

W = linspace(-1,1,256)';
c = abs(W);
R = (1 - c).*[1 1 1] + c.*red2/255;
B = (1 - c).*[1 1 1] + c.*blue1/255;

% Assign negative weights blue, positive weights red:
C = B;
C(W>0,:) = R(W>0,:);

ax1.Children.Color = blue2/255;
ax3(1).Children(1).Color = [1 1 1]*0.5;
ax3(2).Children(1).Color = [1 1 1]*0.5;
ax3(3).Children(1).Color = [1 1 1]*0.5;
ax3(1).Children(2).Color = red2/255;
ax3(2).Children(2).Color = red2/255;
ax3(3).Children(2).Color = red2/255;
ax3(1).Children(2).LineWidth = 0.75;
ax3(2).Children(2).LineWidth = 0.75;
ax3(3).Children(2).LineWidth = 0.75;

fig2.Colormap = C;
fig4.Colormap = colorMap;

%--------------------------------------------------------------------------
% Labels
%--------------------------------------------------------------------------
ax1.XLabel.String = 'Time (us)';

ax2.XLabel.String = 'y (mm)';
ax2.YLabel.String = 'Time (us)';

ax3(1).XLabel.String = '';
ax3(2).XLabel.String = '';
ax3(3).XLabel.String = 'y (mm)';
ax3(1).YLabel.String = '';
ax3(3).YLabel.String = '';

ax4.XLabel.String = 'k (um^{-1})';
ax4.YLabel.String = 'w (MHz)';
ax4.Title.String = '' ;

ax1.FontSize = 8; ax1.FontName = 'Arial';
ax2.FontSize = 8; ax2.FontName = 'Arial';
ax3(1).FontSize = 8; ax3(1).FontName = 'Arial';
ax3(2).FontSize = 8; ax3(2).FontName = 'Arial';
ax3(3).FontSize = 8; ax3(3).FontName = 'Arial';
ax4.FontSize = 8; ax4.FontName = 'Arial';

ax1.XLim = [0 12];
ax1.YLim = [-30 30];

ax2.CLim = [-0.03 0.03];
ax2.YLim = [0 12];

ax3(1).YLim = [-7 7];
ax3(2).YLim = [-7 7];
ax3(3).YLim = [-7 7];

ax3(1).XTick = [];
ax3(2).XTick = [];

ax4.XLim = [0 0.06];
ax4.YLim = [0 35];

%--------------------------------------------------------------------------
% Size settings
%--------------------------------------------------------------------------

figure_height = 1.4; % inch
figure_width  = 2.1; % inch
fig1.Position = compute_figure_position(figure_width,figure_height);

figure_height = 1.7; % inch
figure_width  = 2.1; % inch
fig2.Position = compute_figure_position(figure_width,figure_height);

figure_height = 3.1; % inch
figure_width  = 2.5; % inch
fig3.Position = compute_figure_position(figure_width,figure_height);

figure_height = 1.7; % inch
figure_width  = 2.1; % inch
fig4.Position = compute_figure_position(figure_width,figure_height);

%--------------------------------------------------------------------------
% Export figures
%--------------------------------------------------------------------------
fprintf("Exporting figures ...\n")
print(fig1,'-dsvg',resolution, saveFolder + filesep + "Fig9a")
print(fig2,'-dsvg',resolution, saveFolder + filesep + "Fig9b")
print(fig3,'-dsvg',resolution, saveFolder + filesep + "Fig9c")
print(fig4,'-dsvg',resolution, saveFolder + filesep + "Fig9d")
fprintf("Done.\n")
