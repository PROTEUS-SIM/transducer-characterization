% Figure formatting for Fig. 6.2.
% This file is part of the transducer-characterization project, licensed
% under the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

clear; clc; close all

saveFolder = "exports";
resolution  = '-r600';  % dpi

if not(isfolder(saveFolder))
    mkdir(saveFolder)
end

%--------------------------------------------------------------------------
% Load figure
%--------------------------------------------------------------------------
visibility = "Invisible";
fig = openfig("fig/Fig2.fig",visibility);
fig.WindowStyle = "normal";
ax = fig.Children(2);
lgd = fig.Children(1);

%--------------------------------------------------------------------------
% Color settings
%--------------------------------------------------------------------------
colors
ax.Children(4).LineStyle = '--';
ax.Children(5).LineStyle = '--';
ax.Children(6).LineStyle = '--';

ax.Children(1).Color = red2/255;
ax.Children(2).Color = blue2/255;
ax.Children(3).Color = blue1/255;

ax.Children(4).Color = red2/255;
ax.Children(5).Color = blue2/255;
ax.Children(6).Color = blue1/255;

ax.Children(1).LineWidth = 0.75;
ax.Children(2).LineWidth = 0.75;
ax.Children(3).LineWidth = 0.75;
ax.Children(4).LineWidth = 0.75;
ax.Children(5).LineWidth = 0.75;
ax.Children(6).LineWidth = 0.75;

%--------------------------------------------------------------------------
% Labels
%--------------------------------------------------------------------------
ax.FontSize = 8;
ax.FontName = 'Arial';

ax.XLabel.String = 'Time (us)';

lgd.Location = 'southeast';

%--------------------------------------------------------------------------
% Size settings
%--------------------------------------------------------------------------
figure_height = 2.6; % inch
figure_width  = 5; % inch
fig.Position = compute_figure_position(figure_width,figure_height);

%--------------------------------------------------------------------------
% Export figures
%--------------------------------------------------------------------------
fprintf("Exporting figure ...\n")

print(fig,'-dsvg',resolution, saveFolder + filesep + "Fig2")

lgd.Visible = 'off';
ax.XLim = [5 9];
ax.YLim = ax.YLim/15*4;

print(fig,'-dsvg',resolution, saveFolder + filesep + "Fig2inset")

fprintf("Done.\n")
