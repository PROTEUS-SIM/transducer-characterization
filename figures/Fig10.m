% Figure formatting for Fig. 6.10.
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
% Load figures
%--------------------------------------------------------------------------
visibility = "Invisible";
fig1 = openfig("fig/Fig10a.fig",visibility); ax1 = fig1.Children(2);
fig2 = openfig("fig/Fig10b.fig",visibility); ax2 = fig2.Children(2);
fig3 = openfig("fig/Fig10c.fig",visibility); ax3 = fig3.Children(2);

fig1.WindowStyle = "normal";
fig2.WindowStyle = "normal";
fig3.WindowStyle = "normal";

ax1.Title.String = '';
ax2.Title.String = '';
ax3.Title.String = '';

%--------------------------------------------------------------------------
% Color settings
%--------------------------------------------------------------------------
colors

ax1.Children(1).Color = red2/255;
ax2.Children(1).Color = red2/255;
ax3.Children(1).Color = red2/255;
ax1.Children(2).Color = blue2/255;
ax2.Children(2).Color = blue2/255;
ax3.Children(2).Color = blue2/255;

ax1.Children(1).LineWidth = 0.75;
ax2.Children(1).LineWidth = 0.75;
ax3.Children(1).LineWidth = 0.75;
ax1.Children(2).LineWidth = 0.75;
ax2.Children(2).LineWidth = 0.75;
ax3.Children(2).LineWidth = 0.75;

%--------------------------------------------------------------------------
% Labels
%--------------------------------------------------------------------------
ax1.FontSize = 8; ax1.FontName = 'DejaVu Sans';
ax2.FontSize = 8; ax2.FontName = 'DejaVu Sans';
ax3.FontSize = 8; ax3.FontName = 'DejaVu Sans';

ax2.XLabel.String = 'Time (us)';
ax3.XLabel.String = 'Time (us)';

%--------------------------------------------------------------------------
% Size settings
%--------------------------------------------------------------------------
figure_height = 1.5; % inch
figure_width  = 2.3; % inch
fig1.Position = compute_figure_position(figure_width,figure_height);

figure_height = 1.5; % inch
figure_width  = 2.3; % inch
fig2.Position = compute_figure_position(figure_width,figure_height);
fig3.Position = compute_figure_position(figure_width,figure_height);

%--------------------------------------------------------------------------
% Export figures
%--------------------------------------------------------------------------
fprintf("Exporting figures ...\n")
print(fig1,'-dsvg',resolution, saveFolder + filesep + "Fig10a")
print(fig2,'-dsvg',resolution, saveFolder + filesep + "Fig10b")
print(fig3,'-dsvg',resolution, saveFolder + filesep + "Fig10c")
fprintf("Done.\n")
