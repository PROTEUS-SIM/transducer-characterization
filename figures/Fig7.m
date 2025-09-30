% Figure formatting for Fig. 6.7.
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
fig2 = openfig("fig/Fig7b.fig",visibility); ax2 = fig2.Children;
fig3 = openfig("fig/Fig7c.fig",visibility); ax3 = fig3.Children;
fig4 = openfig("fig/Fig7d.fig",visibility); ax4 = fig4.Children;
fig5 = openfig("fig/Fig7e.fig",visibility); ax5 = fig5.Children(2);
fig6 = openfig("fig/Fig7f.fig",visibility); ax6 = fig6.Children;
lgd = fig5.Children(1);

fig2.WindowStyle = "normal";
fig3.WindowStyle = "normal";
fig4.WindowStyle = "normal";
fig5.WindowStyle = "normal";
fig6.WindowStyle = "normal";

%--------------------------------------------------------------------------
% Color settings
%--------------------------------------------------------------------------
colors

W = linspace(-1,1,256)';
c = abs(W);
R = (1 - c).*[1 1 1] + c.*red2/255;
B = (1 - c).*[1 1 1] + c.*blue1/255;

% Assign negative weights blue, positive weights red:
C = B;
C(W>0,:) = R(W>0,:);

fig2.Colormap = C;
fig4.Colormap = C;
ax2.XLim = [60 67];
ax3.XLim = [60 67];
ax2.CLim = [-750 750];
ax3.YLim = [-750 750];
ax4.CLim = [-26000 26000];
ax4.XLim = [60 67];
ax5.XLim = [0 10];
ax5.YLim = [-80 30];
ax6.XLim = [0 7];
ax6.YLim = [-3 3];

ax2.Children(1).Color = 'black';
ax3.Children(1).Color = red2/255;
ax3.Children(2).Color = blue2/255;
ax5.Children(1).Color = blue1/255;
ax5.Children(2).Color = red2/255;
ax5.Children(3).Color = blue2/255;
ax6.Children(1).Color = blue1/255;
ax6.Children(2).Color = [1 1 1]*0.8;

ax2.Children(1).LineStyle = ':';
ax2.Children(1).LineWidth = 2;

ax3.Children(1).LineWidth = 0.75;
ax3.Children(2).LineWidth = 0.75;
ax5.Children(1).LineWidth = 0.75;
ax5.Children(2).LineWidth = 0.75;
ax5.Children(3).LineWidth = 0.75;

%--------------------------------------------------------------------------
% Labels
%--------------------------------------------------------------------------
ax2.XLabel.String = 'Time (us)';
ax3.XLabel.String = 'Time (us)';

ax3.YLabel.String = 'Voltage (ADCL)';

ax2.Title.String = '';
ax3.Title.String = '';
ax4.Title.String = '';
ax5.Title.String = '';
ax6.Title.String = '';

ax4.XLabel.String = 'Time (us)';
ax5.YLabel.String = 'Amplitude (dB)';
ax6.XLabel.String = 'Time (us)';
ax6.YLabel.String = 'ADCL/Pa/us';

fontName = 'DejaVu Sans';
ax2.FontSize = 8; ax2.FontName = fontName;
ax3.FontSize = 8; ax3.FontName = fontName;
ax4.FontSize = 8; ax4.FontName = fontName;
ax5.FontSize = 8; ax5.FontName = fontName;
ax6.FontSize = 8; ax6.FontName = fontName;

lgd.Location = 'south';
grid(ax5,'on');

%--------------------------------------------------------------------------
% Size settings
%--------------------------------------------------------------------------
figure_height = 2.1; % inch
figure_width  = 2.1; % inch
fig2.Position = compute_figure_position(figure_width,figure_height);
fig4.Position = compute_figure_position(figure_width,figure_height);

figure_height = 1.29; % inch
figure_width  = 2.1;  % inch
fig3.Position = compute_figure_position(figure_width,figure_height);

figure_height = 2.1; % inch
figure_width  = 2.1; % inch
fig5.Position = compute_figure_position(figure_width,figure_height);
fig6.Position = fig5.Position;

ax6.Position(1) = ax5.Position(1);
ax6.Position(3) = ax5.Position(3);

%--------------------------------------------------------------------------
% Export figures
%--------------------------------------------------------------------------
fprintf("Exporting figures ...\n")
print(fig2,'-dsvg',resolution, saveFolder + filesep + "Fig7b")
print(fig3,'-dsvg',resolution, saveFolder + filesep + "Fig7c")
print(fig4,'-dpng',resolution, saveFolder + filesep + "Fig7d")
print(fig4,'-dsvg',resolution, saveFolder + filesep + "Fig7d")
print(fig5,'-dsvg',resolution, saveFolder + filesep + "Fig7e")
print(fig6,'-dsvg',resolution, saveFolder + filesep + "Fig7f")
fprintf("Done.\n")
