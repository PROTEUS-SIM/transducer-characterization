function rayleigh_figure(sensor,source,...
    Grid,medium,SensorPlane,PATHS,options)
% Creates the base plots for Fig. 6.5 of my thesis. Set options.saveFigure
% = true and run twice: with options.plotType = 'abs' and options.plotType
% = 'angle'.
%
% This file is part of the transducer-calibration project, licensed under
% the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

plotType   = options.plotType;
saveFigure = options.saveFigure;
f0         = options.plotFrequency;

x0 = SensorPlane.x0;
y0 = SensorPlane.y0;
z0 = SensorPlane.z0;

Fs = 1/Grid.dt;

%==========================================================================
% SHOW RESULTS
%==========================================================================

fig = figure;
ax = axes;

[~,Iplot] = min(abs(source.frequencies-f0));


if strcmp(plotType,'abs') && strcmp(source.type,'pressure')
    P1 = abs(source.pressures(:,Iplot))/Fs;
elseif strcmp(plotType,'angle') && strcmp(source.type,'pressure')
    P1 = angle(source.pressures(:,Iplot))/pi;
elseif strcmp(plotType,'abs')
    P1 = abs(source.velocities(:,Iplot))/Fs...
        *medium.density*medium.sound_speed;
elseif strcmp(plotType,'angle')
    P1 = angle(source.velocities(:,Iplot))/pi;
end

if strcmp(plotType,'abs')
    P2 = abs(sensor.pressures(:,Iplot))/Fs;
elseif strcmp(plotType,'angle')
    P2 = angle(sensor.pressures(:,Iplot))/pi;
end

displayUnit = 1e-3; % m
x1 = source.points(:,1)/displayUnit; x2 = sensor.points(:,1)/displayUnit;
y1 = source.points(:,2)/displayUnit; y2 = sensor.points(:,2)/displayUnit;
z1 = source.points(:,3)/displayUnit; z2 = sensor.points(:,3)/displayUnit;

scatter3(z1,x1,y1,[],P1,'.');
hold on
scatter3(z2,x2,y2,[],P2,'.');

if strcmp(medium.direction,'backward')
    % Show scan plan centre
    plot3(z0/displayUnit,x0,y0,'r.','MarkerSize',8);
    
    % Show scan plane before rotation
    rmin = min(source.points);
    rmax = max(source.points);
    x1 = [rmin(1) rmax(1) rmax(1) rmin(1) rmin(1)]/displayUnit;
    y1 = [rmin(2) rmin(2) rmax(2) rmax(2) rmin(2)]/displayUnit;
    z1 = [z0 z0 z0 z0 z0]/displayUnit;
    plot3(z1,x1,y1,'--k');
end

ax.DataAspectRatio = [1 1 1];

if strcmp(plotType,'abs')
    ax.Title.String = 'Source and sensor planes (amplitude)';
else
    ax.Title.String = 'Source and sensor planes (phase)';
end

ax.XLabel.String = 'z (mm)';
ax.YLabel.String = 'x (mm)';
ax.ZLabel.String = 'y (mm)';

drawnow
if not(saveFigure) || strcmp(medium.direction,'forward')
    return
end

saveFolder = string(fullfile(PATHS.Figures,'fig'));

if not(isfolder(saveFolder))
    mkdir(saveFolder)
end

if strcmp(plotType,'abs')
    savefig(fig, saveFolder + filesep + "Fig5a.fig")
elseif strcmp(plotType,'angle')
    savefig(fig, saveFolder + filesep + "Fig5b.fig")
end

end