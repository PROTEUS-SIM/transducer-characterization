% This script demonstrates how to find the measurement plane orientation
% (Section 6.2.4 of my thesis).
%
% Creates the base plots for Fig. 6.3 of my thesis.
%
% First create the following file with demo_rayleigh_plane.m:
% data/example_sensor_data.mat
%
% This file is part of the transducer-calibration project, licensed under
% the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

%#ok<*UNRCH>

clear; clc; close all

% Add root directory to path and run path setup:
currentDirectory = fileparts(mfilename('fullpath'));
rootDirectory = fileparts(currentDirectory);
resultsDirectory = fullfile(currentDirectory,'data');
addpath(rootDirectory); clear currentDirectory rootDirectory
PATHS = path_setup_calibration();

savename = 'example_sensor_data.mat';

saveFigure = false;

load(fullfile(resultsDirectory,savename),'P','R','Grid','k')

P = angular_spectrum(P,Grid.x,Grid.y);
Grid = get_kspace(Grid);

fprintf('| X (deg) | Y (deg) | Z (deg) |\n')
options.logscale = true;
[thetaX, thetaY, thetaZ] = find_angles_rectangular(P,k,Grid,options);
axis equal
xlim([-6 6])
ylim([-6 6])
cbar = colorbar;

fprintf('| %7.2f | %7.2f | %7.2f |\n',...
    thetaX*180/pi,thetaY*180/pi,thetaZ*180/pi)

if not(saveFigure)
    return
end

saveFolder = string(fullfile(PATHS.Figures,'fig'));

if not(isfolder(saveFolder))
    mkdir(saveFolder)
end

savefig(gcf, saveFolder + filesep + "Fig3.fig")
