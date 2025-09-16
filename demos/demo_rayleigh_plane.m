% This script creates the simulated sensor data that is used in
% demo_find_angles.m, which, in turn, is used to generate Fig. 6.3 of my
% thesis.
%
% Output:
% data/example_sensor_data.mat
%
% This file is part of the transducer-calibration project, licensed under
% the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

clear; clc; close all

% Add root directory to path and run path setup:
currentDirectory = fileparts(mfilename('fullpath'));
rootDirectory = fileparts(currentDirectory);
resultsDirectory = fullfile(currentDirectory,'data');
addpath(rootDirectory); clear currentDirectory rootDirectory
path_setup_characterization;

savename = 'example_sensor_data.mat';

saveResults = false;

useGPU = true;

%==========================================================================
% SOURCE REPRESENTATION
%==========================================================================

load('GUI_output_parameters.mat', 'Geometry', ...
    'Medium', 'SimulationParameters', 'Transducer', 'Transmit')

f0 = Transmit.CenterFrequency;
ifreq = 1;
d  = 0.05;

Geometry.Axis = 3; % Propagation axis
Grid = define_grid_demos(SimulationParameters,Geometry);

Transducer.Axis = 3; % Propagation axis
Transducer=update_transducer(Transducer,Medium,Grid,SimulationParameters);

Transmit.ContinuousWave = true;
source = define_source_frequency(...
    Transducer,Transmit,Medium,Grid,f0,useGPU);

%==========================================================================
% RAYLEIGH INTEGRAL
%==========================================================================

% Simulate rotation of the measurement plane:
thetaX = 6/180*pi; % Rotation in radians
thetaY = 9/180*pi;
thetaZ = 3/180*pi;
% Get rotation matrices:
Rx = get_rotation_matrix(thetaX,1);
Ry = get_rotation_matrix(thetaY,2);
Rz = get_rotation_matrix(thetaZ,3);
% Rotate about x, y, and z axis respectively:
R = Rz*Ry*Rx;

[X,Y,Z] = ndgrid(Grid.x,Grid.y,d);

medium.density     = Medium.Density;
medium.sound_speed = Medium.SpeedOfSound;
medium.direction   = 'forward';
medium.useGPU      = useGPU;

sensor.points = [X(:) Y(:) Z(:)];

% Apply rotation to the sensor points:
sensor.points(:,3) = sensor.points(:,3) - d;
sensor.points = sensor.points*transpose(R);
sensor.points(:,3) = sensor.points(:,3) + d;

% Get batch size for computing the Rayleigh integral:
gpuDevice(1)
availableMemory = gpuDevice().AvailableMemory;
memoryEstimate = get_memory_estimate(source,sensor)*3;
Nbatch = ceil(memoryEstimate/availableMemory);
BS = floor(size(sensor.points,1)/Nbatch);

[source,sensor] = rayleigh_integral_batched(source, sensor, medium, BS);

P = reshape(sensor.pressures,Grid.Nx,Grid.Ny);

%==========================================================================
% SHOW RESULTS
%==========================================================================
displayUnit = 1e-3; % [m]

figure
imagesc(Grid.x/displayUnit,Grid.y/displayUnit,abs(P)')
set(gca, 'YDir', 'normal');
xlabel('y (mm)')
ylabel('x (mm)')

% 3D plot
figure

X = sensor.points(:,1); X = reshape(X,Grid.Nx,Grid.Ny)/displayUnit;
Y = sensor.points(:,2); Y = reshape(Y,Grid.Nx,Grid.Ny)/displayUnit;
Z = sensor.points(:,3); Z = reshape(Z,Grid.Nx,Grid.Ny)/displayUnit;
p = surf(Z,X,Y,abs(P));
p.EdgeColor = 'none';

% Plot the transducer:
hold on
points = source.points/displayUnit;
plot3(points(:,3),points(:,1),points(:,2),'.','MarkerSize',0.1)

xlabel('z (mm)')
ylabel('x (mm)')
zlabel('y (mm)')

ax = gca;
set(ax,'DataAspectRatio',[1 1 1])

%==========================================================================
% SAVE DATA
%==========================================================================

c0 = Medium.SpeedOfSound;
k = 2*pi*f0/c0; % wavenumber

if saveResults == true
    if not(isfolder(resultsDirectory))
        mkdir(resultsDirectory)
    end
    save(fullfile(resultsDirectory,savename),'P','R','Grid','k')
end
