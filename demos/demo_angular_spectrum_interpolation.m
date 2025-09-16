% This script demonstrates that the angular spectrum method can also be
% applied to an inclined sensor plane through the use of interpolation.
%
% Figure not shown in my thesis, but the result is briefly mentioned in
% Section 6.2.2.
%
% NOTE: If sufficient RAM is available, change the value of maxMemory in
% update_sensor_fast.m (part of the PROTEUS toolbox) to at least 3 GB.
% Otherwise, update_sensor_fast will revert to update_sensor, which will
% increase the computation time by two orders of magnitude.
%
% This file is part of the transducer-calibration project, licensed under
% the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

clear; clc; close all

% Add root directory to path and run path setup:
currentDirectory = fileparts(mfilename('fullpath'));
rootDirectory = fileparts(currentDirectory);
addpath(rootDirectory); clear currentDirectory rootDirectory
PATHS = path_setup_calibration();

load('GUI_output_parameters.mat', 'Geometry', ...
    'Medium', 'SimulationParameters', 'Transducer', 'Transmit')

%==========================================================================
% Get grid that encapsulate sensor points
%==========================================================================

% Simulate rotation of the measurement plane:
theta = -10/180*pi; % Rotation in radians
% Rotation matrix for rotation about the y axis:
R = get_rotation_matrix(theta,2);

% Distance of the measurement plane to the transducer:
d = 0.05;

% Grid0: The original sensor grid in the measurement coordinate system
% Grid1: Grid for the source representation
% Grid2: Sensor grid in the transducer coordinate system that encloses the
% measurement points.

expansion_factor = 2; % Grid expansion factor
Geometry.Axis = 3; % Propagation axis
Grid0 = define_grid_demos(SimulationParameters,Geometry,expansion_factor);

% Rotate grid points:
points = transform_grid(Grid0,[0 0 d],R);

% Sensor grid:
SimulationParameters.MaxPrime = Inf;
SimulationParameters.PMLmin = 0;
Geometry = encapsulate_grid(points, Geometry, SimulationParameters);
Grid2 = define_grid(SimulationParameters,Geometry);

% Source grid:
Geometry.Domain.Zmin = 0;
Geometry.Domain.Zmax = 0;

Grid1 = define_grid(SimulationParameters,Geometry);

Transducer.Axis = 3; % Propagation axis
Transducer=update_transducer(Transducer,Medium,Grid1,SimulationParameters);

%==========================================================================
% SOURCE REPRESENTATION
%==========================================================================

f0 = Transmit.CenterFrequency;
Transmit.ContinuousWave = true;
S = define_source_grid(Transducer, Transmit, Medium, Grid1, f0, false);

%==========================================================================
% ANGULAR SPECTRUM METHOD
%==========================================================================

c0     = Medium.SpeedOfSound;
rho    = Medium.Density;

% Compute the wavenumber of propagation:
k  = 2*pi*f0/c0;

Grid1 = get_kspace(Grid1);
Nx = Grid1.Nx; dkx = Grid1.dkx;
Ny = Grid1.Ny; dky = Grid1.dky;

d = [0 Grid2.z];

G = compute_spectral_propagator(Nx,dkx,Ny,dky,k,d,true,'vp');

% Apply forward propagation:
S = angular_spectrum(S,Grid1.x,Grid1.y);
S = S.*G*rho*c0;
S = angular_spectrum_inverse(S,Grid1.x,Grid1.y);

P = S(:,:,2:end);

sensor.mask = zeros(size(P),'logical');

%==========================================================================
% INTERPOLATE THE SENSOR DATA
%==========================================================================

% Get the indices of the nearest grid points:
[points, ~, points_idx, ~] = voxelize_media_points(points, Grid2);

[sensor, sensor_weights] = update_sensor_fast(...
    sensor, points, points_idx, Grid2, false);

P = P(sensor.mask);
P = full(sensor_weights*P);

P = reshape(P,Grid0.Nx,Grid0.Ny);

%==========================================================================
% SHOW RESULTS
%==========================================================================

displayUnit = 1e-3; % [m]

figure
    
% Plot the pressure field at the transducer surface:
[X,Y,Z] = ndgrid(Grid1.x,Grid1.y,0);
X = X/displayUnit;
Y = Y/displayUnit;
Z = Z/displayUnit;

p = surf(Z,X,Y,abs(S(:,:,1)));
p.EdgeColor = 'none';

hold on

X = points(:,1); X = reshape(X,Grid0.Nx,Grid0.Ny)/displayUnit;
Y = points(:,2); Y = reshape(Y,Grid0.Nx,Grid0.Ny)/displayUnit;
Z = points(:,3); Z = reshape(Z,Grid0.Nx,Grid0.Ny)/displayUnit;
p = surf(Z,X,Y,abs(P));
p.EdgeColor = 'none';

hold off

xlabel('z (mm)')
ylabel('x (mm)')
zlabel('y (mm)')

ax = gca;
set(ax,'DataAspectRatio',[1 1 1])


function sensorPoints = transform_grid(Grid,translation,R)
% Transform the grid of sensor points

[X,Y,Z] = ndgrid(Grid.x,Grid.y,Grid.z);
rotationOrigin = [0 0 0];
sensorPoints = [X(:) Y(:) Z(:)];

sensorPoints = sensorPoints - rotationOrigin;
sensorPoints = R*transpose(sensorPoints);
sensorPoints = transpose(sensorPoints) + rotationOrigin;

sensorPoints = sensorPoints + translation;

end

function Geometry = encapsulate_grid(sensorPoints, ...
    Geometry,SimulationParameters)
% Create a grid around the sensor points that encapsulates all sensor
% points

margin = 4*SimulationParameters.GridSize;

D = Geometry.Domain;

D.Xmin = min(sensorPoints(:,1)) - margin;
D.Xmax = max(sensorPoints(:,1)) + margin;
D.Ymin = min(sensorPoints(:,2)) - margin;
D.Ymax = max(sensorPoints(:,2)) + margin;
D.Zmin = min(sensorPoints(:,3)) - margin;
D.Zmax = max(sensorPoints(:,3)) + margin;

Geometry.Domain = D;

end
