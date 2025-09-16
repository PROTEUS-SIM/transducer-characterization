% Comparison of a rayleigh integral computation and an angular spectrum
% computation (one-dimensional cross section profile).
%
% This demonstration is not included in my thesis.
%
% This file is part of the transducer-characterization project, licensed
% under the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

clear; clc; close all

% Add root directory to path and run path setup:
currentDirectory = fileparts(mfilename('fullpath'));
rootDirectory = fileparts(currentDirectory);
addpath(rootDirectory); clear currentDirectory rootDirectory
path_setup_characterization;

% Distance between source plane and sensor plane [m]:
d = 0.075;

%==========================================================================
% SOURCE REPRESENTATION 1
%==========================================================================

load('GUI_output_parameters.mat', 'Geometry', ...
    'Medium', 'SimulationParameters', 'Transducer', 'Transmit')

Geometry.Axis = 3; % Propagation axis
Grid1 = define_grid_demos(SimulationParameters,Geometry);

Transducer.Axis = 3; % Propagation axis
Transducer=update_transducer(Transducer,Medium,Grid1,SimulationParameters);

useGPU = true;
Transmit.ContinuousWave = true;
source = define_source_frequency(Transducer, Transmit, Medium, Grid1, ...
    Transmit.CenterFrequency,useGPU);

%==========================================================================
% SOURCE REPRESENTATION 2
%==========================================================================

load('GUI_output_parameters.mat', 'Geometry', ...
    'Medium', 'SimulationParameters', 'Transducer', 'Transmit')

expansion_factor = 4; % Grid expansion factor

Geometry.Axis = 3; % Propagation axis
Grid2 = define_grid_demos(SimulationParameters,Geometry,expansion_factor);

Transducer.Axis = 3; % Propagation axis
Transducer=update_transducer(Transducer,Medium,Grid2,SimulationParameters);

% Field distribution representation of transducer:
f0 = Transmit.CenterFrequency;
Transmit.ContinuousWave = true;
S = define_source_grid(Transducer, Transmit, Medium, Grid2, f0, false);

%==========================================================================
% METHOD 1: RAYLEIGH INTEGRAL
%==========================================================================

[X,Y,Z] = ndgrid(0,Grid1.y,d);

medium.density     = Medium.Density;
medium.sound_speed = Medium.SpeedOfSound;
medium.direction   = 'forward';
medium.useGPU      = useGPU;

sensor.points = [X(:) Y(:) Z(:)];

[source,sensor] = rayleigh_integral_frequency(source, sensor, medium);
p1 = sensor.pressures;

%==========================================================================
% METHOD 2: ANGULAR SPECTRUM
%==========================================================================

rho = Medium.Density;
c0  = Medium.SpeedOfSound;

% Compute the wavenumber of propagation:
k  = 2*pi*f0/c0;

Nx = Grid2.Nx; dkx = 2*pi/(Nx*Grid2.dx);
Ny = Grid2.Ny; dky = 2*pi/(Ny*Grid2.dy);

% Apply the 2D FFT, apply the spectral propagator, and apply the inverse 2D
% FFT:
G = compute_spectral_propagator(Nx,dkx,Ny,dky,k,d,true,'vp');
Sfft = angular_spectrum(S,Grid2.x,Grid2.y);
Sfft = Sfft.*G*rho*c0;
S2 = angular_spectrum_inverse(Sfft,Grid2.x,Grid2.y);

% Extract the line x == 0:
p2 = S2(Grid2.x == 0,:);

%==========================================================================
% SHOW RESULTS
%==========================================================================
plot(Grid1.y,abs(p1))
hold on
plot(Grid2.y,abs(p2))
xlabel('Elevation (m)')
ylabel('Pressure (Pa)')
legend('Rayleigh integral', 'Angular spectrum')

figure
plot(Grid1.y,angle(p1))
hold on
plot(Grid2.y,angle(p2))
xlabel('Elevation (m)')
ylabel('Phase (rad)')
legend('Rayleigh integral', 'Angular spectrum')
