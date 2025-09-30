% Simple demonstration of the Rayleigh integral (one-dimensional
% cross-section)
%
% Not included in the arXiv preprint
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

useGPU = true;

%==========================================================================
% SOURCE REPRESENTATION
%==========================================================================

load('GUI_output_parameters.mat', 'Geometry', ...
    'Medium', 'SimulationParameters', 'Transducer', 'Transmit')

f0 = Transmit.CenterFrequency*[0.5 1 1.5 2]; Nf = length(f0);
d  = 0.075;

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

[X,Y,Z] = ndgrid(0,Grid.y,d);

source_type = 'monopole';

% The rigid baffle case can be considered as a sum of monopole sources and
% the soft baffle case as a sum of dipole sources (because of the cosine
% term attached to the integration element):
switch source_type
    case 'monopole'
        baffle = 'rigid';
    case 'dipole'
        baffle = 'soft';
end

source.baffle = baffle;

medium.density     = Medium.Density;
medium.sound_speed = Medium.SpeedOfSound;
medium.direction   = 'forward';
medium.useGPU      = useGPU;

sensor.points = [X(:) Y(:) Z(:)];

[source,sensor] = rayleigh_integral_frequency(source, sensor, medium);
p1 = sensor.pressures;

%==========================================================================
% SHOW RESULTS
%==========================================================================
plot(Grid.y,abs(p1(:,2)))
xlabel('Elevation (m)')
ylabel('Pressure (Pa)')
