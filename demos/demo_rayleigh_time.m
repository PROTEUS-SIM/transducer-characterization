% This script compares a pressure field computed with the time-domain
% Rayleigh integral to a pressure field computed with the frequency-domain
% Rayleigh integral.
%
% The result is not part of the arXiv preprint.
%
% This file is part of the transducer-characterization project, licensed
% under the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

clear; clc; close all

useGPU = true;

% Add root directory to path and run path setup:
currentDirectory = fileparts(mfilename('fullpath'));
rootDirectory = fileparts(currentDirectory);
addpath(rootDirectory); clear currentDirectory rootDirectory
path_setup_characterization;

load('GUI_output_parameters.mat', 'Geometry', ...
    'Medium', 'SimulationParameters', 'Transducer', 'Transmit')

sensorpoint = [-1e-3 0 5e-3];

Geometry.Axis = 3; % Propagation axis
Grid = define_grid_demos(SimulationParameters,Geometry);

% Divide the transducer surfaces into integration elements:
Transducer.Axis = 3; % Propagation axis
Transducer=update_transducer(Transducer,Medium,Grid,SimulationParameters);

%==========================================================================
% COMPUTE THE RAYLEIGH INTEGRAL IN TIME DOMAIN
%==========================================================================

medium.density     = Medium.Density;
medium.sound_speed = Medium.SpeedOfSound;

source = define_source_time(Transducer, Transmit, Medium, Grid, 'rigid');

sensor.points = sensorpoint;

p1 = rayleigh_integral_time(source, sensor, medium);

% Corresponding time array
N = size(p1,2);
t1 = (0:(N-1))/Transmit.SamplingRate;


%==========================================================================
% COMPUTE THE RAYLEIGH INTEGRAL IN FREQUENCY DOMAIN
%==========================================================================

% Downsample signal (Nyquist sampled is sufficient)
downsample = 25;
Transmit.PressureSignal = Transmit.PressureSignal(1:downsample:end);
Transmit.VoltageSignal  = Transmit.VoltageSignal( 1:downsample:end);
Transmit.SamplingRate   = Transmit.SamplingRate/downsample;

% Make the signal longer to accommodate maximum travel times
source.points = reshape(Transducer.integration_points,[],3);
sensor.points = sensorpoint;

[rmin,rmax] = find_min_max_distance(source,sensor);
M = ceil((rmax-rmin)/Medium.SpeedOfSound*Transmit.SamplingRate);
Transmit.PressureSignal = [Transmit.PressureSignal, zeros(1,M)];
sensor.offset = rmin;

% Define medium
medium.direction   = 'forward';
medium.useGPU      = useGPU;

% Define source
source = define_source_frequency(...
    Transducer,Transmit,Medium,Grid,[],useGPU);

% Get batch size for computing the Rayleigh integral:
gpuDevice(1)
availableMemory = gpuDevice().AvailableMemory;
memoryEstimate = get_memory_estimate(source,sensor)*3;
Nbatch = ceil(memoryEstimate/availableMemory);
BS = floor(size(sensor.points,1)/Nbatch);

% Compute the Rayleigh integral
[~,sensor] = rayleigh_integral_batched(source, sensor, medium, BS);

% Convert result to time domain
N = length(Transmit.PressureSignal);
p2 = sensor.pressures;
p2 = ifft(conj(p2),N,2,'symmetric');

% Shifted time array
t2 = (0:(N-1))/Transmit.SamplingRate + rmin/Medium.SpeedOfSound;

%==========================================================================
% SHOW RESULTS
%==========================================================================

% Compute the values of p2 at times t1
p2 = sinc_interpolation(t2,p2,t1);

plot(t1,p1)
hold on
plot(t1,p2)
