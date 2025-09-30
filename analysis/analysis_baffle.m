% This script tries to answer the question whether a monopole (pressure) or
% dipole (velocity) representation gives the best representation of the
% experimentally observed data.
%
% The result is inconclusive, as stated in Section IV C of the arXiv
% preprint:
% "However, the edge waves in the hydrophone data were too obscured by
% spurious waves and noise to provide a definitive answer to this
% question."
%
% Run this script three times: with simtypes 'pressure', 'velocity', and
% 'data'.
%
% The resulting figure is not included in the arXiv preprint.
%
% Dependencies:
% data/TransmitData.mat
% data/sensitivity.mat
% data/ScanData.mat
% results/transformation_scan_plane.mat
% results/Transducer
% results/IR_transmit.mat
%
% This file is part of the transducer-characterization project, licensed
% under the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

%#ok<*UNRCH>

clear
clc

simtype = 'pressure';

% Set up the sensor
sensor.points = [10 0 5]*1e-3;

% Add root directory to path and run path setup:
currentDirectory = fileparts(mfilename('fullpath'));
rootDirectory = fileparts(currentDirectory);
addpath(rootDirectory); clear currentDirectory rootDirectory
PATHS = path_setup_characterization();

load(fullfile(PATHS.Results,'IR_transmit.mat'),'IR','Fs')
Tmax = 2e-6;
N = round(Fs*Tmax);
IR = IR(1:N);
t1 = (0:(N-1))/Fs;

load(fullfile(PATHS.Data,'ScanData.mat'),'medium','Grid')

% Make grid spacing identical in all directions:
Grid.dz = Grid.dx;

if strcmp(simtype,'data')

    %======================================================================
    % LOAD HYDROPHONE DATA
    %======================================================================

    load(fullfile(PATHS.Data,'ScanData.mat'),'voltage_mV')

    Fs = 1/Grid.dt;
    f = (0:(Grid.Nt-1))/Grid.Nt*Fs;

    P = voltage_to_pressure(voltage_mV,Fs,3,...
        fullfile(PATHS.Data,'sensitivity.mat'));

    % Outside the circle k^2 = kx^2 + ky^2 in k-space, all nonzero
    % components can be considered to be noise:
    k = 2*pi*f/medium.sound_speed; % wavenumber of propagation
    P = conj(fft(P,[],3));
    P = angular_spectrum(P,Grid.x,Grid.y);
    P = filter_evanescent_waves(P,Grid,k);
    P = angular_spectrum_inverse(P,Grid.x,Grid.y);
    P = ifft(conj(P),[],3,'symmetric');

    [X,Y,Z] = ndgrid(Grid.x,Grid.y,Grid.z);

    load(fullfile(PATHS.Results,'transformation_scan_plane.mat'),'ScanPlane')

    x0 = ScanPlane.x0;
    y0 = ScanPlane.y0;
    z0 = ScanPlane.z0;

    thetaX = ScanPlane.thetaX/180*pi;
    thetaY = ScanPlane.thetaY/180*pi;
    thetaZ = ScanPlane.thetaZ/180*pi;

    % Get rotation matrices:
    Rx = get_rotation_matrix(thetaX,1);
    Ry = get_rotation_matrix(thetaY,2);
    Rz = get_rotation_matrix(thetaZ,3);
    % Rotate about x, y, and z axis respectively:
    R = Rz*Ry*Rx;

else
    
    %======================================================================
    % COMPUTE PRESSURE PULSE
    %======================================================================
    
    % Resample the impulse response
    Fs = single(1/Grid.dt);
    t = t1(1):1/Fs:t1(end);
    IR = sinc_interpolation(t1,IR,t);
    
    % TRANSMIT DATA
    % Load transmit data
    load(fullfile(PATHS.Data,'TransmitData.mat'),'Transmit');
    
    % Integrated value of the drive pulse [V*s]:
    A = sum(Transmit.VoltageSignal)/Transmit.SamplingRate;

    P  = IR*A;
    
    Medium.SpeedOfSound = medium.sound_speed;
    Medium.Density      = medium.density;

    load(fullfile(PATHS.Results,'Transducer.mat'),'Transducer')

    SimulationParameters.IntegrationDensity = 1.5;
    SimulationParameters.TransducerOnGrid = false;
    Transducer.Axis = 3;
    Transducer = update_transducer(...
        Transducer,Medium,Grid,SimulationParameters);
    
    X = Transducer.integration_points(:,:,1); X = X(:);
    Y = Transducer.integration_points(:,:,2); Y = Y(:);
    Z = Transducer.integration_points(:,:,3); Z = Z(:);
    x0 = 0; y0 = 0; z0 = 0;
    R = eye(3);
    
end

source.points = [X(:) Y(:) Z(:)]*transpose(R);
source.points = source.points + [x0 y0 z0];
source.type = 'pressure';
source.normal = zeros(size(source.points));
source.normal(:,3) = 1;
source.weights = Grid.dx*Grid.dy*ones(size(source.points,1),1);

rmin = min(vecnorm(source.points-sensor.points,2,2));
rmax = max(vecnorm(source.points-sensor.points,2,2));

sensor.offset = rmin;

M = ceil((rmax-rmin)/medium.sound_speed*Fs);

if strcmp(simtype,'pressure') || strcmp(simtype,'velocity')
    N = length(P);
    P = [P zeros(1,M+N)];
    N = length(P);
else
    P = cat(3, P, zeros(size(P,1),size(P,2),M));
    N = size(P,3);
end

M = floor(N/2) + 1;
f = (0:(N-1))/N*Fs;

if strcmp(simtype,'pressure') || strcmp(simtype,'velocity')
    Transmit.PressureSignal = P;
    Transmit.Delays         = zeros(1,Transducer.NumberOfElements);
    Transmit.Apodization    =  ones(1,Transducer.NumberOfElements);
    Transmit.SamplingRate   = Fs;

    source = define_source_frequency(Transducer, Transmit, Medium, ...
        Grid, [], true);
    
    % % Keep only the frequencies with nonzero, independent content
    source.velocities  = source.velocities(:,2:M);
    source.frequencies = source.frequencies(2:M);
else
    % % Keep only the frequencies with nonzero, independent content
    frequencies = f(2:M);
    P = conj(fft(P,[],3));
    p = P(:,:,2:M);

    source.frequencies = frequencies;
    source.pressures = reshape(p,[],length(frequencies));
end

if strcmp(simtype,'pressure')
    source.type = 'pressure';
    source.pressures = source.velocities*medium.density*medium.sound_speed;
    source = rmfield(source,'velocities');
end

medium.direction   = 'forward';
medium.useGPU      = true;

gpuDevice(1)
availableMemory = gpuDevice().AvailableMemory;
memoryEstimate = get_memory_estimate(source,sensor)*3;
Nbatch = ceil(memoryEstimate/availableMemory);
BS = floor(size(source.points,1)/Nbatch);

[source,sensor] = rayleigh_integral_batched(source, sensor, medium, BS);

% Convert frequency domain data to time domain data
P = zeros(1,N);
P(2:M) = sensor.pressures;
P = ifft(conj(P),'symmetric');
t = (0:(N-1))/Fs + rmin/medium.sound_speed;

% Show results
if strcmp(simtype,'data')
    plot(t,P,'Color',[1 1 1]*0.8)
else
    plot(t,P)
end
hold on
