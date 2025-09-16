function [Y, rmin] = rayleigh_propagation(...
    P,Grid,medium,SourcePlane,SensorPlane,PATHS,options)
% This function sets up source and sensor structures and computes the
% rayleigh integral.
%
% Input:
% - P: frequency-domain pressure data (Nx-by-Ny-by-Nsamples)
%
% Output:
% - Y:    frequency-domain pressure or velocity data (depending on
%         options.sourceType and options.direction) (Nx-by-Ny-by-Nsamples)
% - rmin: offset distance used in Rayleigh integral, see
%         rayleigh_integral_frequency
%
% This file is part of the transducer-characterization project, licensed
% under the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

%==========================================================================
% SET UP SOURCE AND SENSOR STRUCTURES
%==========================================================================

[X,Y,Z] = ndgrid(Grid.x,Grid.y,Grid.z);

if isempty(SourcePlane)
    thetaX = 0; thetaY = 0; thetaZ = 0; x0 = 0; y0 = 0; z0 = 0;
else
    thetaX = SourcePlane.thetaX;
    thetaY = SourcePlane.thetaY;
    thetaZ = SourcePlane.thetaZ;
    x0 = SourcePlane.x0;
    y0 = SourcePlane.y0;
    z0 = SourcePlane.z0;
end

% Get rotation matrices:
Rx = get_rotation_matrix(thetaX/180*pi,1);
Ry = get_rotation_matrix(thetaY/180*pi,2);
Rz = get_rotation_matrix(thetaZ/180*pi,3);
% Rotate about x, y, and z axis respectively:
R = Rz*Ry*Rx;

source.type = options.sourceType;
source.points = [X(:) Y(:) Z(:)]*transpose(R);
source.points = source.points + [x0 y0 z0];
source.normal = zeros(size(source.points));
source.normal(:,3) = 1;
source.normal = source.normal*transpose(R);
source.weights = Grid.dx*Grid.dy*ones(size(source.points,1),1);

if isempty(SensorPlane)
    thetaX = 0; thetaY = 0; thetaZ = 0; x0 = 0; y0 = 0; z0 = 0;
else
    thetaX = SensorPlane.thetaX;
    thetaY = SensorPlane.thetaY;
    thetaZ = SensorPlane.thetaZ;
    x0 = SensorPlane.x0;
    y0 = SensorPlane.y0;
    z0 = SensorPlane.z0;
end

% Get rotation matrices:
Rx = get_rotation_matrix(thetaX/180*pi,1);
Ry = get_rotation_matrix(thetaY/180*pi,2);
Rz = get_rotation_matrix(thetaZ/180*pi,3);
% Rotate about x, y, and z axis respectively:
R = Rz*Ry*Rx;

sensor.points = [X(:) Y(:) Z(:)]*transpose(R);
sensor.points = sensor.points + [x0 y0 z0];
sensor.normal = zeros(size(sensor.points));
sensor.normal(:,3) = -1;
sensor.normal = sensor.normal*transpose(R);
sensor.weights = Grid.dx*Grid.dy*ones(size(sensor.points,1),1);

Fs = 1/Grid.dt;
if strcmp(options.direction,'forward')
    % Adjust the time domain of the signal. Define an offset distance.
    batchMemory = 1024^3; % Memory per batch in bytes
    [rmin,rmax] = find_min_max_distance_gpu(source,sensor,...
        batchMemory,options.useGPU);

    N = size(P,3) + ceil((rmax-rmin)/medium.sound_speed*Fs);
    P = ifft(conj(P),[],3,'symmetric');
    P = conj(fft(P,N,3));

    sensor.offset = rmin;
else
    rmin = 0;
    N = size(P,3);
end
f = (0:(N-1))/N*Fs;

% Real signals have a conjugate symmetric Fourier transform. Only keep
% the first M frequency components, where M is the frequency component at
% or just before the folding frequency Fs/2. The other frequency components
% can be recovered using the symmetry property.
M = floor(N/2) + 1;

% Keep only the frequencies with nonzero, independent content
frequencies = f(2:M);
P = P(:,:,2:M);

source.frequencies = frequencies;

if strcmp(options.direction,'forward')
    source.pressures = reshape(P,[],length(frequencies));
else
    sensor.pressures = reshape(P,[],length(frequencies));
end

medium.direction   = options.direction;
medium.useGPU      = options.useGPU;

%==========================================================================
% COMPUTE RAYLEIGH INTEGRAL
%==========================================================================

if options.useGPU
    gpuDevice(1);
    availableMemory = gpuDevice().AvailableMemory;
else
    availableMemory = 1024^3;
end
memoryEstimate = get_memory_estimate(source,sensor)*3;
Nbatch = ceil(memoryEstimate/availableMemory);
BS = floor(size(source.points,1)/Nbatch);

[source,sensor] = rayleigh_integral_batched(source,sensor,medium,BS);

% Return data in same format as input data:
Y = zeros(Grid.Nx,Grid.Ny,N);
if strcmp(options.direction,'forward')
    Y(:,:,2:M) = reshape(sensor.pressures, Grid.Nx,Grid.Ny,M-1);
elseif strcmp(source.type,'pressure')
    Y(:,:,2:M) = reshape(source.pressures, Grid.Nx,Grid.Ny,M-1);
else
    Y(:,:,2:M) = reshape(source.velocities,Grid.Nx,Grid.Ny,M-1);
end

%==========================================================================
% VISUALIZATION
%==========================================================================

if options.showFigure
    options.plotType = 'abs';
    rayleigh_figure(sensor,source,Grid,medium,SensorPlane,PATHS,options)

    options.plotType = 'angle';
    rayleigh_figure(sensor,source,Grid,medium,SensorPlane,PATHS,options)
end

end
