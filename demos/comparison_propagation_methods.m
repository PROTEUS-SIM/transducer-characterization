% This script shows a comparison of propagation methods and source
% presentations as described in Section 6.2.3 of my thesis.
%
%
% Run this script twice, once with sourceType = 'monopole' and once with
% sourceType = 'dipole'.
%
% Run comparison_propagation_methods_plot.m to recreate Fig. 6.2 of my
% thesis.
%
% NOTE: this scripts requires the installation of k-Wave.
% Adjust the variable PATHS.kWave in path_setup_calibration to match the
% location of k-Wave.
% Adjust the variable PATHS.PROTEUS to match the location of the root
% directory of the PROTEUS repository.
%
% NOTE: this script only works correctly with the latest changes to
% PROTEUS, included in the branch transducer-calibration
%
% NOTE: If sufficient RAM is available, change the value of maxMemory in
% update_sensor_fast.m (part of the PROTEUS toolbox) to at least 4 GB.
% Otherwise, update_sensor_fast will revert to update_sensor, which will
% increase the computation time by two orders of magnitude.
%
% This file is part of the transducer-calibration project, licensed under
% the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

%#ok<*UNRCH>

clear; clc; close all

sourceType = 'monopole';

% Add root directory to path and run path setup:
currentDirectory = fileparts(mfilename('fullpath'));
rootDirectory = fileparts(currentDirectory);
resultsDirectory = fullfile(currentDirectory,'data');
addpath(rootDirectory); clear currentDirectory rootDirectory
PATHS = path_setup_calibration();

if strcmp(sourceType,'dipole')
    savename1 = 'sensor_data_dipole_rayleigh.mat';
    savename2 = 'sensor_data_dipole_angular.mat';
    savename3 = 'sensor_data_dipole_k_Wave.mat';
else
    savename1 = 'sensor_data_monopole_rayleigh.mat';
    savename2 = 'sensor_data_monopole_angular.mat';
    savename3 = 'sensor_data_monopole_k_Wave.mat';
end

savename1 = fullfile(resultsDirectory,savename1);
savename2 = fullfile(resultsDirectory,savename2);
savename3 = fullfile(resultsDirectory,savename3);

addpath(PATHS.PROTEUS)
addpath(PATHS.kWave)
load('GUI_output_parameters.mat', 'Geometry', ...
    'Medium', 'SimulationParameters', 'Transducer', 'Transmit')

%==========================================================================
% K-WAVE SIMULATION
%==========================================================================

% Depth of the k-Wave simulation domain
Geometry.Domain.Xmax = 10e-3;

% Do not include the vessel structure in the medium
Geometry.EmbedVessel = false;

% Run the k-Wave C++ code on GPU
SimulationParameters.Solver = '3DG';

% Set up the k-Wave grid
Grid = define_grid(SimulationParameters,Geometry);
kgrid = kWaveGrid(Grid.Nx,Grid.dx,Grid.Ny,Grid.dy,Grid.Nz,Grid.dz);
kgrid.dt = Grid.dt;

% In k-Wave, dipole sources are defined on the points between grid points.
% Shift the grid by dx/2, such that the dipole source plane coincides with
% the transducer surface.
if strcmp(sourceType,'dipole')
    Grid.x = Grid.x - Grid.dx/2;
end

% Position of the sensor [m]
x = 5e-3; y = 0; z = 0;

[~,ix] = min(abs(Grid.x-x));
[~,iy] = min(abs(Grid.y-y));
[~,iz] = min(abs(Grid.z-z));

sensorpoint = [x y z];
sensoridx   = sub2ind([Grid.Nx,Grid.Ny,Grid.Nz],ix,iy,iz);


% Divide the transducer surfaces into integration elements:
SimulationParameters.IntegrationDensity = 1.5;
SimulationParameters.TransducerOnGrid = false;
Transducer=update_transducer(Transducer,Medium,Grid,SimulationParameters);

% Filter and resample transmit signal:
Transmit = preprocess_transmit(Transmit,Medium,kgrid);

% Simulate the field with all elements on:
Transmit.SeqPulse = 'full';

Transducer.SourceType = sourceType;

% Make the signal longer to accommodate maximum travel times
source.points = reshape(Transducer.integration_points,[],3);
sensor.points = sensorpoint;

[~,rmax] = find_min_max_distance(source,sensor);
rmin = 0;
M = ceil((rmax-rmin)/Medium.SpeedOfSound*Transmit.SamplingRate);
Transmit.PressureSignal = [Transmit.PressureSignal, zeros(1,M)];

clear source sensor

% simulation settings
run_param = sim_setup(SimulationParameters);

% Location of the geometry data:
Geometry.GeometriesPath = run_param.GeometriesPath;

run_param.PML = Grid.PML;

% Define the k-Wave medium:
disp('Creating k-Wave medium ...')
medium = define_medium(Grid, Medium, Geometry);

% Record signals long enough for back and forth pass of the wave
run_param = compute_travel_times(run_param, ...
    Geometry,Medium,Transducer,Transmit);

% Create the time array
kgrid.Nt = floor(run_param.tr(1) / kgrid.dt) + 1;

% define the transducer source
disp('Creating k-Wave sensor object for transducer.')
[sensor_transducer, sensor_weights] = define_sensor_transducer(...
    Transducer, Grid);

mask_idx_trans = find(logical(sensor_transducer.mask));

disp('Creating k-Wave source object for transducer.')
source = define_source_transducer(Transducer, Transmit, ...
    Medium, Grid, transpose(sensor_weights), mask_idx_trans);

% Compute the sensor grid points and sensor weigths for off-grid sensors.
Grid.threshold = 24; % Truncation threshold
sensor.mask = zeros(Grid.Nx,Grid.Ny,Grid.Nz,'logical');
[sensor, sensor_weights] = update_sensor_fast(sensor, ...
    sensorpoint, sensoridx, Grid, false);
sensor.record={'p'};  

% Run the k-Wave simulation
sensor_data = run_simulation(run_param, kgrid, medium, source, sensor);

sensor_data.p = full(sensor_weights*double(sensor_data.p));

%==========================================================================
% COMPUTE THE RAYLEIGH INTEGRAL IN FREQUENCY DOMAIN
%==========================================================================

clear source sensor medium

sensor.points = sensorpoint;
sensor.offset = rmin;

useGPU = true;

medium.density     = Medium.Density;
medium.sound_speed = Medium.SpeedOfSound;
medium.direction   = 'forward';
medium.useGPU      = useGPU;

source = define_source_frequency(...
    Transducer,Transmit,Medium,Grid,[],useGPU);

if strcmp(sourceType,'dipole')
    rho0 = medium.density;
    c0   = medium.sound_speed;
    source.type = 'pressure';
    source.pressures = source.velocities*rho0*c0;
    source = rmfield(source,'velocities');
end

% Get batch size for computing the Rayleigh integral:
gpuDevice(1);
availableMemory = gpuDevice().AvailableMemory;
memoryEstimate = get_memory_estimate(source,sensor)*3;
Nbatch = ceil(memoryEstimate/availableMemory);
BS = floor(size(sensor.points,1)/Nbatch);

% Compute the Rayleigh integral
[~,sensor] = rayleigh_integral_batched(source, sensor, medium,BS);

N = length(Transmit.PressureSignal);

p1 = sensor.pressures;
p1 = ifft(conj(p1),N,2,'symmetric');

% Resample the signal at the original sampling rate
Fs = Transducer.SamplingRate;
t1 = (0:(N-1))/Transmit.SamplingRate + rmin/Medium.SpeedOfSound;
t2 = t1(1):1/Fs:t1(end);
p2 = sinc_interpolation(t1,p1,t2);

%==========================================================================
% ANGULAR SPECTRUM METHOD
%==========================================================================

% Change coordinate system:
x = sensorpoint(2);
y = sensorpoint(3);
z = sensorpoint(1);

load('GUI_output_parameters.mat','Transducer')

% Grid for the angular spectrum method
expansion_factor = 4; % Grid expansion factor

Geometry.Axis = 3; % Propagation axis
Grid = define_grid_demos(SimulationParameters,Geometry,expansion_factor);

% Source representation for angular spectrum method
SimulationParameters.IntegrationDensity = 1;
SimulationParameters.TransducerOnGrid = false;
Transducer.Axis = 3; % Propagation axis
Transducer=update_transducer(Transducer,Medium,Grid,SimulationParameters);

% Field distribution representation of transducer:
disp('Computing grid based source presentation ...')
[S, frequencies] = define_source_grid(...
    Transducer, Transmit, Medium, Grid, [], true);

d = z - Grid.z;

% Compute the wavenumber of propagation:
k  = 2*pi*frequencies/medium.sound_speed;

Grid = get_kspace(Grid);
Nx = Grid.Nx; dkx = Grid.dkx;
Ny = Grid.Ny; dky = Grid.dky;

% Forward propagator:
if strcmp(sourceType,'dipole')
    G = compute_spectral_propagator(Nx,dkx,Ny,dky,k,d,true,'pp');
else
    G = compute_spectral_propagator(Nx,dkx,Ny,dky,k,d,true,'vp');
end

% Apply forward propagation:
disp('Computing angular spectrum ...')
S = reshape(S,[size(S,[1 2]),1,size(S,3)]);
S = angular_spectrum(S,Grid.x,Grid.y);
disp('Applying spectral propagator ...')
S = S.*G*medium.density*medium.sound_speed;
disp('Computing inverse angular spectrum ...')
S = angular_spectrum_inverse(S,Grid.x,Grid.y);
S = squeeze(S);

% Revert to time domain:
P = zeros([size(S,[1 2]),length(Transmit.PressureSignal)]);
P(:,:,1:length(frequencies)) = S;
P = ifft(conj(P),[],3,'symmetric');

p3 = squeeze(P(Grid.x==x,Grid.y==y,:));
t3 = (0:(N-1))/Transmit.SamplingRate;
p3 = sinc_interpolation(t3,p3,t2);

p4 = sinc_interpolation(kgrid.t_array + Grid.dt,sensor_data.p,t2);

%==========================================================================
% SHOW RESULTS
%==========================================================================

plot(t2,p2)
hold on
plot(t2,p3)
plot(t2,p4)

if not(isfolder(resultsDirectory))
    mkdir(resultsDirectory)
end

t = t2;
p = p2;
save(savename1,'t','p')
p = p3;
save(savename2,'t','p')
p = p4;
save(savename3,'t','p')

rmpath(PATHS.PROTEUS)
rmpath(PATHS.kWave)
